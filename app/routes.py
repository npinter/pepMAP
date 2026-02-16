from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
import json
import os
import re

import numpy as np
import pandas as pd
from flask import Blueprint, current_app, jsonify, render_template, request

from .parsers import detect_parser
from .services.data import (
    apply_charge_state_mode,
    apply_q_value_cutoff,
    apply_sample_name_cleanup,
    apply_sample_name_cleanup_value,
    get_sample_name_cleanup_preview,
    find_peptide_positions,
    find_uniprot_id_by_accession,
    find_uniprot_id_by_gene_symbol,
    normalize_organism,
    normalize_search_input,
    parse_custom_features_tsv,
    parse_fasta,
    parse_search_inputs,
    parse_search_labels,
)
from .services.plots import (
    build_peptide_summary_table,
    fetch_protein_features,
    plot_features,
    plot_peptides,
    render_tabs,
)
from .services.storage import get_request_session_id, get_store, normalize_session_id

bp = Blueprint('main', __name__)

def _build_peptide_cache_key(
    proteotypic_only,
    charge_state_mode,
    q_value_cutoff,
    sample_name_cleanup,
    sample_name_custom_pattern
):
    payload = {
        'proteotypic_only': bool(proteotypic_only),
        'charge_state_mode': charge_state_mode,
        'q_value_cutoff': float(q_value_cutoff),
        'sample_name_cleanup': sample_name_cleanup or '',
        'sample_name_custom_pattern': sample_name_custom_pattern or ''
    }
    return json.dumps(payload, sort_keys=True, separators=(',', ':'))


def _get_runtime_cache(session_id):
    runtime = current_app.extensions.setdefault('pepmap_runtime_cache', {})
    session_cache = runtime.get(session_id)
    if not isinstance(session_cache, dict):
        session_cache = {}
        runtime[session_id] = session_cache
    return session_cache


def _load_positions_cache(store, session_id):
    session_cache = _get_runtime_cache(session_id)
    cache = session_cache.get('peptide_positions_cache')
    if not isinstance(cache, dict):
        cache = {}
        session_cache['peptide_positions_cache'] = cache
    return cache


def _deserialize_dataframe(data):
    if data is None:
        return None
    try:
        return pd.read_json(StringIO(data), orient='records')
    except ValueError:
        return None


def _serialize_dataframe(df):
    return df.to_json(orient='records')


def _compute_positions_cached(
    store,
    session_id,
    cache_key,
    report_df,
    fasta_df,
    protein_ids,
    proteotypic_only,
    ep_column,
    protein_index=None
):
    cache = _load_positions_cache(store, session_id)
    bucket = cache.get(cache_key)
    if not isinstance(bucket, dict):
        bucket = {}

    positions = {}
    errors = {}
    to_compute = []

    for protein_id in protein_ids:
        cached_df = _deserialize_dataframe(bucket.get(protein_id))
        if cached_df is not None:
            positions[protein_id] = cached_df
        else:
            to_compute.append(protein_id)

    if to_compute:
        max_workers = min(len(to_compute), max(1, min(4, os.cpu_count() or 2)))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(
                    find_peptide_positions,
                    report_df,
                    fasta_df,
                    protein_id,
                    proteotypic_only,
                    ep_column,
                    protein_index
                ): protein_id
                for protein_id in to_compute
            }
            for future in as_completed(future_map):
                protein_id = future_map[future]
                try:
                    positions_df = future.result()
                except Exception as exc:
                    errors[protein_id] = exc
                    continue
                positions[protein_id] = positions_df
                bucket[protein_id] = _serialize_dataframe(positions_df)

        cache[cache_key] = bucket
    return positions, errors


def _build_report_cache_key(
    sample_name_cleanup,
    sample_name_custom_pattern,
    q_value_cutoff,
    charge_state_mode
):
    payload = {
        'sample_name_cleanup': sample_name_cleanup or '',
        'sample_name_custom_pattern': sample_name_custom_pattern or '',
        'q_value_cutoff': float(q_value_cutoff),
        'charge_state_mode': charge_state_mode
    }
    return json.dumps(payload, sort_keys=True, separators=(',', ':'))


def _load_report_cache(store, session_id):
    session_cache = _get_runtime_cache(session_id)
    cache = session_cache.get('report_df_cache')
    if not isinstance(cache, dict):
        cache = {}
        session_cache['report_df_cache'] = cache
    return cache


def _get_report_df_base(store, session_id, report_data):
    session_cache = _get_runtime_cache(session_id)
    base_df = session_cache.get('report_df_base')
    if isinstance(base_df, pd.DataFrame):
        return base_df, True
    base_df = pd.read_json(StringIO(report_data))
    session_cache['report_df_base'] = base_df
    return base_df, False


def _get_protein_index(store, session_id, report_data, protein_ids=None):
    session_cache = _get_runtime_cache(session_id)
    protein_index = session_cache.get('protein_index')
    if not isinstance(protein_index, dict):
        protein_index = {}
        session_cache['protein_index'] = protein_index

    if not protein_ids:
        return protein_index, True

    targets = [pid for pid in protein_ids if pid and pid not in protein_index]
    if not targets:
        return protein_index, True

    base_df, _ = _get_report_df_base(store, session_id, report_data)
    ids_series = base_df.get('Protein.Ids')
    if ids_series is None:
        for pid in targets:
            protein_index.setdefault(pid, [])
        return protein_index, False

    escaped = [re.escape(pid) for pid in targets]
    pattern = r'(^|;)(' + '|'.join(escaped) + r')(;|$)'
    mask = ids_series.fillna('').str.contains(pattern, na=False, regex=True)
    target_set = set(targets)
    matched_rows = 0
    for row_index, raw_value in ids_series[mask].items():
        matched_rows += 1
        for token in str(raw_value).split(';'):
            token = token.strip()
            if token in target_set:
                protein_index.setdefault(token, []).append(row_index)

    for pid in targets:
        protein_index.setdefault(pid, [])
    return protein_index, False


def _filter_report_df_by_protein_ids(report_df, protein_ids, protein_index=None):
    if report_df.empty or not protein_ids:
        return report_df.iloc[0:0]
    if 'Protein.Ids' not in report_df.columns:
        return report_df.iloc[0:0]
    if isinstance(protein_index, dict):
        indices = set()
        for pid in protein_ids:
            indices.update(protein_index.get(pid, []))
        if not indices:
            return report_df.iloc[0:0]
        return report_df.loc[report_df.index.intersection(indices)]

    escaped = [re.escape(pid) for pid in protein_ids if pid]
    if not escaped:
        return report_df.iloc[0:0]
    pattern = r'(^|;)(' + '|'.join(escaped) + r')(;|$)'
    mask = report_df['Protein.Ids'].fillna('').str.contains(pattern, na=False, regex=True)
    return report_df.loc[mask]


def _parse_selected_runs(raw_selected_runs):
    if raw_selected_runs is None:
        return None
    if raw_selected_runs == '':
        return None
    try:
        parsed = json.loads(raw_selected_runs)
    except Exception:
        return None
    if not isinstance(parsed, list):
        return None
    runs = [str(item).strip() for item in parsed if str(item).strip()]
    return runs


def _filter_peptide_positions_by_runs(peptide_positions_df, selected_runs):
    if peptide_positions_df.empty:
        return peptide_positions_df
    if selected_runs is None:
        return peptide_positions_df
    if not selected_runs:
        return peptide_positions_df.iloc[0:0]
    run_set = {str(run) for run in selected_runs}
    return peptide_positions_df.loc[peptide_positions_df['Run'].astype(str).isin(run_set)]


def _compute_global_log2_range(peptide_positions_map, protein_ids, selected_runs):
    if not protein_ids:
        return 0.0, 0.0
    run_set = None
    if selected_runs is not None:
        run_set = {str(run) for run in selected_runs}
        if not run_set:
            return 0.0, 0.0
    values = []
    for protein_id in protein_ids:
        peptide_positions_df = peptide_positions_map.get(protein_id)
        if peptide_positions_df is None or peptide_positions_df.empty:
            continue
        if run_set is not None:
            peptide_positions_df = peptide_positions_df.loc[
                peptide_positions_df['Run'].astype(str).isin(run_set)
            ]
            if peptide_positions_df.empty:
                continue
        log2_values = np.log2(peptide_positions_df['Intensity'])
        log2_values = log2_values[np.isfinite(log2_values)]
        if log2_values.size:
            values.append(np.asarray(log2_values))
    if not values:
        return 0.0, 0.0
    combined = np.concatenate(values)
    return float(combined.min()), float(combined.max())


def _clear_runtime_cache(session_id):
    runtime = current_app.extensions.get('pepmap_runtime_cache')
    if not isinstance(runtime, dict):
        return
    runtime.pop(session_id, None)


def _prefetch_uniprot_features(protein_ids):
    if not protein_ids:
        return
    unique_ids = sorted({pid for pid in protein_ids if pid})
    if len(unique_ids) <= 1:
        return
    max_workers = min(len(unique_ids), max(1, min(6, os.cpu_count() or 2)))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(fetch_protein_features, protein_id) for protein_id in unique_ids]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception:
                continue


def _get_report_df_cached(
    store,
    session_id,
    report_data,
    cache_key,
    sample_name_cleanup,
    sample_name_custom_pattern,
    q_value_cutoff,
    charge_state_mode
):
    cache = _load_report_cache(store, session_id)
    cached = cache.get(cache_key)
    if isinstance(cached, pd.DataFrame):
        return cached, _find_ep_column(cached)

    base_df, _ = _get_report_df_base(store, session_id, report_data)
    report_df = base_df
    ep_column = _find_ep_column(report_df)
    if not ep_column:
        return report_df, None
    report_df = apply_sample_name_cleanup(report_df, sample_name_cleanup, sample_name_custom_pattern)
    report_df = apply_q_value_cutoff(report_df, q_value_cutoff, ep_column)
    report_df = apply_charge_state_mode(report_df, charge_state_mode, ep_column)
    cache[cache_key] = report_df
    return report_df, ep_column


def _find_ep_column(report_df):
    return next((col for col in report_df.columns if col.startswith('EP')), None)


@bp.route('/', methods=['GET'])
def index():
    return render_template('index.html', app_version=current_app.config['APP_VERSION'])


@bp.route('/resolve_labels', methods=['POST'])
def resolve_labels():
    payload = request.get_json(silent=True) or {}
    session_id = normalize_session_id(payload.get('session_id'))
    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400
    identifiers = payload.get('identifiers') or []
    if not isinstance(identifiers, list):
        return jsonify({'error': 'identifiers must be a list.'}), 400

    store = get_store()
    fasta_data = store.read(session_id, 'fasta_data')
    if not fasta_data:
        return jsonify({'error': 'Session does not contain FASTA data.'}), 400
    fasta_df = pd.read_json(StringIO(fasta_data))

    resolved = []
    for raw in identifiers:
        token = normalize_search_input(raw)
        if not token:
            resolved.append({'input': raw, 'uniprot_id': None, 'label': '', 'found': False})
            continue
        uniprot_id = find_uniprot_id_by_accession(fasta_df, token)
        if uniprot_id is None:
            uniprot_id = find_uniprot_id_by_gene_symbol(fasta_df, token)
        if uniprot_id is None:
            resolved.append({'input': token, 'uniprot_id': token, 'label': token, 'found': False})
            continue
        gene_symbol = fasta_df.loc[fasta_df['uniprot_id'] == uniprot_id, 'gene_symbol'].iloc[0]
        label = gene_symbol if isinstance(gene_symbol, str) and gene_symbol.strip() else token
        resolved.append({'input': token, 'uniprot_id': uniprot_id, 'label': label, 'found': True})

    return jsonify({'resolved': resolved}), 200


@bp.route('/plot_peptides', methods=['POST'])
def plot_peptides_route():
    search_input = normalize_search_input(request.form.get('search_input'))
    custom_title = request.form.get('custom_title', '')
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400

    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    store = get_store()
    fasta_data = store.read(session_id, 'fasta_data')
    report_data = store.read(session_id, 'report_data')
    if not fasta_data or not report_data or not search_input:
        return jsonify({'error': 'All fields must be provided (and session must contain uploaded data).'}), 400

    settings = store.read(session_id, 'settings') or {}
    proteotypic_only = request.form.get('proteotypic_checkbox')
    if proteotypic_only is None:
        proteotypic_only = bool(settings.get('proteotypic_only', False))
    else:
        proteotypic_only = proteotypic_only == 'true'
    charge_state_mode = request.form.get('charge_state_mode') or settings.get('charge_state_mode', 'all')
    raw_q_value_cutoff = request.form.get('q_value_cutoff') or str(settings.get('q_value_cutoff', '0.01'))
    sample_name_cleanup = request.form.get('sample_name_cleanup') or settings.get('sample_name_cleanup', 'none')
    sample_name_custom_pattern = request.form.get('sample_name_custom_pattern')
    if sample_name_custom_pattern is None:
        sample_name_custom_pattern = settings.get('sample_name_custom_pattern', '')
    summary_mode = request.form.get('summary_mode')
    if summary_mode is None:
        summary_mode = settings.get('summary_mode', 'per_sample')
    selected_runs = _parse_selected_runs(request.form.get('selected_runs'))

    try:
        q_value_cutoff = float(raw_q_value_cutoff)
    except (TypeError, ValueError):
        q_value_cutoff = 0.01
    q_value_cutoff = max(0.0, min(0.05, q_value_cutoff))
    charge_state_mode = charge_state_mode if charge_state_mode in {'all', 'unique', 'overlap'} else 'all'
    if summary_mode not in {'per_sample', 'per_topology', 'per_domain'}:
        summary_mode = 'per_sample'

    fasta_df = pd.read_json(StringIO(fasta_data))

    search_inputs = parse_search_inputs(search_input)
    if not search_inputs:
        return jsonify({'error': 'No valid search input provided.'}), 400
    search_labels = parse_search_labels(request.form.get('search_labels'), len(search_inputs))
    is_multi = len(search_inputs) > 1

    token_items = []
    for idx, token in enumerate(search_inputs):
        label_override = search_labels[idx] if idx < len(search_labels) else ''
        label_override = label_override.strip()
        fallback_label = label_override if label_override and label_override != token else token

        selected_protein_id = find_uniprot_id_by_gene_symbol(fasta_df, token)
        if selected_protein_id is None:
            selected_protein_id = find_uniprot_id_by_accession(fasta_df, token)

        token_items.append({
            'token': token,
            'label_override': label_override,
            'fallback_label': fallback_label,
            'protein_id': selected_protein_id
        })

    protein_ids = [item['protein_id'] for item in token_items if item['protein_id']]
    report_charge_mode = 'all' if charge_state_mode == 'overlap' else charge_state_mode
    report_cache_key = _build_report_cache_key(
        sample_name_cleanup,
        sample_name_custom_pattern,
        q_value_cutoff,
        report_charge_mode
    )

    report_df, ep_column = _get_report_df_cached(
        store,
        session_id,
        report_data,
        report_cache_key,
        sample_name_cleanup,
        sample_name_custom_pattern,
        q_value_cutoff,
        report_charge_mode
    )
    if not ep_column:
        return jsonify({'error': 'No EP column found in the report.'}), 400

    protein_index = None
    if charge_state_mode == 'overlap':
        protein_index, _ = _get_protein_index(store, session_id, report_data, protein_ids)
        report_df = _filter_report_df_by_protein_ids(report_df, protein_ids, protein_index)
        report_df = apply_charge_state_mode(report_df, 'overlap', ep_column)

    if 'Run' in report_df.columns:
        all_runs = sorted(report_df['Run'].astype(str).unique())
    else:
        all_runs = []
    if selected_runs is not None:
        selected_set = {str(run) for run in selected_runs}
        all_runs = [run for run in all_runs if run in selected_set]
    ep_name = ep_column

    tab_items = []
    counts_by_protein = {}
    feature_counts_by_protein = {}
    feature_labels_by_protein = {}
    table_labels = []

    cache_key = _build_peptide_cache_key(
        proteotypic_only,
        charge_state_mode,
        q_value_cutoff,
        sample_name_cleanup,
        sample_name_custom_pattern
    )
    if protein_index is None:
        protein_index, _ = _get_protein_index(store, session_id, report_data, protein_ids)
    peptide_positions_map, peptide_errors = _compute_positions_cached(
        store,
        session_id,
        cache_key,
        report_df,
        fasta_df,
        protein_ids,
        proteotypic_only,
        ep_column,
        protein_index
    )
    if summary_mode != 'per_sample':
        _prefetch_uniprot_features(protein_ids)

    filtered_positions_map = {
        protein_id: _filter_peptide_positions_by_runs(
            peptide_positions_map.get(protein_id, pd.DataFrame()),
            selected_runs
        )
        for protein_id in protein_ids
    }
    global_log2_min, global_log2_max = _compute_global_log2_range(
        filtered_positions_map,
        protein_ids,
        None
    )

    def build_protein_result(item):
        selected_protein_id = item['protein_id']
        label_override = item['label_override']
        peptide_positions_df = filtered_positions_map.get(selected_protein_id, pd.DataFrame())

        protein_label = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'gene_symbol'].iloc[0]
        display_label = protein_label if isinstance(protein_label, str) and protein_label.strip() else selected_protein_id
        override_label = label_override if label_override else ''
        if override_label:
            display_label = override_label

        plot_title = display_label if is_multi else (override_label or custom_title or display_label)
        plot_html = plot_peptides(
            peptide_positions_df,
            fasta_df,
            selected_protein_id,
            global_log2_min,
            global_log2_max,
            ep_column,
            ep_name,
            plot_title,
            all_runs=all_runs
        )

        if charge_state_mode == 'all':
            unique_peptides = peptide_positions_df.drop_duplicates(subset=['Run', 'Peptide', 'Charge'])
        else:
            unique_peptides = peptide_positions_df.drop_duplicates(subset=['Run', 'Peptide'])
        counts_by_run = unique_peptides.groupby('Run')['Peptide'].size().to_dict()

        feature_labels = None
        feature_counts = None
        if summary_mode != 'per_sample':
            protein_features = fetch_protein_features(selected_protein_id)
            if summary_mode == 'per_domain':
                feature_candidates = [f for f in protein_features if f.get('type') == 'DOMAIN']
                fallback_label = 'Domain'
            else:
                feature_candidates = [
                    f for f in protein_features if f.get('type') in {'TOPO_DOM', 'TRANSMEM'}
                ]
                fallback_label = 'Topology'

            feature_defs = []
            for feature in feature_candidates:
                start = feature.get('start')
                end = feature.get('end')
                if start is None or end is None:
                    continue
                description = (feature.get('description') or '').strip()
                if not description:
                    if feature.get('type') == 'TRANSMEM':
                        description = 'Transmembrane'
                    elif feature.get('type') == 'TOPO_DOM':
                        description = 'Topological domain'
                    else:
                        description = fallback_label
                if start == end:
                    label = f"{description} ({start})"
                else:
                    label = f"{description} ({start}-{end})"
                feature_defs.append({
                    'label': label,
                    'start': int(start),
                    'end': int(end)
                })
            feature_defs.sort(key=lambda item: (item['start'], item['end']))

            feature_labels = [item['label'] for item in feature_defs]
            feature_counts = {}
            for feature in feature_defs:
                label = feature['label']
                overlap = (
                    (peptide_positions_df['Start'] <= feature['end'])
                    & (peptide_positions_df['End'] >= feature['start'])
                )
                if charge_state_mode == 'all':
                    overlap_df = peptide_positions_df.loc[
                        overlap, ['Run', 'Peptide', 'Charge']
                    ].drop_duplicates()
                else:
                    overlap_df = peptide_positions_df.loc[
                        overlap, ['Run', 'Peptide']
                    ].drop_duplicates()
                feature_counts[label] = overlap_df.groupby('Run').size().to_dict()

        result = {
            'display_label': display_label,
            'plot_html': plot_html,
            'counts_by_run': counts_by_run,
            'feature_labels': feature_labels,
            'feature_counts': feature_counts
        }
        return result

    results_by_index = {}
    max_workers = min(len(token_items), max(1, min(4, os.cpu_count() or 2)))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {}
        for index, item in enumerate(token_items):
            selected_protein_id = item['protein_id']
            if selected_protein_id is None:
                continue
            if selected_protein_id in peptide_errors:
                continue
            peptide_positions_df = filtered_positions_map.get(selected_protein_id, pd.DataFrame())
            if peptide_positions_df.empty:
                continue
            future = executor.submit(build_protein_result, item)
            future_map[future] = index

        for future in as_completed(future_map):
            index = future_map[future]
            try:
                results_by_index[index] = future.result()
            except Exception as exc:
                results_by_index[index] = exc

    for index, item in enumerate(token_items):
        fallback_label = item['fallback_label']
        label_override = item['label_override']
        selected_protein_id = item['protein_id']

        if selected_protein_id is None:
            if not is_multi:
                return jsonify({'error': 'No protein found for the given search input.'}), 400
            tab_items.append({
                'label': fallback_label,
                'content': f'<div class="text-red-600">No protein found for {fallback_label}.</div>'
            })
            table_labels.append(fallback_label)
            counts_by_protein.setdefault(fallback_label, {})
            continue

        if selected_protein_id in peptide_errors:
            error = peptide_errors[selected_protein_id]
            if not is_multi:
                return jsonify({'error': str(error)}), 400
            tab_items.append({
                'label': fallback_label,
                'content': f'<div class="text-red-600">{str(error)}</div>'
            })
            continue

        peptide_positions_df = filtered_positions_map.get(selected_protein_id, pd.DataFrame())
        if peptide_positions_df.empty:
            message = f'No peptide positions found for {label_override or selected_protein_id}.'
            if not is_multi:
                return jsonify({'error': message}), 400
            tab_items.append({
                'label': fallback_label,
                'content': f'<div class="text-red-600">{message}</div>'
            })
            continue

        result = results_by_index.get(index)
        if isinstance(result, Exception) or result is None:
            if not is_multi:
                return jsonify({'error': str(result) if result else 'Failed to build plot.'}), 400
            tab_items.append({
                'label': fallback_label,
                'content': f'<div class="text-red-600">{str(result) if result else "Failed to build plot."}</div>'
            })
            continue

        display_label = result['display_label']
        tab_items.append({
            'label': display_label,
            'content': result['plot_html']
        })

        table_label = display_label
        table_labels.append(table_label)
        counts_by_protein.setdefault(table_label, {})
        counts_by_protein[table_label] = result['counts_by_run']

        if result['feature_labels'] is not None:
            feature_labels_by_protein[table_label] = result['feature_labels']
        if result['feature_counts'] is not None:
            feature_counts_by_protein[table_label] = result['feature_counts']

    summary_html = build_peptide_summary_table(
        summary_mode,
        counts_by_protein,
        all_runs,
        column_order=table_labels,
        feature_counts_by_protein=feature_counts_by_protein,
        feature_labels_by_protein=feature_labels_by_protein
    )
    if is_multi:
        tabs_html = render_tabs('pepmap-peptides-tabs', tab_items)
        return f'{tabs_html}{summary_html}'

    return f'{tab_items[0]["content"]}{summary_html}'


@bp.route('/plot_features', methods=['POST'])
def plot_features_route():
    search_input = normalize_search_input(request.form.get('search_input'))
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400

    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    store = get_store()
    fasta_data = store.read(session_id, 'fasta_data')
    report_data = store.read(session_id, 'report_data')
    if not fasta_data or not report_data or not search_input:
        return jsonify({'error': 'All fields must be provided (and session must contain uploaded data).'}), 400
    fasta_df = pd.read_json(StringIO(fasta_data))
    settings = store.read(session_id, 'settings') or {}
    proteotypic_only = request.form.get('proteotypic_checkbox')
    if proteotypic_only is None:
        proteotypic_only = bool(settings.get('proteotypic_only', False))
    else:
        proteotypic_only = proteotypic_only == 'true'
    charge_state_mode = request.form.get('charge_state_mode') or settings.get('charge_state_mode', 'all')
    raw_q_value_cutoff = request.form.get('q_value_cutoff') or str(settings.get('q_value_cutoff', '0.01'))
    sample_name_cleanup = request.form.get('sample_name_cleanup') or settings.get('sample_name_cleanup', 'none')
    sample_name_custom_pattern = request.form.get('sample_name_custom_pattern')
    if sample_name_custom_pattern is None:
        sample_name_custom_pattern = settings.get('sample_name_custom_pattern', '')
    selected_runs = _parse_selected_runs(request.form.get('selected_runs'))
    try:
        q_value_cutoff = float(raw_q_value_cutoff)
    except (TypeError, ValueError):
        q_value_cutoff = 0.01
    q_value_cutoff = max(0.0, min(0.05, q_value_cutoff))
    charge_state_mode = charge_state_mode if charge_state_mode in {'all', 'unique', 'overlap'} else 'all'

    custom_features_df = None
    custom_features_label = None
    custom_features_data = store.read(session_id, 'custom_features_data')
    if custom_features_data:
        custom_features_df = pd.read_json(StringIO(custom_features_data))
        custom_features_label = store.read(session_id, 'custom_features_label') or 'Custom Features'

    search_inputs = parse_search_inputs(search_input)
    if not search_inputs:
        return jsonify({'error': 'No valid search input provided.'}), 400
    is_multi = len(search_inputs) > 1

    protein_ids = []
    for token in search_inputs:
        selected_protein_id = find_uniprot_id_by_gene_symbol(fasta_df, token)
        if selected_protein_id is None:
            selected_protein_id = find_uniprot_id_by_accession(fasta_df, token)
        if selected_protein_id:
            protein_ids.append(selected_protein_id)

    report_charge_mode = 'all' if charge_state_mode == 'overlap' else charge_state_mode
    report_cache_key = _build_report_cache_key(
        sample_name_cleanup,
        sample_name_custom_pattern,
        q_value_cutoff,
        report_charge_mode
    )
    report_df, ep_column = _get_report_df_cached(
        store,
        session_id,
        report_data,
        report_cache_key,
        sample_name_cleanup,
        sample_name_custom_pattern,
        q_value_cutoff,
        report_charge_mode
    )
    if not ep_column:
        return jsonify({'error': 'No EP column found in the report.'}), 400

    protein_index = None
    if charge_state_mode == 'overlap':
        protein_index, _ = _get_protein_index(store, session_id, report_data, protein_ids)
        report_df = _filter_report_df_by_protein_ids(report_df, protein_ids, protein_index)
        report_df = apply_charge_state_mode(report_df, 'overlap', ep_column)

    cache_key = _build_peptide_cache_key(
        proteotypic_only,
        charge_state_mode,
        q_value_cutoff,
        sample_name_cleanup,
        sample_name_custom_pattern
    )
    if protein_index is None:
        protein_index, _ = _get_protein_index(store, session_id, report_data, protein_ids)
    peptide_positions_map, peptide_errors = _compute_positions_cached(
        store,
        session_id,
        cache_key,
        report_df,
        fasta_df,
        protein_ids,
        proteotypic_only,
        ep_column,
        protein_index
    )
    _prefetch_uniprot_features(protein_ids)

    filtered_positions_map = {
        protein_id: _filter_peptide_positions_by_runs(
            peptide_positions_map.get(protein_id, pd.DataFrame()),
            selected_runs
        )
        for protein_id in protein_ids
    }

    feature_panels = []
    for idx, token in enumerate(search_inputs):
        selected_protein_id = find_uniprot_id_by_gene_symbol(fasta_df, token)
        if selected_protein_id is None:
            selected_protein_id = find_uniprot_id_by_accession(fasta_df, token)
        if selected_protein_id is None:
            if not is_multi:
                return jsonify({'error': 'No protein found for the given search input.'}), 400
            feature_panels.append(
                f'<div class="feature-panel" data-index="{idx}">'
                f'<div class="text-red-600">No protein found for {token}.</div>'
                '</div>'
            )
            continue

        if selected_protein_id in peptide_errors:
            peptide_positions_df = pd.DataFrame()
        else:
            peptide_positions_df = filtered_positions_map.get(selected_protein_id, pd.DataFrame())
        if peptide_positions_df.empty:
            if is_multi:
                feature_panels.append(
                    f'<div class="feature-panel empty-feature-panel" data-index="{idx}"></div>'
                )
                continue
            return ''

        try:
            features_html = plot_features(fasta_df, selected_protein_id, custom_features_df, custom_features_label)
        except Exception as exc:
            if not is_multi:
                return jsonify({'error': str(exc)}), 400
            feature_panels.append(
                f'<div class="feature-panel" data-index="{idx}">'
                f'<div class="text-red-600">{str(exc)}</div>'
                '</div>'
            )
            continue

        feature_panels.append(
            f'<div class="feature-panel" data-index="{idx}">{features_html}</div>'
        )

    if is_multi:
        return ''.join(feature_panels)

    return feature_panels[0]




@bp.route('/summary_counts', methods=['POST'])
def summary_counts_route():
    search_input = normalize_search_input(request.form.get('search_input'))
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400
    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    store = get_store()
    fasta_data = store.read(session_id, 'fasta_data')
    report_data = store.read(session_id, 'report_data')
    if not fasta_data or not report_data or not search_input:
        return jsonify({'error': 'All fields must be provided (and session must contain uploaded data).'}), 400

    settings = store.read(session_id, 'settings') or {}
    proteotypic_only = request.form.get('proteotypic_checkbox')
    if proteotypic_only is None:
        proteotypic_only = bool(settings.get('proteotypic_only', False))
    else:
        proteotypic_only = proteotypic_only == 'true'
    charge_state_mode = request.form.get('charge_state_mode') or settings.get('charge_state_mode', 'all')
    raw_q_value_cutoff = request.form.get('q_value_cutoff') or str(settings.get('q_value_cutoff', '0.01'))
    sample_name_cleanup = request.form.get('sample_name_cleanup') or settings.get('sample_name_cleanup', 'none')
    sample_name_custom_pattern = request.form.get('sample_name_custom_pattern')
    if sample_name_custom_pattern is None:
        sample_name_custom_pattern = settings.get('sample_name_custom_pattern', '')
    selected_runs = _parse_selected_runs(request.form.get('selected_runs'))
    summary_mode = request.form.get('summary_mode')
    if summary_mode is None:
        summary_mode = settings.get('summary_mode', 'per_sample')

    try:
        q_value_cutoff = float(raw_q_value_cutoff)
    except (TypeError, ValueError):
        q_value_cutoff = 0.01
    q_value_cutoff = max(0.0, min(0.05, q_value_cutoff))
    charge_state_mode = charge_state_mode if charge_state_mode in {'all', 'unique', 'overlap'} else 'all'
    if summary_mode not in {'per_sample', 'per_topology', 'per_domain'}:
        summary_mode = 'per_sample'

    fasta_df = pd.read_json(StringIO(fasta_data))

    search_inputs = parse_search_inputs(search_input)
    if not search_inputs:
        return jsonify({'error': 'No valid search input provided.'}), 400

    protein_ids = []
    for token in search_inputs:
        selected_protein_id = find_uniprot_id_by_gene_symbol(fasta_df, token)
        if selected_protein_id is None:
            selected_protein_id = find_uniprot_id_by_accession(fasta_df, token)
        if selected_protein_id:
            protein_ids.append(selected_protein_id)

    report_charge_mode = 'all' if charge_state_mode == 'overlap' else charge_state_mode
    report_cache_key = _build_report_cache_key(
        sample_name_cleanup,
        sample_name_custom_pattern,
        q_value_cutoff,
        report_charge_mode
    )

    report_df, ep_column = _get_report_df_cached(
        store,
        session_id,
        report_data,
        report_cache_key,
        sample_name_cleanup,
        sample_name_custom_pattern,
        q_value_cutoff,
        report_charge_mode
    )
    if not ep_column:
        return jsonify({'error': 'No EP column found in the report.'}), 400

    protein_index = None
    if charge_state_mode == 'overlap':
        protein_index, _ = _get_protein_index(store, session_id, report_data, protein_ids)
        report_df = _filter_report_df_by_protein_ids(report_df, protein_ids, protein_index)
        report_df = apply_charge_state_mode(report_df, 'overlap', ep_column)

    if 'Run' in report_df.columns:
        all_runs = sorted(report_df['Run'].astype(str).unique())
    else:
        all_runs = []
    if selected_runs is not None:
        selected_set = {str(run) for run in selected_runs}
        all_runs = [run for run in all_runs if run in selected_set]

    cache_key = _build_peptide_cache_key(
        proteotypic_only,
        charge_state_mode,
        q_value_cutoff,
        sample_name_cleanup,
        sample_name_custom_pattern
    )
    if protein_index is None:
        protein_index, _ = _get_protein_index(store, session_id, report_data, protein_ids)
    peptide_positions_map, peptide_errors = _compute_positions_cached(
        store,
        session_id,
        cache_key,
        report_df,
        fasta_df,
        protein_ids,
        proteotypic_only,
        ep_column,
        protein_index
    )

    filtered_positions_map = {
        protein_id: _filter_peptide_positions_by_runs(
            peptide_positions_map.get(protein_id, pd.DataFrame()),
            selected_runs
        )
        for protein_id in protein_ids
    }

    results = []
    for token in search_inputs:
        selected_protein_id = find_uniprot_id_by_gene_symbol(fasta_df, token)
        if selected_protein_id is None:
            selected_protein_id = find_uniprot_id_by_accession(fasta_df, token)

        entry = {
            'input': token,
            'uniprot_id': selected_protein_id,
            'found': bool(selected_protein_id),
            'label': token,
            'counts_by_run': {},
            'features': []
        }

        if selected_protein_id is None:
            results.append(entry)
            continue

        try:
            protein_label = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'gene_symbol'].iloc[0]
            if isinstance(protein_label, str) and protein_label.strip():
                entry['label'] = protein_label
        except Exception:
            pass

        if selected_protein_id in peptide_errors:
            entry['error'] = str(peptide_errors[selected_protein_id])
            results.append(entry)
            continue

        peptide_positions_df = filtered_positions_map.get(selected_protein_id, pd.DataFrame())

        if not peptide_positions_df.empty:
            if charge_state_mode == 'all':
                unique_peptides = peptide_positions_df.drop_duplicates(subset=['Run', 'Peptide', 'Charge'])
            else:
                unique_peptides = peptide_positions_df.drop_duplicates(subset=['Run', 'Peptide'])
            entry['counts_by_run'] = unique_peptides.groupby('Run')['Peptide'].size().to_dict()

            if summary_mode != 'per_sample':
                protein_features = fetch_protein_features(selected_protein_id)
                if summary_mode == 'per_domain':
                    feature_candidates = [f for f in protein_features if f.get('type') == 'DOMAIN']
                    fallback_label = 'Domain'
                else:
                    feature_candidates = [
                        f for f in protein_features if f.get('type') in {'TOPO_DOM', 'TRANSMEM'}
                    ]
                    fallback_label = 'Topology'

                feature_defs = []
                for feature in feature_candidates:
                    start = feature.get('start')
                    end = feature.get('end')
                    if start is None or end is None:
                        continue
                    description = (feature.get('description') or '').strip()
                    if not description:
                        if feature.get('type') == 'TRANSMEM':
                            description = 'Transmembrane'
                        elif feature.get('type') == 'TOPO_DOM':
                            description = 'Topological domain'
                        else:
                            description = fallback_label
                    if start == end:
                        label = f"{description} ({start})"
                    else:
                        label = f"{description} ({start}-{end})"
                    feature_defs.append({
                        'label': label,
                        'start': int(start),
                        'end': int(end),
                        'type': feature.get('type'),
                        'description': description
                    })
                feature_defs.sort(key=lambda item: (item['start'], item['end']))

                for feature in feature_defs:
                    overlap = (
                        (peptide_positions_df['Start'] <= feature['end'])
                        & (peptide_positions_df['End'] >= feature['start'])
                    )
                    if charge_state_mode == 'all':
                        overlap_df = peptide_positions_df.loc[
                            overlap, ['Run', 'Peptide', 'Charge']
                        ].drop_duplicates()
                    else:
                        overlap_df = peptide_positions_df.loc[
                            overlap, ['Run', 'Peptide']
                        ].drop_duplicates()
                    counts = overlap_df.groupby('Run').size().to_dict()
                    entry['features'].append({
                        'label': feature['label'],
                        'type': feature.get('type'),
                        'description': feature.get('description'),
                        'start': feature['start'],
                        'end': feature['end'],
                        'counts_by_run': counts
                    })

        results.append(entry)

    return jsonify({'summary_mode': summary_mode, 'runs': all_runs, 'proteins': results}), 200

@bp.route('/upload', methods=['POST'])
def upload_files():
    report_file = request.files.get('report_file')
    fasta_file = request.files.get('fasta_file')
    custom_organism = request.form.get('custom_organism')
    organism = normalize_organism(request.form.get('organism'), custom_organism)
    custom_features_file = request.files.get('custom_features_file')
    custom_features_label = normalize_search_input(request.form.get('custom_features_label'))
    settings = {
        'proteotypic_only': request.form.get('proteotypic_checkbox') == 'true',
        'charge_state_mode': request.form.get('charge_state_mode', 'all'),
        'q_value_cutoff': request.form.get('q_value_cutoff', '0.01'),
        'sample_name_cleanup': request.form.get('sample_name_cleanup', 'none'),
        'sample_name_custom_pattern': request.form.get('sample_name_custom_pattern', ''),
        'summary_mode': request.form.get('summary_mode', 'per_sample')
    }
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400

    store = get_store()
    if report_file and fasta_file:
        custom_features_uploaded = False
        parser_name, parser = detect_parser(report_file.stream, filename=report_file.filename)
        report_df, report_logo = parser(report_file.stream, filename=report_file.filename)
        run_names = []
        if 'Run' in report_df.columns:
            run_names = sorted(report_df['Run'].dropna().astype(str).unique())
        if not session_id:
            session_id = store.create(meta={'user_agent': request.headers.get('User-Agent', '')})
        elif store.get_payload(session_id) is None:
            store.create_with_id(session_id, meta={'user_agent': request.headers.get('User-Agent', '')})

        fasta_df = parse_fasta(fasta_file.stream, organism or '')
        gene_symbols = fasta_df['gene_symbol'].astype(str)
        uniprot_ids = fasta_df['uniprot_id'].astype(str)
        candidates = pd.concat([gene_symbols, uniprot_ids], ignore_index=True).dropna()
        candidates = pd.Series(candidates.unique())
        autocomplete_candidates = candidates.tolist()
        autocomplete_candidates_lower = [candidate.lower() for candidate in autocomplete_candidates]

        update_kwargs = {
            'fasta_data': fasta_df.to_json(),
            'report_data': report_df.to_json(),
            'organism': organism,
            'report_filename': report_file.filename,
            'fasta_filename': fasta_file.filename,
            'custom_features_filename': custom_features_file.filename if custom_features_file else None,
            'report_type': parser_name,
            'report_logo': report_logo,
            'run_names': run_names,
            'autocomplete_candidates': autocomplete_candidates,
            'autocomplete_candidates_lower': autocomplete_candidates_lower,
            'settings': settings
        }
        if custom_features_file:
            if not custom_features_label:
                return jsonify({'error': 'Custom features label is required when uploading a custom features file.'}), 400
            update_kwargs['custom_features_data'] = parse_custom_features_tsv(custom_features_file.stream).to_json()
            update_kwargs['custom_features_label'] = custom_features_label
            custom_features_uploaded = True
        else:
            update_kwargs['custom_features_data'] = None
            update_kwargs['custom_features_label'] = None

        store.update(session_id, **update_kwargs)
        _clear_runtime_cache(session_id)
        return jsonify({
            'message': 'Files uploaded successfully',
            'custom_features_uploaded': custom_features_uploaded,
            'custom_features_label': custom_features_label,
            'report_type': parser_name,
            'report_logo': report_logo,
            'runs': run_names,
            'settings': settings,
            'session_id': session_id
        }), 200

    if custom_features_file:
        if not session_id:
            return jsonify({'error': 'session_id is required to upload custom features.'}), 400
        if store.get_payload(session_id) is None:
            return jsonify({'error': 'Unknown or expired session_id.'}), 404
        if not custom_features_label:
            return jsonify({'error': 'Custom features label is required when uploading a custom features file.'}), 400

        store.update(
            session_id,
            custom_features_data=parse_custom_features_tsv(custom_features_file.stream).to_json(),
            custom_features_label=custom_features_label,
            custom_features_filename=custom_features_file.filename
        )
        _clear_runtime_cache(session_id)
        return jsonify({
            'message': 'Custom features uploaded successfully',
            'custom_features_uploaded': True,
            'custom_features_label': custom_features_label,
            'session_id': session_id
        }), 200

    return jsonify({'error': 'Missing files'}), 400


@bp.route('/sample_cleanup_preview', methods=['POST'])
def sample_cleanup_preview():
    payload = request.get_json(silent=True) or {}
    session_id = normalize_session_id(payload.get('session_id'))
    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    cleanup_mode = payload.get('cleanup_mode') or 'none'
    custom_pattern = payload.get('custom_pattern') or ''

    store = get_store()
    run_names = store.read(session_id, 'run_names') or []
    run_names = [str(run) for run in run_names]

    if not run_names:
        return jsonify({'preview': '', 'cleaned_runs': [], 'error': ''}), 200

    preview_value, error = get_sample_name_cleanup_preview(
        run_names[0],
        cleanup_mode,
        custom_pattern
    )

    cleaned_runs = []
    seen = set()
    for run_name in run_names:
        cleaned = apply_sample_name_cleanup_value(run_name, cleanup_mode, custom_pattern)
        if cleaned in seen:
            continue
        seen.add(cleaned)
        cleaned_runs.append(cleaned)

    return jsonify({
        'preview': {
            'original': run_names[0],
            'cleaned': preview_value
        },
        'cleaned_runs': cleaned_runs,
        'error': error or ''
    }), 200


@bp.route('/autocomplete', methods=['GET'])
def autocomplete():
    query = normalize_search_input(request.args.get('query'))
    session_id, session_error = get_request_session_id()
    if session_error or not query:
        return jsonify({'suggestions': []}), 200
    if not session_id:
        return jsonify({'suggestions': []}), 200

    store = get_store()
    candidates = store.read(session_id, 'autocomplete_candidates')
    candidates_lower = store.read(session_id, 'autocomplete_candidates_lower')
    if candidates and candidates_lower:
        query_lower = query.lower()
        matches = [
            candidate for candidate, lower_value in zip(candidates, candidates_lower)
            if query_lower in lower_value
        ]
        suggestions = matches[:50]
        return jsonify({'suggestions': suggestions}), 200

    fasta_data = store.read(session_id, 'fasta_data')
    if not fasta_data:
        return jsonify({'suggestions': []}), 200
    fasta_df = pd.read_json(StringIO(fasta_data))

    gene_symbols = fasta_df['gene_symbol'].astype(str)
    uniprot_ids = fasta_df['uniprot_id'].astype(str)
    candidates = pd.concat([gene_symbols, uniprot_ids], ignore_index=True).dropna()
    candidates = pd.Series(candidates.unique())
    matches = candidates[candidates.str.contains(query, case=False, regex=False)]
    suggestions = matches.head(50).tolist()
    return jsonify({'suggestions': suggestions}), 200


@bp.route('/session_info', methods=['GET'])
def session_info():
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400
    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    store = get_store()
    payload = store.get_payload(session_id)
    if payload is None:
        return jsonify({'error': 'Unknown or expired session_id.'}), 404

    data = payload.get('data', {})
    settings = data.get('settings') or {}
    run_names = data.get('run_names')
    if run_names is None:
        report_data = data.get('report_data')
        if report_data:
            try:
                report_df = pd.read_json(StringIO(report_data))
                if 'Run' in report_df.columns:
                    run_names = sorted(report_df['Run'].dropna().astype(str).unique())
            except Exception:
                run_names = []
    return jsonify({
        'session_id': session_id,
        'report_filename': data.get('report_filename'),
        'fasta_filename': data.get('fasta_filename'),
        'custom_features_filename': data.get('custom_features_filename'),
        'custom_features_label': data.get('custom_features_label'),
        'report_type': data.get('report_type'),
        'report_logo': data.get('report_logo'),
        'settings': settings,
        'runs': run_names or [],
        'organism': data.get('organism'),
        'has_custom_features': bool(data.get('custom_features_data'))
    }), 200


@bp.route('/flush', methods=['POST'])
def flush_session():
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400
    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    store = get_store()
    store.delete(session_id)
    _clear_runtime_cache(session_id)
    return jsonify({'message': 'Session cleared'}), 200
