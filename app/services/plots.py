import html
from functools import lru_cache
from collections import OrderedDict
import numpy as np
import plotly.graph_objs as go
import plotly.io as pio
import requests

_session = requests.Session()


def generate_dynamic_ticks(protein_length):
    ranges = [
        (0, 500, 10),
        (501, 1000, 20),
        (1001, 2000, 50),
        (2001, float('inf'), 100)
    ]

    for start, end, step in ranges:
        if start <= protein_length < end:
            chosen_step = step
            break
    else:
        chosen_step = 100

    tickvals = list(range(0, protein_length + 1, chosen_step))
    ticktext = [str(tick) for tick in tickvals]
    tickvals[0] = 1
    ticktext[0] = '1'

    if protein_length % chosen_step != 0:
        tickvals.append(protein_length)
        ticktext.append('')

    return tickvals, ticktext


def _fetch_protein_features_impl(uniprot_id):
    url = (
        'https://www.ebi.ac.uk/proteins/api/proteins'
        f'?offset=0&size=100&accession={uniprot_id}'
    )
    response = _session.get(url, headers={'Accept': 'application/json'})
    if response.status_code != 200:
        return []

    def normalize_position(value):
        if value is None:
            return None
        if isinstance(value, (int, float)):
            return int(value)
        text = str(value).strip()
        if not text:
            return None
        if text.startswith(('>', '<')):
            text = text[1:]
        if text.isdigit():
            return int(text)
        digits = ''.join(ch for ch in text if ch.isdigit())
        return int(digits) if digits else None

    data = response.json()
    features = []
    for entry in data:
        if 'features' not in entry:
            continue
        for feature in entry['features']:
            feature_type = feature.get('type')
            start = feature.get('begin')
            end = feature.get('end')
            start_pos = normalize_position(start)
            end_pos = normalize_position(end)
            if start_pos is None or end_pos is None:
                continue
            if feature_type == 'DOMAIN':
                features.append({
                    'group': 'Domains',
                    'type': feature_type,
                    'description': feature.get('description') or 'Domain',
                    'start': start_pos,
                    'end': end_pos
                })
            elif feature_type == 'BINDING':
                features.append({
                    'group': 'Binding Sites',
                    'type': feature_type,
                    'description': feature.get('description') or 'Binding',
                    'molecule': feature.get('molecule'),
                    'ligand': feature.get('ligand', {}).get('name'),
                    'start': start_pos,
                    'end': end_pos
                })
            elif feature_type == 'SITE':
                features.append({
                    'group': 'Sites',
                    'type': feature_type,
                    'description': feature.get('description') or 'Site',
                    'start': start_pos,
                    'end': end_pos
                })
            elif feature_type == 'MOD_RES':
                features.append({
                    'group': 'Modified Residues',
                    'type': feature_type,
                    'description': feature.get('description') or 'Modified residue',
                    'start': start_pos,
                    'end': end_pos
                })
            elif feature_type == 'VARIANT':
                features.append({
                    'group': 'Variants',
                    'type': feature_type,
                    'ftID': feature.get('ftId'),
                    'description': feature.get('description') or 'Variant',
                    'alternativeSequence': feature.get('alternativeSequence', ''),
                    'start': start_pos,
                    'end': end_pos
                })
            elif feature_type == 'TOPO_DOM':
                features.append({
                    'group': 'Topology',
                    'type': feature_type,
                    'description': feature.get('description') or 'Topological domain',
                    'start': start_pos,
                    'end': end_pos
                })
            elif feature_type == 'TRANSMEM':
                features.append({
                    'group': 'Transmembrane',
                    'type': feature_type,
                    'description': feature.get('description') or 'Transmembrane',
                    'start': start_pos,
                    'end': end_pos
                })
    return features


@lru_cache(maxsize=256)
def _fetch_protein_features_cached(uniprot_id):
    return tuple(_fetch_protein_features_impl(uniprot_id))


def fetch_protein_features(uniprot_id):
    return [dict(feature) for feature in _fetch_protein_features_cached(uniprot_id)]


def plot_peptides(
    peptide_positions_df,
    fasta_df,
    selected_protein_id,
    global_log2_min,
    global_log2_max,
    p_value_column,
    p_value_name,
    custom_title,
    all_runs=None
):
    protein_sequence = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'sequence'].iloc[0]
    protein_length = len(protein_sequence)

    peptide_positions_df['log2_intensity'] = np.log2(peptide_positions_df['Intensity'])
    intensity_span = global_log2_max - global_log2_min
    if intensity_span == 0:
        intensity_span = 1
    peptide_positions_df['normalized_intensity'] = (
        (peptide_positions_df['log2_intensity'] - global_log2_min) / intensity_span
    )

    if all_runs is None:
        unique_runs = sorted(peptide_positions_df['Run'].unique())
    else:
        unique_runs = list(all_runs)

    peptide_bar_height = 20
    peptide_bar_margin = 5
    peptide_bar_line_width = 1
    run_label_size = 14
    traces = []
    shapes = []
    min_height = 130
    current_y = 0
    y_tickvals = []
    y_ticktext = []

    min_run_height = peptide_bar_height + peptide_bar_margin

    for run in unique_runs:
        group = peptide_positions_df[peptide_positions_df['Run'] == run].copy()
        group.sort_values('Start', inplace=True)

        end_positions_per_y_offset = []
        y_offsets = []

        for _, row in group.iterrows():
            placed = False
            for y_offset_idx, end_time in enumerate(end_positions_per_y_offset):
                if row['Start'] >= end_time:
                    end_positions_per_y_offset[y_offset_idx] = row['End']
                    y_offsets.append(y_offset_idx)
                    placed = True
                    break
            if not placed:
                end_positions_per_y_offset.append(row['End'])
                y_offsets.append(len(end_positions_per_y_offset) - 1)

        group['overlap_offset'] = [offset * (peptide_bar_height + peptide_bar_margin) for offset in y_offsets]
        max_y_offset = max(y_offsets) if y_offsets else 0

        y_base = current_y
        group['y_base'] = y_base

        run_height = max((max_y_offset + 1) * (peptide_bar_height + peptide_bar_margin), min_run_height)
        y_tickvals.append(y_base + run_height / 2)
        y_ticktext.append(run)

        shapes.append({
            'type': 'rect',
            'x0': 0,
            'y0': y_base,
            'x1': protein_length,
            'y1': y_base + run_height,
            'fillcolor': 'rgba(255,255,255,0)',
            'line': {'width': 0},
            'layer': 'below'
        })

        for _, row in group.iterrows():
            color = f"rgba(255,{255 - row['normalized_intensity'] * 255},0,0.8)"
            trace = go.Bar(
                x=[row['End'] - row['Start']],
                y=[row['y_base'] + row['overlap_offset']],
                width=peptide_bar_height,
                base=row['Start'],
                orientation='h',
                marker=dict(
                    color=color,
                    line=dict(color='black', width=peptide_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=(
                    f"<b>{row['Peptide']}</b>"
                    f"<br>Position: {row['Start']}-{row['End']}"
                    f"<br>Log2 Intensity: {row['log2_intensity']:.2f}"
                    f"<br>Charge: {row['Charge']}"
                    f"<br>{p_value_name}: {row[p_value_column]}"
                    f"<br>Proteotypic: {'Yes' if row['Proteotypic'] else 'No'}"
                ),
                showlegend=False
            )
            traces.append(trace)

        shapes.append({
            'type': 'line',
            'x0': 0,
            'y0': y_base + run_height,
            'x1': protein_length,
            'y1': y_base + run_height,
            'line': {
                'color': 'black',
                'width': 1
            }
        })

        current_y += run_height + peptide_bar_margin * 4

    tickvals, ticktext = generate_dynamic_ticks(protein_length)

    selected_protein_name = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'gene_symbol'].iloc[0]
    final_height = max(current_y + 50, min_height)
    display_title = custom_title.strip() if custom_title else selected_protein_name

    layout = go.Layout(
        title=f'Peptide Mapping for {display_title}',
        xaxis=dict(
            range=[1, protein_length],
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext,
            tickangle=0,
            tickfont=dict(size=9),
            ticks='outside',
            fixedrange=True,
            automargin=True
        ),
        yaxis=dict(
            tickmode='array',
            tickvals=y_tickvals,
            ticktext=y_ticktext,
            tickfont=dict(size=run_label_size),
            fixedrange=True,
            autorange='reversed'
        ),
        barmode='overlay',
        showlegend=False,
        plot_bgcolor='white',
        margin=dict(l=250, r=100, t=40, b=50),
        shapes=shapes,
        height=final_height
    )

    config = {
        'displayModeBar': False,
        'scrollZoom': False,
        'staticPlot': False,
        'doubleClick': 'reset'
    }
    fig = go.Figure(data=traces, layout=layout)
    html_output = pio.to_html(fig, full_html=False, config=config, validate=False)
    return html_output


def plot_features(fasta_df, selected_protein_id, custom_features_df=None, custom_label=None):
    protein_features = fetch_protein_features(selected_protein_id)
    protein_sequence = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'sequence'].iloc[0]
    protein_length = len(protein_sequence)

    if not protein_features:
        raise ValueError('No features found for the selected protein.')

    feature_groups = OrderedDict()
    feature_traces = []
    render_features = []

    feature_bar_height = 0.5
    feature_bar_line_width = 1
    feature_label_size = 14
    global_features_height = 50

    for feature in protein_features:
        feature_length = feature['end'] - feature['start']
        feature_length_offset = 0
        feature_position = f"{feature['start']} - {feature['end']}"

        if feature_length == 0:
            feature_length_offset = 1
            feature_position = f"{feature['start']}"

        if feature['type'] == 'VARIANT' and feature.get('alternativeSequence', '') == '':
            continue

        render_features.append({
            'group': feature['group'],
            'type': feature['type'],
            'description': feature.get('description', ''),
            'molecule': feature.get('molecule'),
            'ligand': feature.get('ligand'),
            'ftID': feature.get('ftID'),
            'alternativeSequence': feature.get('alternativeSequence'),
            'start': feature['start'],
            'end': feature['end'],
            'length': feature_length + feature_length_offset,
            'position': feature_position
        })

    if custom_features_df is not None and custom_label:
        filtered_df = custom_features_df[
            custom_features_df['uniprot_id'].astype(str).str.upper() == selected_protein_id.upper()
        ]
        if not filtered_df.empty:
            custom_group = custom_label
            for _, feature in filtered_df.iterrows():
                feature_length = feature['end'] - feature['start']
                feature_length_offset = 0
                feature_position = f"{feature['start']} - {feature['end']}"

                if feature_length == 0:
                    feature_length_offset = 1
                    feature_position = f"{feature['start']}"

                render_features.append({
                    'group': custom_group,
                    'type': 'CUSTOM',
                    'description': feature['description'] or 'Custom Feature',
                    'literature': feature['literature'],
                    'start': feature['start'],
                    'end': feature['end'],
                    'length': feature_length + feature_length_offset,
                    'position': feature_position
                })

    if not render_features:
        raise ValueError('No features found for the selected protein.')

    group_priority = {
        'Topology': 0,
        'Transmembrane': 1,
        'Domains': 2,
        'Binding Sites': 3,
        'Sites': 4,
        'Modified Residues': 5,
        'Variants': 6
    }
    unique_groups = []
    seen_groups = set()
    for feature in render_features:
        group = feature['group']
        if group in seen_groups:
            continue
        seen_groups.add(group)
        unique_groups.append(group)

    ordered_groups = sorted(unique_groups, key=lambda g: group_priority.get(g, 50))
    for group in ordered_groups:
        feature_groups[group] = len(feature_groups)

    for feature in render_features:
        group = feature['group']
        feature_position = feature['position']

        if feature['type'] == 'DOMAIN':
            feature_trace = go.Bar(
                x=[feature['length']],
                y=[feature_groups[group]],
                base=feature['start'],
                orientation='h',
                width=feature_bar_height,
                marker=dict(
                    color='lightblue',
                    line=dict(color='blue', width=feature_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=(
                    f"{feature['description']}"
                    f"<br>Position: {feature_position}"
                ),
                hoverlabel=dict(align='left')
            )
        elif feature['type'] == 'BINDING':
            feature_trace = go.Bar(
                x=[feature['length']],
                y=[feature_groups[group]],
                base=feature['start'],
                orientation='h',
                width=feature_bar_height,
                marker=dict(
                    color='lightgreen',
                    line=dict(color='green', width=feature_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=(
                    f"{feature['description']}"
                    f"<br>Position: {feature_position}"
                    f"<br>Molecule: {feature['molecule']}"
                    f"<br>Ligand: {feature['ligand']}"
                ),
                hoverlabel=dict(align='left')
            )
        elif feature['type'] == 'MOD_RES':
            feature_trace = go.Bar(
                x=[feature['length']],
                y=[feature_groups[group]],
                base=feature['start'],
                orientation='h',
                width=feature_bar_height,
                marker=dict(
                    color='darkred',
                    line=dict(color='red', width=feature_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=(
                    f"{feature['description']}"
                    f"<br>Position: {feature_position}"
                ),
                hoverlabel=dict(align='left')
            )
        elif feature['type'] == 'SITE':
            feature_trace = go.Bar(
                x=[feature['length']],
                y=[feature_groups[group]],
                base=feature['start'],
                orientation='h',
                width=feature_bar_height,
                marker=dict(
                    color='lightcoral',
                    line=dict(color='red', width=feature_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=(
                    f"{feature['description']}"
                    f"<br>Position: {feature_position}"
                ),
                hoverlabel=dict(align='left')
            )
        elif feature['type'] == 'VARIANT':
            consensus_aa = protein_sequence[feature['start'] - 1:feature['end']]
            feature_trace = go.Bar(
                x=[feature['length']],
                y=[feature_groups[group]],
                base=feature['start'],
                orientation='h',
                width=feature_bar_height,
                marker=dict(
                    color='lightgray',
                    line=dict(color='gray', width=feature_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=(
                    f"{feature['description']}"
                    f"<br>Position: {feature_position}"
                    f"<br>SAAV: {consensus_aa}->{feature['alternativeSequence']}"
                    f"<br>Feature ID: {feature['ftID']}"
                ),
                hoverlabel=dict(align='left')
            )
        elif feature['type'] == 'TOPO_DOM':
            feature_trace = go.Bar(
                x=[feature['length']],
                y=[feature_groups[group]],
                base=feature['start'],
                orientation='h',
                width=feature_bar_height,
                marker=dict(
                    color='lightyellow',
                    line=dict(color='goldenrod', width=feature_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=(
                    f"{feature['description']}"
                    f"<br>Position: {feature_position}"
                ),
                hoverlabel=dict(align='left')
            )
        elif feature['type'] == 'TRANSMEM':
            feature_trace = go.Bar(
                x=[feature['length']],
                y=[feature_groups[group]],
                base=feature['start'],
                orientation='h',
                width=feature_bar_height,
                marker=dict(
                    color='paleturquoise',
                    line=dict(color='teal', width=feature_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=(
                    f"{feature['description']}"
                    f"<br>Position: {feature_position}"
                ),
                hoverlabel=dict(align='left')
            )
        elif feature['type'] == 'CUSTOM':
            literature = feature.get('literature')
            hover_lines = [feature['description'], f"Position: {feature_position}"]
            if literature and str(literature).lower() != 'nan':
                hover_lines.append(f"Literature: {literature}")
            hovertext = '<br>'.join(hover_lines)
            feature_trace = go.Bar(
                x=[feature['length']],
                y=[feature_groups[group]],
                base=feature['start'],
                orientation='h',
                width=feature_bar_height,
                marker=dict(
                    color='wheat',
                    line=dict(color='saddlebrown', width=feature_bar_line_width)
                ),
                hoverinfo='text',
                hovertext=hovertext,
                hoverlabel=dict(align='left')
            )
        else:
            continue
        feature_traces.append(feature_trace)

    layout = go.Layout(
        xaxis=dict(
            range=[1, protein_length],
            tickvals=[''],
            ticktext=[''],
            fixedrange=True
        ),
        yaxis=dict(
            tickmode='array',
            tickvals=list(feature_groups.values()),
            ticktext=list(feature_groups.keys()),
            tickfont=dict(size=feature_label_size),
            fixedrange=True
        ),
        barmode='stack',
        showlegend=False,
        plot_bgcolor='white',
        margin=dict(l=250, r=100, t=0, b=0),
        height=global_features_height * len(feature_groups)
    )

    config = {
        'displayModeBar': False,
        'scrollZoom': False,
        'staticPlot': False,
        'doubleClick': 'reset'
    }
    fig = go.Figure(data=feature_traces, layout=layout)
    html_output = pio.to_html(fig, full_html=False, config=config, validate=False)
    return html_output


def render_tabs(tab_id, items):
    parts = [f'<div class="pepmap-tabs" data-tab-group="{tab_id}">']
    parts.append('<div class="flex flex-wrap gap-2 pr-24">')
    for idx, item in enumerate(items):
        label = html.escape(item['label'])
        active = idx == 0
        base_classes = (
            'tab-btn rounded-full border px-3 py-1 text-xs font-semibold transition '
            'border-slate-300'
        )
        if active:
            classes = f'{base_classes} bg-slate-300 text-slate-900'
            active_attr = ' data-active="true"'
        else:
            classes = f'{base_classes} bg-white text-slate-700 hover:border-slate-400'
            active_attr = ''
        parts.append(
            f'<button type="button" class="{classes}" data-tab-target="{tab_id}-panel-{idx}" '
            f'data-index="{idx}"{active_attr}>{label}</button>'
        )
    parts.append('</div>')
    parts.append('<div class="tab-panels mt-3">')
    for idx, item in enumerate(items):
        panel_class = '' if idx == 0 else 'hidden'
        parts.append(
            f'<div id="{tab_id}-panel-{idx}" class="tab-panel {panel_class}">{item["content"]}</div>'
        )
    parts.append('</div></div>')
    return ''.join(parts)


def _render_summary_table(header_cells, rows):
    header_html = ''.join(
        f'<th class="border border-slate-200 bg-slate-50 px-3 py-2 text-left font-semibold">'
        f'{html.escape(str(cell))}</th>'
        for cell in header_cells
    )

    body_rows = []
    for row in rows:
        cells = []
        for cell_key in header_cells:
            value = row.get(cell_key, '')
            cells.append(
                f'<td class="border border-slate-200 px-3 py-2 text-left">{html.escape(str(value))}</td>'
            )
        body_rows.append(f'<tr>{"".join(cells)}</tr>')

    return (
        '<div class="overflow-x-auto">'
        '<table class="min-w-full border border-slate-200">'
        f'<thead><tr>{header_html}</tr></thead>'
        f'<tbody>{"".join(body_rows)}</tbody>'
        '</table>'
        '</div>'
    )


def build_peptide_summary_table(
    summary_mode,
    counts_by_protein,
    all_runs,
    column_order=None,
    feature_counts_by_protein=None,
    feature_labels_by_protein=None
):
    if not counts_by_protein:
        return ''
    all_runs = list(all_runs or [])
    if not all_runs:
        return ''

    if column_order:
        ordered_labels = [label for label in column_order if label in counts_by_protein]
        remaining = [label for label in counts_by_protein.keys() if label not in ordered_labels]
        label_iter = ordered_labels + remaining
    else:
        label_iter = list(counts_by_protein.keys())

    if summary_mode == 'per_sample':
        rows = []
        for run in all_runs:
            row = {'Run': run}
            total = 0
            for label in label_iter:
                count = counts_by_protein.get(label, {}).get(run, 0)
                row[label] = int(count)
                total += int(count)
            row['Total'] = total
            rows.append(row)

        header_cells = ['Run'] + label_iter + ['Total']
        table_html = _render_summary_table(header_cells, rows)
        return (
            '<div class="peptide-summary mt-6 text-xs text-slate-700">'
            '<div class="text-sm font-semibold text-center text-slate-600 mb-2">Peptides per sample</div>'
            f'{table_html}'
            '</div>'
        )

    feature_counts_by_protein = feature_counts_by_protein or {}
    feature_labels_by_protein = feature_labels_by_protein or {}
    title = (
        'Peptides per topology per sample'
        if summary_mode == 'per_topology'
        else 'Peptides per domain per sample'
    )

    parts = [f'<div class="text-sm font-semibold text-center text-slate-600 mb-2">{title}</div>']
    multiple_proteins = len(label_iter) > 1
    for label in label_iter:
        feature_labels = feature_labels_by_protein.get(label, [])
        if multiple_proteins:
            parts.append(
                f'<div class="text-sm font-semibold text-slate-600 mt-4 mb-2">{html.escape(label)}</div>'
            )
        if not feature_labels:
            parts.append('<div class="text-slate-500">No features found for this protein</div>')
            continue

        rows = []
        for run in all_runs:
            row = {'Run': run}
            for feature_label in feature_labels:
                count = feature_counts_by_protein.get(label, {}).get(feature_label, {}).get(run, 0)
                row[feature_label] = int(count)
            rows.append(row)

        header_cells = ['Run'] + feature_labels
        parts.append(_render_summary_table(header_cells, rows))

    return (
        '<div class="peptide-summary mt-6 text-xs text-slate-700">'
        f'{"".join(parts)}'
        '</div>'
    )
