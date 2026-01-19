import os
import re
import json
import time
import uuid
import tempfile
import numpy as np
import shutil
import pandas as pd
import plotly.io as pio
import plotly.graph_objs as go
import requests
from io import StringIO
from pathlib import Path
from flask_caching import Cache
from flask import Flask, render_template, request, jsonify
from apscheduler.schedulers.background import BackgroundScheduler
from collections import OrderedDict

app = Flask(__name__)
app.config['CACHE_TYPE'] = 'SimpleCache'
app.config['SECRET_KEY'] = os.urandom(24)
cache = Cache(app)

__version__ = "1.0.0"

PEPMAP_STORE_DIR = Path(os.environ.get('PEPMAP_STORE_DIR', 'storage'))
PEPMAP_STORE_DIR.mkdir(parents=True, exist_ok=True)
PEPMAP_STORE_TTL_SECONDS = int(os.environ.get('PEPMAP_STORE_TTL_SECONDS', '1800'))


class FileSessionStore:
    def __init__(self, root, ttl_seconds=1800):
        self.root = root
        self.ttl_seconds = ttl_seconds

    def _path(self, sid):
        return self.root / f"{sid}.json"

    def _write(self, sid, payload):
        tmp = tempfile.NamedTemporaryFile('w', delete=False, dir=str(self.root))
        try:
            json.dump(payload, tmp)
            tmp.flush()
            os.fsync(tmp.fileno())
        finally:
            tmp.close()
        os.replace(tmp.name, self._path(sid))

    def create(self, meta=None):
        sid = uuid.uuid4().hex
        payload = {
            'created_at': time.time(),
            'updated_at': time.time(),
            'ttl': self.ttl_seconds,
            'meta': meta or {},
            'data': {}
        }
        self._write(sid, payload)
        return sid

    def create_with_id(self, sid, meta=None):
        payload = {
            'created_at': time.time(),
            'updated_at': time.time(),
            'ttl': self.ttl_seconds,
            'meta': meta or {},
            'data': {}
        }
        self._write(sid, payload)
        return sid

    def get_payload(self, sid):
        path = self._path(sid)
        if not path.exists():
            return None
        try:
            payload = json.loads(path.read_text())
        except Exception:
            return None
        ttl = int(payload.get('ttl', self.ttl_seconds))
        updated_at = float(payload.get('updated_at', 0))
        if time.time() - updated_at > ttl:
            self.delete(sid)
            return None
        return payload

    def update(self, sid, **kv):
        payload = self.get_payload(sid)
        if payload is None:
            return False
        payload['updated_at'] = time.time()
        payload.setdefault('data', {}).update(kv)
        self._write(sid, payload)
        return True

    def read(self, sid, key):
        payload = self.get_payload(sid)
        if payload is None:
            return None
        return payload.get('data', {}).get(key)

    def delete(self, sid):
        try:
            self._path(sid).unlink()
        except FileNotFoundError:
            pass

    def purge_expired(self):
        now = time.time()
        for f in self.root.glob('*.json'):
            try:
                payload = json.loads(f.read_text())
                ttl = int(payload.get('ttl', self.ttl_seconds))
                updated_at = float(payload.get('updated_at', 0))
                if now - updated_at > ttl:
                    f.unlink()
            except Exception:
                pass


store = FileSessionStore(PEPMAP_STORE_DIR, ttl_seconds=PEPMAP_STORE_TTL_SECONDS)


def normalize_session_id(raw_session_id):
    if raw_session_id is None:
        return None
    session_id = raw_session_id.strip()
    if not session_id:
        return None
    if not re.match(r'^[A-Za-z0-9_-]{8,128}$', session_id):
        return None
    return session_id


def get_request_session_id():
    if request.method == 'GET':
        raw_session_id = request.args.get('session_id')
    else:
        raw_session_id = request.form.get('session_id')
    if raw_session_id is None:
        return None, None
    session_id = normalize_session_id(raw_session_id)
    if session_id is None:
        return None, 'Invalid session_id format.'
    return session_id, None


def generate_dynamic_ticks(protein_length):
    # define ranges and steps
    ranges = [
        (0, 500, 10),
        (501, 1000, 20),
        (1001, 2000, 50),
        (2001, float('inf'), 100)
    ]

    # determine the appropriate step based on protein length
    for start, end, step in ranges:
        if start <= protein_length < end:
            chosen_step = step
            break
    else:
        chosen_step = 100

    # generate tick values and labels using the chosen step
    tickvals = list(range(0, protein_length + 1, chosen_step))
    ticktext = [str(tick) for tick in tickvals]

    # change first tick value to 1
    tickvals[0] = 1
    ticktext[0] = '1'

    # add the last tick value
    if protein_length % chosen_step != 0:
        tickvals.append(protein_length)
        ticktext.append('')

    return tickvals, ticktext


def fetch_protein_features(uniprot_id):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession={uniprot_id}"
    response = requests.get(url, headers={"Accept": "application/json"})
    if response.status_code == 200:
        data = response.json()
        features = []
        for entry in data:
            if 'features' in entry:
                for feature in entry['features']:
                    if feature['type'] == 'DOMAIN':
                        features.append({
                            'group': 'Domains',
                            'type': feature['type'],
                            'description': feature['description'],
                            'start': int(feature['begin']),
                            'end': int(feature['end']),
                        })
                    elif feature['type'] == 'BINDING':
                        features.append({
                            'group': 'Binding Sites',
                            'type': feature['type'],
                            'description': feature['description'],
                            'molecule': feature['molecule'],
                            'ligand': feature['ligand']['name'] if 'ligand' in feature else None,
                            'start': int(feature['begin']),
                            'end': int(feature['end']),
                        })
                    elif feature['type'] == 'SITE':
                        features.append({
                            'group': 'Sites',
                            'type': feature['type'],
                            'description': feature['description'],
                            'start': int(feature['begin']),
                            'end': int(feature['end']),
                        })
                    elif feature['type'] == 'MOD_RES':
                        features.append({
                            'group': 'Modified Residues',
                            'type': feature['type'],
                            'description': feature['description'],
                            'start': int(feature['begin']),
                            'end': int(feature['end']),
                        })
                    elif feature['type'] == 'VARIANT':
                        features.append({
                            'group': 'Variants',
                            'type': feature['type'],
                            'ftID': feature['ftId'],
                            'description': feature['description'],
                            'alternativeSequence': feature['alternativeSequence'],
                            'start': int(feature['begin']),
                            'end': int(feature['end']),
                        })
        return features
    else:
        return []


def parse_fasta(fasta_stream, organism):
    entries = []
    sequence = ''
    uniprot_id = ''
    gene_symbol = ''
    for line in fasta_stream:
        line = line.decode('utf-8').strip()
        # remove entries for reverse and contaminant sequences of FragPipe FASTAs
        if line.startswith('>rev_sp') or line.startswith('>contam_sp'):
            continue
        elif line.startswith('>'):
            if sequence:
                entries.append({'uniprot_id': uniprot_id, 'gene_symbol': gene_symbol, 'sequence': sequence})
            sequence = ''
            parts = line[1:].split('|')
            try:
                uniprot_id = parts[1]
            except IndexError:
                # Using full header instead.
                uniprot_id = parts
            try:
                gene_symbol = parts[2].split(' ')[0].split(f"_{organism}")[0]
            except IndexError:
                # Using full gene name instead.
                gene_symbol = parts
        else:
            sequence += line
    if sequence:
        entries.append({'uniprot_id': uniprot_id, 'gene_symbol': gene_symbol, 'sequence': sequence})

    fasta_df = pd.DataFrame(entries)
    return fasta_df


def normalize_search_input(search_input):
    if search_input is None:
        return ''
    return search_input.strip()


def parse_search_inputs(search_input):
    if search_input is None:
        return []
    tokens = [token.strip() for token in re.split(r'[\s,;]+', search_input) if token.strip()]
    seen = set()
    results = []
    for token in tokens:
        key = token.upper()
        if key in seen:
            continue
        seen.add(key)
        results.append(token)
    return results


def parse_search_labels(raw_labels, expected_count):
    if not raw_labels:
        return [''] * expected_count
    try:
        labels = json.loads(raw_labels)
    except Exception:
        return [''] * expected_count
    if not isinstance(labels, list):
        return [''] * expected_count
    normalized = []
    for value in labels:
        label = ''
        if isinstance(value, str):
            label = value.strip()
        normalized.append(label)
    if len(normalized) < expected_count:
        normalized.extend([''] * (expected_count - len(normalized)))
    return normalized[:expected_count]


def normalize_organism(organism):
    if organism is None:
        return None
    organism_value = organism.strip().upper()
    if organism_value in {'HUMAN', 'MOUSE'}:
        return organism_value
    return None


@app.route('/resolve_labels', methods=['POST'])
def resolve_labels():
    payload = request.get_json(silent=True) or {}
    session_id = normalize_session_id(payload.get('session_id'))
    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400
    identifiers = payload.get('identifiers') or []
    if not isinstance(identifiers, list):
        return jsonify({'error': 'identifiers must be a list.'}), 400

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


def find_uniprot_id_by_gene_symbol(fasta_df, gene_symbol):
    if not gene_symbol:
        return None
    matches = fasta_df[fasta_df['gene_symbol'].astype(str).str.upper() == gene_symbol.upper()]
    if matches.empty:
        return None
    return matches['uniprot_id'].iloc[0]


def find_uniprot_id_by_accession(fasta_df, accession):
    if not accession:
        return None
    matches = fasta_df[fasta_df['uniprot_id'].astype(str).str.upper() == accession.upper()]
    if matches.empty:
        return None
    return matches['uniprot_id'].iloc[0]


def apply_sample_name_cleanup(report_df, cleanup_mode, custom_split_pattern):
    if cleanup_mode not in {'split_underscore', 'custom'}:
        return report_df

    cleaned_df = report_df.copy()
    if cleanup_mode == 'split_underscore':
        cleaned_df['Run'] = cleaned_df['Run'].astype(str).apply(
            lambda value: value.split('_', 1)[0]
        )
    else:
        pattern = (custom_split_pattern or '').strip()
        if pattern:
            cleaned_df['Run'] = cleaned_df['Run'].astype(str).apply(
                lambda value: value.split(pattern, 1)[0]
            )
    return cleaned_df


def detect_file_type(tsv_stream):
    header = tsv_stream.readline().decode('utf-8').strip().split('\t')
    tsv_stream.seek(0)

    if 'Run' in header and 'Protein.Ids' in header and 'Precursor.Normalised' in header:
        return 'diann'
    elif 'Spectrum' in header and 'Peptide' in header and 'Protein ID' in header:
        return 'fragpipe'
    else:
        raise ValueError("Unable to determine file type. Please ensure it's a valid DIA-NN or FragPipe report.")


def parse_report_tsv(tsv_stream):
    file_type = detect_file_type(tsv_stream)

    if file_type == 'diann':
        report_df = pd.read_csv(tsv_stream, delimiter='\t')
        report_df = report_df.rename(columns={'Q.Value': 'P.Value (Q.Value)'})
        report_df = report_df[['Run', 'Protein.Ids', 'Precursor.Normalised', 'Stripped.Sequence',
                               'Precursor.Charge', 'P.Value (Q.Value)', 'Proteotypic']]
    elif file_type == 'fragpipe':
        report_df = pd.read_csv(tsv_stream, delimiter='\t')
        report_df = report_df.rename(columns={
            'Spectrum': 'Run',
            'Peptide': 'Stripped.Sequence',
            'Charge': 'Precursor.Charge',
            'Expectation': 'P.Value (Expectation)',
            'Intensity': 'Precursor.Normalised',
            'Is Unique': 'Proteotypic',
            'Protein ID': 'Protein.Ids'
        })

        # strip . from behind (3 times) and keep first (to get proper Sample names)
        report_df['Run'] = report_df['Run'].str.rsplit('.', n=3).str[0]

        # drop Purity < 0.5 (in line with TMT-Integrator)
        report_df = report_df[report_df['Purity'] >= 0.5]

        # sum up all intensities for the same peptide in the same Run with the same Charge
        def aggregate_q_value(group):
            min_val = group.min()
            max_val = group.max()
            count = len(group)

            if count == 1:
                return f"{min_val:.2E}"
            else:
                return f"{min_val:.2E} - {max_val:.2E} (count: {count})"

        report_df = report_df.groupby(['Run', 'Stripped.Sequence', 'Precursor.Charge'], as_index=False).agg({
            'Precursor.Normalised': 'sum',
            'P.Value (Expectation)': aggregate_q_value,
            'Proteotypic': 'first',
            'Protein.Ids': 'first',
            'Mapped Proteins': 'first'
        })

        # extract primary protein ID from 'Protein ID' column
        report_df['Protein.Ids'] = report_df['Protein.Ids'].apply(lambda x: x.split('|')[1] if '|' in x else x)

        # format 'Protein.Ids' column to match DIA-NN format
        def extract_protein_ids(mapped_proteins):
            return ';'.join([p.split('|')[1] for p in mapped_proteins.split(',') if '|' in p])

        report_df['Protein.Ids'] = report_df.apply(
            lambda row: f"{row['Protein.Ids']};{extract_protein_ids(row['Mapped Proteins'])}"
            if pd.notna(row['Mapped Proteins']) else row['Protein.Ids'],
            axis=1
        )

        report_df = report_df[['Run', 'Protein.Ids', 'Precursor.Normalised', 'Stripped.Sequence',
                               'Precursor.Charge', 'P.Value (Expectation)', 'Proteotypic']]

    else:
        raise ValueError("Invalid file type")

    # drop rows with 0 in Precursor.Normalised
    report_df = report_df[report_df['Precursor.Normalised'] != 0]

    return report_df


def parse_custom_features_tsv(tsv_stream):
    custom_df = pd.read_csv(tsv_stream, delimiter='\t')
    custom_df.columns = [col.strip() for col in custom_df.columns]
    column_map = {col.strip().lower(): col for col in custom_df.columns}
    required = {'uniprot', 'position', 'description', 'literature'}
    missing = sorted(req for req in required if req not in column_map)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    custom_df = custom_df.rename(columns={
        column_map['uniprot']: 'uniprot_id',
        column_map['position']: 'position',
        column_map['description']: 'description',
        column_map['literature']: 'literature'
    })

    def parse_position_value(value):
        if pd.isna(value):
            return None, None
        text = str(value).strip()
        if not text:
            return None, None
        range_match = re.match(r'^(\d+)\s*-\s*(\d+)$', text)
        if range_match:
            return int(range_match.group(1)), int(range_match.group(2))
        single_match = re.match(r'^(\d+)$', text)
        if single_match:
            position = int(single_match.group(1))
            return position, position
        return None, None

    def normalize_text(value):
        if pd.isna(value):
            return ''
        return str(value).strip()

    rows = []
    for index, row in custom_df.iterrows():
        start, end = parse_position_value(row['position'])
        if start is None or end is None:
            raise ValueError(f"Invalid position at row {index + 2}: {row['position']}")
        if start > end:
            raise ValueError(f"Start position exceeds end at row {index + 2}: {row['position']}")
        uniprot_id = normalize_text(row['uniprot_id'])
        if not uniprot_id:
            raise ValueError(f"Missing UniProt ID at row {index + 2}.")
        rows.append({
            'uniprot_id': uniprot_id,
            'start': start,
            'end': end,
            'description': normalize_text(row['description']),
            'literature': normalize_text(row['literature'])
        })

    parsed_df = pd.DataFrame(rows)
    if parsed_df.empty:
        raise ValueError("No valid custom features found in the uploaded file.")
    return parsed_df


def plot_peptides(peptide_positions_df, fasta_df, selected_protein_id, global_log2_min, global_log2_max, p_value_column, p_value_name, custom_title):
    # get the protein sequence and length
    protein_sequence = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'sequence'].iloc[0]
    protein_length = len(protein_sequence)

    # calculate normalized intensity
    peptide_positions_df['log2_intensity'] = np.log2(peptide_positions_df['Precursor.Normalised'])
    peptide_positions_df['normalized_intensity'] = (peptide_positions_df['log2_intensity'] - global_log2_min) / (global_log2_max - global_log2_min)

    # get unique runs
    unique_runs = sorted(peptide_positions_df['Run'].unique())

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

        # assign peptides to y-offsets to avoid overlaps
        for idx, row in group.iterrows():
            placed = False
            for y_offset_idx, end_time in enumerate(end_positions_per_y_offset):
                if row['Start'] >= end_time:
                    # Place peptide at this y-offset
                    end_positions_per_y_offset[y_offset_idx] = row['End']
                    y_offsets.append(y_offset_idx)
                    placed = True
                    break
            if not placed:
                # create a new y-offset
                end_positions_per_y_offset.append(row['End'])
                y_offsets.append(len(end_positions_per_y_offset) - 1)

        # calculate overlap offsets
        group['overlap_offset'] = [offset * (peptide_bar_height + peptide_bar_margin) for offset in y_offsets]
        max_y_offset = max(y_offsets) if y_offsets else 0

        # define y_base for the current run
        y_base = current_y
        group['y_base'] = y_base

        # add y-tick values and labels
        # position the label at the center of the run's vertical space
        run_height = max((max_y_offset + 1) * (peptide_bar_height + peptide_bar_margin), min_run_height)
        y_tickvals.append(y_base + run_height / 2)
        y_ticktext.append(run)

        for idx, row in group.iterrows():
            color = f'rgba(255,{255 - row["normalized_intensity"] * 255},0,0.8)'
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
                hovertext=f'<b>{row["Peptide"]}</b>'
                          f'<br>Position: {row["Start"]}-{row["End"]}'
                          f'<br>Log2 Intensity: {row["log2_intensity"]:.2f}'
                          f'<br>Charge: {row["Precursor.Charge"]}'
                          f'<br>{p_value_name}: {row[p_value_column]}'
                          f'<br>Proteotypic: {"Yes" if row["Proteotypic"] else "No"}',
                showlegend=False
            )
            traces.append(trace)

        # add separator line after each run
        shapes.append({
            'type': 'line',
            'x0': 0,
            'y0': y_base + run_height,
            'x1': protein_length,
            'y1': y_base + run_height,
            'line': {
                'color': 'black',
                'width': 1,
            },
        })

        # update current_y for the next run
        current_y += run_height + peptide_bar_margin * 4 # Add extra space between runs

    # generate dynamic ticks for the x-axis
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

    fig = go.Figure(data=traces, layout=layout)

    config = {
        'displayModeBar': False,
        'scrollZoom': False,
        'staticPlot': False,
        'doubleClick': 'reset',
    }

    plot_peptides_html = pio.to_html(fig, full_html=False, config=config)
    return plot_peptides_html


def plot_features(fasta_df, selected_protein_id, custom_features_df=None, custom_label=None):
    protein_features = fetch_protein_features(selected_protein_id)
    protein_sequence = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'sequence'].iloc[0]
    protein_length = len(protein_sequence)

    feature_groups = OrderedDict()

    feature_bar_height = 0.5
    feature_bar_line_width = 1
    feature_label_size = 14
    global_features_height = 50

    feature_traces = []

    if not protein_features:
        return jsonify({'error': 'No features found for the selected protein.'}), 400
    else:
        for feature in protein_features:
            feature_length = feature['end'] - feature['start']
            feature_length_offset = 0
            feature_position = f"{feature['start']} - {feature['end']}"

            if feature_length == 0:
                feature_length_offset = 1
                feature_position = f"{feature['start']}"

            group = feature['group']
            if group not in feature_groups:
                feature_groups[group] = len(feature_groups)

            if feature['type'] == 'DOMAIN':
                # create a bar for each domain feature
                feature_trace = go.Bar(
                    x=[feature_length + feature_length_offset],
                    y=[feature_groups[group]],
                    base=feature['start'],
                    orientation='h',
                    width=feature_bar_height,
                    marker=dict(
                        color='lightblue',
                        line=dict(color='blue', width=feature_bar_line_width)
                    ),
                    hoverinfo='text',
                    hovertext=f"{feature['description']}"
                              f"<br>Position: {feature_position}",
                    hoverlabel=dict(align='left')
                )
            elif feature['type'] == 'BINDING':
                # create a bar for each binding site
                feature_trace = go.Bar(
                    x=[feature_length + feature_length_offset],
                    y=[feature_groups[group]],
                    base=feature['start'],
                    orientation='h',
                    width=feature_bar_height,
                    marker=dict(
                        color='lightgreen',
                        line=dict(color='green', width=feature_bar_line_width)
                    ),
                    hoverinfo='text',
                    hovertext=f"{feature['description']}"
                              f"<br>Position: {feature_position}"
                              f"<br>Molecule: {feature['molecule']}"
                              f"<br>Ligand: {feature['ligand']}",
                    hoverlabel=dict(align='left')
                )
            elif feature['type'] == 'MOD_RES':
                # create a bar for each modified residue
                feature_trace = go.Bar(
                    x=[feature_length + feature_length_offset],
                    y=[feature_groups[group]],
                    base=feature['start'],
                    orientation='h',
                    width=feature_bar_height,
                    marker=dict(
                        color='darkred',
                        line=dict(color='red', width=feature_bar_line_width)
                    ),
                    hoverinfo='text',
                    hovertext=f"{feature['description']}"
                              f"<br>Position: {feature_position}",
                    hoverlabel=dict(align='left')
                )
            elif feature['type'] == 'SITE':
                # create a bar for each site feature
                feature_trace = go.Bar(
                    x=[feature_length + feature_length_offset],
                    y=[feature_groups[group]],
                    base=feature['start'],
                    orientation='h',
                    width=feature_bar_height,
                    marker=dict(
                        color='lightcoral',
                        line=dict(color='red', width=feature_bar_line_width)
                    ),
                    hoverinfo='text',
                    hovertext=f"{feature['description']}"
                              f"<br>Position: {feature_position}",
                    hoverlabel=dict(align='left')
                )
            elif feature['type'] == 'VARIANT':
                if feature['alternativeSequence'] == "":
                    continue
                # create a bar for each variant feature
                consensus_aa = protein_sequence[feature['start']-1:feature['end']]

                feature_trace = go.Bar(
                    x=[feature_length + feature_length_offset],
                    y=[feature_groups[group]],
                    base=feature['start'],
                    orientation='h',
                    width=feature_bar_height,
                    marker=dict(
                        color='lightgray',
                        line=dict(color='gray', width=feature_bar_line_width)
                    ),
                    hoverinfo='text',
                    hovertext=f"{feature['description']}"
                              f"<br>Position: {feature_position}"
                              f"<br>SAAV: {consensus_aa}->{feature['alternativeSequence']}"
                              f"<br>Feature ID: {feature['ftID']}",
                    hoverlabel=dict(align='left')
                )
            feature_traces.append(feature_trace)

    if custom_features_df is not None and custom_label:
        filtered_df = custom_features_df[
            custom_features_df['uniprot_id'].astype(str).str.upper() == selected_protein_id.upper()
        ]
        if not filtered_df.empty:
            custom_group = custom_label
            if custom_group not in feature_groups:
                feature_groups[custom_group] = len(feature_groups)

            for _, feature in filtered_df.iterrows():
                feature_length = feature['end'] - feature['start']
                feature_length_offset = 0
                feature_position = f"{feature['start']} - {feature['end']}"

                if feature_length == 0:
                    feature_length_offset = 1
                    feature_position = f"{feature['start']}"

                description = feature['description'] or 'Custom Feature'
                literature = feature['literature']
                hover_lines = [description, f"Position: {feature_position}"]
                if literature and literature.lower() != 'nan':
                    hover_lines.append(f"Literature: {literature}")
                hovertext = "<br>".join(hover_lines)

                feature_traces.append(
                    go.Bar(
                        x=[feature_length + feature_length_offset],
                        y=[feature_groups[custom_group]],
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
                )

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
    return pio.to_html(fig, full_html=False, config=config)


def render_tabs(tab_id, items):
    tab_rules = []
    for idx in range(len(items)):
        input_id = f"{tab_id}-tab-{idx}"
        panel_id = f"{tab_id}-panel-{idx}"
        tab_rules.append(f"#{tab_id} #{input_id}:checked ~ .tab-panels #{panel_id}{{display:block;}}")
    tab_style = (
        "<style>"
        f"#{tab_id}{{margin-top:8px;}}"
        f"#{tab_id} input[type=radio]{{display:none;}}"
        f"#{tab_id} .tab-label{{display:inline-block;padding:6px 10px;border:1px solid #ccc;"
        "border-radius:4px;background:#f3f3f3;cursor:pointer;font-size:12px;margin:0 6px 6px 0;}}"
        f"#{tab_id} input[type=radio]:checked + .tab-label{{background:#e7e7e7;font-weight:600;}}"
        f"#{tab_id} .tab-panels{{border:none;border-radius:0;padding:6px;background:#fff;}}"
        f"#{tab_id} .tab-panel{{display:none;}}"
        + "".join(tab_rules) +
        "</style>"
    )
    parts = [tab_style, f'<div class="pepmap-tabs" id="{tab_id}">']
    for idx, item in enumerate(items):
        safe_label = item['label']
        input_id = f"{tab_id}-tab-{idx}"
        checked = ' checked' if idx == 0 else ''
        parts.append(
            f'<input type="radio" name="{tab_id}-tabs" id="{input_id}" data-index="{idx}"{checked}>'
        )
        parts.append(f'<label class="tab-label" for="{input_id}">{safe_label}</label>')
    parts.append('<div class="tab-panels">')
    for idx, item in enumerate(items):
        panel_id = f"{tab_id}-panel-{idx}"
        parts.append(f'<div class="tab-panel" id="{panel_id}">{item["content"]}</div>')
    parts.append('</div></div>')
    return ''.join(parts)


def build_peptide_summary_table(counts_by_protein, all_runs, column_order=None):
    if not counts_by_protein:
        return ''
    all_runs = list(all_runs or [])
    if not all_runs:
        return ''
    summary_df = pd.DataFrame(index=all_runs)
    if column_order:
        ordered_labels = [label for label in column_order if label in counts_by_protein]
        remaining = [label for label in counts_by_protein.keys() if label not in ordered_labels]
        label_iter = ordered_labels + remaining
    else:
        label_iter = counts_by_protein.keys()
    for label in label_iter:
        counts = counts_by_protein.get(label, {})
        summary_df[label] = [counts.get(run, 0) for run in all_runs]
    summary_df['Total'] = summary_df.sum(axis=1)
    summary_df.index.name = ''
    summary_df = summary_df.fillna(0).astype(int)
    table_html = summary_df.to_html(classes='peptide-summary-table', border=1, index_names=False)
    return (
        "<style>"
        ".peptide-summary{margin-top:12px;font-size:12px;text-align:left;}"
        ".peptide-summary-title{font-weight:600;margin-bottom:6px;text-align:center;}"
        ".peptide-summary-table{border-collapse:collapse;width:90%;margin:0 auto;}"
        ".peptide-summary-table th,.peptide-summary-table td{border:1px solid #ddd;padding:4px 6px;text-align:left;}"
        ".peptide-summary-table th{background:#f5f5f5;}"
        "</style>"
        '<div class="peptide-summary">'
        '<div class="peptide-summary-title">Peptides per sample</div>'
        f'{table_html}'
        '</div>'
    )


def find_peptide_positions(report_df, fasta_df, selected_protein_id, proteotypic_only, p_value_column):
    try:
        protein_sequence = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'sequence'].iloc[0]
    except IndexError:
        raise ValueError(f"No sequence found for Protein.Ids: {selected_protein_id}")

    # filter the report DataFrame for the selected protein
    protein_report_df = report_df[report_df['Protein.Ids'].str.contains(selected_protein_id, na=False)]
    if proteotypic_only:
        protein_report_df = protein_report_df[protein_report_df['Proteotypic'] == 1]

    if protein_report_df.empty:
        raise ValueError(f"No peptides found for Protein.Ids: {selected_protein_id}")

    peptide_data = []

    for _, row in protein_report_df.iterrows():
        peptide_sequence = row['Stripped.Sequence']
        start_positions = [i for i in range(len(protein_sequence)) if protein_sequence.startswith(peptide_sequence, i)]

        if not start_positions:
            print(f"Peptide {peptide_sequence} not found in protein sequence.")
            continue

        for start_position in start_positions:
            peptide_data.append({
                'Run': row['Run'],
                'Peptide': peptide_sequence,
                'Start': start_position + 1,
                'End': start_position + len(peptide_sequence) - 1 + 1,
                'Precursor.Normalised': row['Precursor.Normalised'],
                'Precursor.Charge': row['Precursor.Charge'],
                p_value_column: row[p_value_column],
                'Proteotypic': row['Proteotypic']
            })

    peptide_positions_df = pd.DataFrame(peptide_data)
    return peptide_positions_df


def start_scheduler():
    scheduler = BackgroundScheduler()
    scheduler.add_job(
        store.purge_expired,
        'cron',
        hour=1
    )
    scheduler.start()


def clear_store_dir(store_dir):
    for filename in os.listdir(store_dir):
        file_path = os.path.join(store_dir, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


@app.route('/', methods=['GET'])
def index():
    return render_template('plot.html', app_version=__version__)


@app.route('/plot_peptides', methods=['POST'])
def plot_peptides_route():
    search_input = normalize_search_input(request.form.get('search_input'))
    proteotypic_only = request.form.get('proteotypic_checkbox') == 'true'
    sample_name_cleanup = request.form.get('sample_name_cleanup', 'none')
    sample_name_custom_pattern = request.form.get('sample_name_custom_pattern', '')
    custom_title = request.form.get('custom_title', '')
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400

    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    fasta_data = store.read(session_id, 'fasta_data')
    report_data = store.read(session_id, 'report_data')
    if not fasta_data or not report_data or not search_input:
        return jsonify({'error': 'All fields must be provided (and session must contain uploaded data).'}), 400
    fasta_df = pd.read_json(StringIO(fasta_data))
    report_df = pd.read_json(StringIO(report_data))
    report_df = apply_sample_name_cleanup(report_df, sample_name_cleanup, sample_name_custom_pattern)

    # find the P.Value column
    p_value_column = next((col for col in report_df.columns if re.match(r'^P\.Value', col)), None)

    # extract the value in brackets for hover text
    p_value_name = re.search(r'\((.*?)\)', p_value_column)
    p_value_name = p_value_name.group(1) if p_value_name else p_value_column

    # calculate global log2 intensities
    report_df['log2_intensity'] = np.log2(report_df['Precursor.Normalised'])
    global_log2_min = report_df[np.isfinite(report_df['log2_intensity'])]['log2_intensity'].min()
    global_log2_max = report_df['log2_intensity'].max()

    search_inputs = parse_search_inputs(search_input)
    if not search_inputs:
        return jsonify({'error': 'No valid search input provided.'}), 400
    search_labels = parse_search_labels(request.form.get('search_labels'), len(search_inputs))
    is_multi = len(search_inputs) > 1

    try:
        tab_items = []
        counts_by_protein = {}
        table_labels = []
        for idx, token in enumerate(search_inputs):
            label_override = search_labels[idx] if idx < len(search_labels) else ''
            table_label = label_override or token
            table_labels.append(table_label)
            counts_by_protein.setdefault(table_label, {})
            selected_protein_id = find_uniprot_id_by_gene_symbol(fasta_df, token)
            if selected_protein_id is None:
                selected_protein_id = find_uniprot_id_by_accession(fasta_df, token)
            if selected_protein_id is None:
                if not is_multi:
                    return jsonify({'error': 'No protein found for the given search input.'}), 400
                tab_items.append({
                    'label': table_label,
                    'content': f'<div class="error-message">No protein found for {table_label}.</div>'
                })
                continue

            try:
                peptide_positions_df = find_peptide_positions(
                    report_df, fasta_df, selected_protein_id, proteotypic_only, p_value_column
                )
            except Exception as e:
                if not is_multi:
                    return jsonify({'error': str(e)}), 400
                tab_items.append({
                    'label': table_label,
                    'content': f'<div class="error-message">{str(e)}</div>'
                })
                continue

            if peptide_positions_df.empty:
                message = f'No peptide positions found for {label_override or selected_protein_id}.'
                if not is_multi:
                    return jsonify({'error': message}), 400
                tab_items.append({
                    'label': table_label,
                    'content': f'<div class="error-message">{message}</div>'
                })
                continue

            if is_multi:
                plot_title = label_override or display_label
            else:
                plot_title = label_override or custom_title
            plot_html = plot_peptides(
                peptide_positions_df,
                fasta_df,
                selected_protein_id,
                global_log2_min,
                global_log2_max,
                p_value_column,
                p_value_name,
                plot_title
            )
            protein_label = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'gene_symbol'].iloc[0]
            display_label = protein_label if isinstance(protein_label, str) and protein_label.strip() else selected_protein_id
            display_label = label_override or display_label
            tab_items.append({
                'label': display_label,
                'content': plot_html
            })
            unique_peptides = peptide_positions_df.drop_duplicates(subset=['Run', 'Peptide'])
            counts_by_protein[table_label] = unique_peptides.groupby('Run')['Peptide'].size().to_dict()

        if is_multi:
            tabs_html = render_tabs('pepmap-peptides-tabs', tab_items)
            all_runs = sorted(report_df['Run'].astype(str).unique())
            summary_html = build_peptide_summary_table(counts_by_protein, all_runs, column_order=table_labels)
            return f'{tabs_html}{summary_html}'

        return tab_items[0]['content']

    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/plot_features', methods=['POST'])
def plot_features_route():
    search_input = normalize_search_input(request.form.get('search_input'))
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400

    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    fasta_data = store.read(session_id, 'fasta_data')
    if not fasta_data or not search_input:
        return jsonify({'error': 'All fields must be provided (and session must contain uploaded data).'}), 400
    fasta_df = pd.read_json(StringIO(fasta_data))
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

    try:
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
                    f'<div class="error-message">No protein found for {token}.</div>'
                    '</div>'
                )
                continue

            features_html = plot_features(fasta_df, selected_protein_id, custom_features_df, custom_features_label)
            if isinstance(features_html, tuple):
                message = features_html[0].json.get('error', 'Failed to load features plot.')
                if not is_multi:
                    return jsonify({'error': message}), 400
                feature_panels.append(
                    f'<div class="feature-panel" data-index="{idx}">'
                    f'<div class="error-message">{message}</div>'
                    '</div>'
                )
                continue

            feature_panels.append(
                f'<div class="feature-panel" data-index="{idx}">{features_html}</div>'
            )

        if is_multi:
            return ''.join(feature_panels)

        return feature_panels[0]

    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/upload', methods=['POST'])
def upload_files():
    report_file = request.files.get('report_file')
    fasta_file = request.files.get('fasta_file')
    organism = normalize_organism(request.form.get('organism'))
    custom_features_file = request.files.get('custom_features_file')
    custom_features_label = normalize_search_input(request.form.get('custom_features_label'))
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400

    if report_file and fasta_file:
        custom_features_uploaded = False
        if not session_id:
            session_id = store.create(meta={'user_agent': request.headers.get('User-Agent', '')})
        elif store.get_payload(session_id) is None:
            store.create_with_id(session_id, meta={'user_agent': request.headers.get('User-Agent', '')})

        update_kwargs = {
            'fasta_data': parse_fasta(fasta_file.stream, organism or '').to_json(),
            'report_data': parse_report_tsv(report_file.stream).to_json(),
            'organism': organism,
            'report_filename': report_file.filename,
            'fasta_filename': fasta_file.filename,
            'custom_features_filename': custom_features_file.filename if custom_features_file else None
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
        return jsonify({
            'message': 'Files uploaded successfully',
            'custom_features_uploaded': custom_features_uploaded,
            'custom_features_label': custom_features_label,
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
        return jsonify({
            'message': 'Custom features uploaded successfully',
            'custom_features_uploaded': True,
            'custom_features_label': custom_features_label,
            'session_id': session_id
        }), 200
    return jsonify({'error': 'Missing files'}), 400


@app.route('/autocomplete', methods=['GET'])
def autocomplete():
    query = normalize_search_input(request.args.get('query'))
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'suggestions': []}), 200
    if not query:
        return jsonify({'suggestions': []}), 200

    if not session_id:
        return jsonify({'suggestions': []}), 200

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


@app.route('/session_info', methods=['GET'])
def session_info():
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400
    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400

    payload = store.get_payload(session_id)
    if payload is None:
        return jsonify({'error': 'Unknown or expired session_id.'}), 404

    data = payload.get('data', {})
    return jsonify({
        'session_id': session_id,
        'report_filename': data.get('report_filename'),
        'fasta_filename': data.get('fasta_filename'),
        'custom_features_filename': data.get('custom_features_filename'),
        'custom_features_label': data.get('custom_features_label'),
        'organism': data.get('organism'),
        'has_custom_features': bool(data.get('custom_features_data'))
    }), 200


@app.route('/flush', methods=['POST'])
def flush_session():
    session_id, session_error = get_request_session_id()
    if session_error:
        return jsonify({'error': session_error}), 400
    if not session_id:
        return jsonify({'error': 'session_id is required.'}), 400
    store.delete(session_id)
    return jsonify({'message': 'Session cleared'}), 200


if __name__ == '__main__':
    clear_store_dir(PEPMAP_STORE_DIR)
    start_scheduler()
    app.run(port=7007, debug=False)
