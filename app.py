import os
import re
import numpy as np
import shutil
import pandas as pd
import plotly.io as pio
import plotly.graph_objs as go
import requests
from io import StringIO
from flask_caching import Cache
from flask import Flask, render_template, request, jsonify, session
from apscheduler.schedulers.background import BackgroundScheduler
from collections import OrderedDict
from flask_session import Session
from datetime import timedelta

app = Flask(__name__)
app.config['CACHE_TYPE'] = 'SimpleCache'
app.config['SECRET_KEY'] = os.urandom(24)
app.config['SESSION_TYPE'] = 'filesystem'
app.config['SESSION_PERMANENT'] = False
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(minutes=30)
app.config['SESSION_FILE_DIR'] = 'flask_session'
cache = Cache(app)
Session(app)


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


def plot_peptides(peptide_positions_df, fasta_df, selected_protein_id, global_log2_min, global_log2_max, p_value_column, p_value_name):
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

    layout = go.Layout(
        title=f'Peptide Mapping for {selected_protein_name} ({selected_protein_id})',
        xaxis=dict(
            range=[1, protein_length],
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext,
            tickangle=0,
            tickfont=dict(size=9),
            ticks='outside',
            fixedrange=True
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
        margin=dict(l=250, r=100, t=40, b=0),
        shapes=shapes,
        height=current_y + 50
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


def plot_features(fasta_df, selected_protein_id):
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


def clear_filesystem_sessions(session_dir):
    for filename in os.listdir(session_dir):
        file_path = os.path.join(session_dir, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def start_scheduler():
    scheduler = BackgroundScheduler()
    # trigger 'clear_filesystem_sessions' at 1 AM every day
    scheduler.add_job(
        clear_filesystem_sessions,
        'cron',
        hour=1,
        args=[app.config['SESSION_FILE_DIR']]
    )
    scheduler.start()


@app.before_request
def make_session_permanent():
    session.permanent = False


@app.route('/', methods=['GET'])
def index():
    return render_template('plot.html')


@app.route('/plot_peptides', methods=['POST'])
def plot_peptides_route():
    search_input = request.form.get('search_input')
    proteotypic_only = request.form.get('proteotypic_checkbox') == 'true'

    if 'fasta_data' in session and 'report_data' in session and search_input is not None:
        fasta_df = pd.read_json(StringIO(session['fasta_data']))
        report_df = pd.read_json(StringIO(session['report_data']))
        selected_protein_id = fasta_df.loc[fasta_df['uniprot_id'] == search_input, 'uniprot_id']

        # find the P.Value column
        p_value_column = next((col for col in report_df.columns if re.match(r'^P\.Value', col)), None)

        # extract the value in brackets for hover text
        p_value_name = re.search(r'\((.*?)\)', p_value_column)
        p_value_name = p_value_name.group(1) if p_value_name else p_value_column

        # calculate global log2 intensities
        report_df['log2_intensity'] = np.log2(report_df['Precursor.Normalised'])
        global_log2_min = report_df[np.isfinite(report_df['log2_intensity'])]['log2_intensity'].min()
        global_log2_max = report_df['log2_intensity'].max()

        if len(selected_protein_id) == 0:
            selected_protein_id = fasta_df.loc[fasta_df['gene_symbol'] == search_input, 'uniprot_id']
            if len(selected_protein_id) == 0:
                return jsonify({'error': 'No protein found for the given search input.'}), 400

        selected_protein_id = selected_protein_id.iloc[0]

    else:
        return jsonify({'error': 'All fields must be provided.'}), 400

    try:
        # find peptide positions
        peptide_positions_df = find_peptide_positions(report_df, fasta_df, selected_protein_id, proteotypic_only, p_value_column)
        if peptide_positions_df.empty:
            return jsonify({'error': 'No peptide positions found.'}), 400

        return plot_peptides(peptide_positions_df, fasta_df, selected_protein_id, global_log2_min, global_log2_max, p_value_column, p_value_name)

    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/plot_features', methods=['POST'])
def plot_features_route():
    search_input = request.form.get('search_input')

    if 'fasta_data' in session and search_input is not None:
        fasta_df = pd.read_json(StringIO(session['fasta_data']))

        selected_protein_id = fasta_df.loc[fasta_df['uniprot_id'] == search_input, 'uniprot_id']
        if len(selected_protein_id) == 0:
            selected_protein_id = fasta_df.loc[fasta_df['gene_symbol'] == search_input, 'uniprot_id']
            if len(selected_protein_id) == 0:
                return jsonify({'error': 'No protein found for the given search input.'}), 400

        selected_protein_id = selected_protein_id.iloc[0]
    else:
        return jsonify({'error': 'All fields must be provided.'}), 400

    try:
        return plot_features(fasta_df, selected_protein_id)

    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/upload', methods=['POST'])
def upload_files():
    report_file = request.files.get('report_file')
    fasta_file = request.files.get('fasta_file')
    organism = request.form.get('organism')

    if report_file and fasta_file:
        session['fasta_data'] = parse_fasta(fasta_file.stream, organism).to_json()
        session['report_data'] = parse_report_tsv(report_file.stream).to_json()
        return jsonify({'message': 'Files uploaded successfully'}), 200
    else:
        return jsonify({'error': 'Missing files'}), 400


@app.route('/flush', methods=['POST'])
def flush_session():
    session.clear()
    return jsonify({'message': 'Session cleared'}), 200


if __name__ == '__main__':
    clear_filesystem_sessions(app.config['SESSION_FILE_DIR'])
    start_scheduler()
    app.run(port=7007, debug=True)
