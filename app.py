from flask import Flask, render_template, request
import pandas as pd
import plotly.io as pio
import plotly.graph_objs as go
import requests
from flask import jsonify
from collections import OrderedDict

app = Flask(__name__)


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
        ticktext.append(str(protein_length))

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


def parse_fasta(fasta_stream):
    entries = []
    sequence = ''
    uniprot_id = ''
    gene_symbol = ''
    for line in fasta_stream:
        line = line.decode('utf-8').strip()
        if line.startswith('>'):
            if sequence:
                entries.append({'uniprot_id': uniprot_id, 'gene_symbol': gene_symbol, 'sequence': sequence})
            sequence = ''
            parts = line[1:].split('|')
            uniprot_id = parts[1]
            gene_symbol = parts[2].split(' ')[0]
        else:
            sequence += line
    if sequence:
        entries.append({'uniprot_id': uniprot_id, 'gene_symbol': gene_symbol, 'sequence': sequence})

    fasta_df = pd.DataFrame(entries)
    return fasta_df


def parse_report_tsv(tsv_stream):
    report_df = pd.read_csv(tsv_stream, delimiter='\t')

    report_df = report_df[['Run', 'Protein.Ids', 'Precursor.Normalised', 'Stripped.Sequence',
                           'Precursor.Charge', 'Q.Value']]

    return report_df


def plot_peptides(peptide_positions_df, fasta_df, selected_protein_id):
    # get the protein sequence for the selected protein ID
    protein_sequence = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'sequence'].iloc[0]
    protein_length = len(protein_sequence)

    # normalize the values to use them for color mapping
    pep_int_max = peptide_positions_df['Precursor.Normalised'].max()
    peptide_positions_df['color_intensity'] = (
            peptide_positions_df['Precursor.Normalised'] / pep_int_max
    ) if pep_int_max else 0

    # create a sorted list of unique runs for the y-axis categories
    unique_runs = sorted(peptide_positions_df['Run'].unique())
    run_categories = {run: i for i, run in enumerate(unique_runs, start=1)}

    # assign y positions for each peptide based on its run category
    peptide_positions_df['y_pos'] = peptide_positions_df['Run'].map(run_categories)

    # sort peptides within each run by their start position
    peptide_positions_df.sort_values(['Run', 'Start'], ascending=[True, True], inplace=True)

    # initialize an empty dictionary to keep track of the last end position for each run
    last_end_positions = {run: 0 for run in unique_runs}

    peptide_bar_height = 0.5
    peptide_bar_margin = 0.2
    peptide_bar_line_width = 0.7
    run_gap = peptide_bar_height  # initial gap before the first run
    run_offset = 0  # initial offset for the first run
    run_height = 40
    run_label_size = 14

    global_max_sub_row_count = 0
    traces = []
    temp_traces = []
    shapes = []

    for run, group in peptide_positions_df.groupby('Run'):
        group = group.reset_index()
        max_offset_within_run = 1
        for i, (idx, row) in enumerate(group.iterrows()):
            # determine if the current peptide overlaps with the previous one
            overlap_offset = 0
            if row['Start'] < last_end_positions[row['Run']]:
                overlap_offset = max_offset_within_run + peptide_bar_height + peptide_bar_margin

            # update the last end position for the run
            last_end_positions[row['Run']] = row['End']
            max_offset_within_run = overlap_offset

            color = f'rgba(255,{255 - row["color_intensity"] * 255},0,0.8)'
            trace = go.Bar(
                x=[row['End'] - row['Start']],
                y=[row['y_pos'] + overlap_offset + run_offset],
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
                          f'<br>Intensity: {row["Precursor.Normalised"]:.2f}'
                          f'<br>Charge: {row["Precursor.Charge"]}'
                          f'<br>Q Value: {row["Q.Value"]:.7f}',
                showlegend=False
            )
            traces.append(trace)
            temp_traces.append(trace)

        sub_row_count = [trace['base'] for trace in temp_traces]
        # clear the temporary traces list
        temp_traces = []
        # count all identical start positions of list
        sub_row_count_dict = {i: sub_row_count.count(i) for i in sub_row_count}
        # get the maximum count of identical start positions
        max_sub_row_count = max(sub_row_count_dict.values())

        if max_sub_row_count > global_max_sub_row_count:
            global_max_sub_row_count = max_sub_row_count

        # line after each sample
        shapes.append({
            'type': 'line',
            'x0': 0,
            'y0': run_gap,
            'x1': protein_length,
            'y1': run_gap,
            'line': {
                'color': 'black',
                'width': 1,
            },
        })

        # get the run gap position after the last peptide in the run
        run_gap = max([trace['y'][0] for trace in traces]) + peptide_bar_height

        # increase the group offset for the next run
        if global_max_sub_row_count > 1:
            run_offset += 2 * (peptide_bar_height + peptide_bar_margin)
        else:
            run_offset = 0

    tickvals, ticktext = generate_dynamic_ticks(protein_length)

    layout = go.Layout(
        title=f'Peptide Mapping for {selected_protein_id}',
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
            tickvals=[shape["y0"] + peptide_bar_height for shape in shapes],
            ticktext=unique_runs,
            tickfont=dict(size=run_label_size),
            fixedrange=True
        ),
        barmode='stack',
        showlegend=False,
        plot_bgcolor='white',
        margin=dict(l=250, r=20, t=40, b=0),
        shapes=shapes,
        # calculate max peptide per row and multiply it with the total run count
        height=(global_max_sub_row_count * run_height) * len(unique_runs)
        if global_max_sub_row_count > 1 else (run_height * 2) * len(unique_runs),
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
            # create a bar for each variant feature
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
                          f"<br>SAAV: {protein_sequence[int(feature_position)-1]}->{feature['alternativeSequence']}"
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
        margin=dict(l=250, r=20, t=0, b=0),
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


@app.route('/', methods=['GET'])
def index():
    return render_template('plot.html')


def find_peptide_positions(report_df, fasta_df, selected_protein_id):
    try:
        protein_sequence = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'sequence'].iloc[0]
    except IndexError:
        raise ValueError(f"No sequence found for Protein.Ids: {selected_protein_id}")

    # filter the report DataFrame for the selected protein
    protein_report_df = report_df[report_df['Protein.Ids'] == selected_protein_id]
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
                'Start': start_position,
                'End': start_position + len(peptide_sequence) - 1,
                'Precursor.Normalised': row['Precursor.Normalised'],
                'Precursor.Charge': row['Precursor.Charge'],
                'Q.Value': row['Q.Value']
            })

    peptide_positions_df = pd.DataFrame(peptide_data)
    return peptide_positions_df


@app.route('/plot_peptides', methods=['POST'])
def plot_peptides_route():
    report_file = request.files.get('report_file')
    fasta_file = request.files.get('fasta_file')
    selected_protein_id = request.form.get('uniprot_id')

    if not report_file or not fasta_file or not selected_protein_id:
        return jsonify({'error': 'All fields must be provided.'}), 400

    try:
        # parse the FASTA file into a DataFrame
        fasta_df = parse_fasta(fasta_file.stream)

        # parse the report.tsv file
        report_df = parse_report_tsv(report_file.stream)

        # find peptide positions
        peptide_positions_df = find_peptide_positions(report_df, fasta_df, selected_protein_id)
        if peptide_positions_df.empty:
            return jsonify({'error': 'No peptide positions found.'}), 400

        return plot_peptides(peptide_positions_df, fasta_df, selected_protein_id)

    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/plot_features', methods=['POST'])
def plot_features_route():
    fasta_file = request.files.get('fasta_file')
    selected_protein_id = request.form.get('uniprot_id')

    if not fasta_file or not selected_protein_id:
        return jsonify({'error': 'All fields must be provided.'}), 400

    try:
        # parse the FASTA file into a DataFrame
        fasta_df = parse_fasta(fasta_file.stream)

        return plot_features(fasta_df, selected_protein_id)

    except Exception as e:
        return jsonify({'error': str(e)}), 400


if __name__ == '__main__':
    app.run(port=7007, debug=True)
