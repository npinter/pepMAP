import json
import re

import pandas as pd


def _extract_gene_symbol(header, organism):
    match = re.search(r'\bGN=([^\s]+)', header)
    if match:
        return match.group(1)
    parts = header.split('|')
    try:
        candidate = parts[2].split(' ')[0]
        if organism:
            return candidate.split(f"_{organism}")[0]
        return candidate
    except IndexError:
        return header


def parse_fasta(fasta_stream, organism):
    entries = []
    sequence = ''
    uniprot_id = ''
    gene_symbol = ''
    for line in fasta_stream:
        line = line.decode('utf-8').strip()
        if line.startswith('>rev_sp') or line.startswith('>contam_sp') or line.startswith('>rev_') or line.startswith('>contam_'):
            continue
        if line.startswith('>'):
            if sequence:
                entries.append({'uniprot_id': uniprot_id, 'gene_symbol': gene_symbol, 'sequence': sequence})
            sequence = ''
            header = line[1:]
            parts = header.split('|')
            try:
                uniprot_id = parts[1]
            except IndexError:
                uniprot_id = header
            gene_symbol = _extract_gene_symbol(header, organism)
        else:
            sequence += line
    if sequence:
        entries.append({'uniprot_id': uniprot_id, 'gene_symbol': gene_symbol, 'sequence': sequence})
    return pd.DataFrame(entries)


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


def normalize_organism(organism, custom_organism=None):
    if organism is None:
        return None
    organism_value = organism.strip().upper()
    if organism_value == 'CUSTOM':
        custom_value = (custom_organism or '').strip()
        return custom_value.upper() if custom_value else None
    if organism_value in {'HUMAN', 'MOUSE'}:
        return organism_value
    return None


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


def _get_sample_name_cleanup_regex(cleanup_mode, custom_split_pattern):
    if cleanup_mode == 'split_underscore':
        pattern = r'^([^_]+)'
    elif cleanup_mode == 'custom':
        pattern = (custom_split_pattern or '').strip()
    else:
        return None, None
    if not pattern:
        return None, None
    try:
        return re.compile(pattern), None
    except re.error as exc:
        return None, str(exc)


def _apply_sample_name_regex(value, regex):
    text = str(value)
    if regex is None:
        return text
    match = regex.search(text)
    if not match:
        return text
    if match.lastindex:
        group_value = match.group(1)
        return group_value if group_value is not None else match.group(0)
    return match.group(0)


def apply_sample_name_cleanup_value(value, cleanup_mode, custom_split_pattern):
    regex, _ = _get_sample_name_cleanup_regex(cleanup_mode, custom_split_pattern)
    return _apply_sample_name_regex(value, regex)


def get_sample_name_cleanup_preview(value, cleanup_mode, custom_split_pattern):
    regex, error = _get_sample_name_cleanup_regex(cleanup_mode, custom_split_pattern)
    cleaned = _apply_sample_name_regex(value, regex)
    return cleaned, error


def apply_sample_name_cleanup(report_df, cleanup_mode, custom_split_pattern):
    if cleanup_mode not in {'split_underscore', 'custom'}:
        return report_df
    cleaned_df = report_df.copy()
    if 'Run' not in cleaned_df.columns:
        return cleaned_df
    regex, _ = _get_sample_name_cleanup_regex(cleanup_mode, custom_split_pattern)
    if regex is None:
        return cleaned_df
    cleaned_df['Run'] = cleaned_df['Run'].astype(str).apply(lambda value: _apply_sample_name_regex(value, regex))
    return cleaned_df


def apply_charge_state_mode(report_df, mode, p_value_column):
    if mode == 'unique':
        idx = report_df.groupby(['Run', 'Sequence'])[p_value_column].idxmin()
        result = report_df.loc[idx]
        return result

    if mode == 'overlap':
        group_cols = ['Run', 'Protein.Ids'] if 'Protein.Ids' in report_df.columns else ['Run']
        kept_indices = []
        for _, group in report_df.groupby(group_cols, sort=False):
            if group.empty:
                continue
            idx = group.groupby('Sequence', sort=False)[p_value_column].idxmin()
            group = group.loc[idx].copy()
            group = group.sort_values(p_value_column, kind='mergesort')
            kept_sequences = []
            for row_idx, seq in group['Sequence'].astype(str).items():
                if any(seq in kept or kept in seq for kept in kept_sequences):
                    continue
                kept_indices.append(row_idx)
                kept_sequences.append(seq)

        result = report_df.loc[kept_indices]
        return result

    return report_df


def apply_q_value_cutoff(report_df, cutoff, p_value_column):
    if cutoff == 0:
        return report_df
    if p_value_column not in report_df.columns:
        return report_df
    return report_df[report_df[p_value_column] <= cutoff]


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
        raise ValueError('No valid custom features found in the uploaded file.')
    return parsed_df


def find_peptide_positions(
    report_df,
    fasta_df,
    selected_protein_id,
    proteotypic_only,
    p_value_column,
    protein_index=None
):
    try:
        protein_sequence = fasta_df.loc[fasta_df['uniprot_id'] == selected_protein_id, 'sequence'].iloc[0]
    except IndexError:
        raise ValueError(f"No sequence found for Protein.Ids: {selected_protein_id}")

    protein_report_df = None
    if isinstance(protein_index, dict):
        indices = protein_index.get(selected_protein_id, [])
        if indices:
            protein_report_df = report_df.loc[report_df.index.intersection(indices)]
        else:
            protein_report_df = report_df.iloc[0:0]
    if protein_report_df is None:
        protein_report_df = report_df[report_df['Protein.Ids'].str.contains(selected_protein_id, na=False, regex=False)]
    if proteotypic_only:
        protein_report_df = protein_report_df[protein_report_df['Proteotypic'] == 1]

    if protein_report_df.empty:
        raise ValueError(f"No peptides found for Protein.Ids: {selected_protein_id}")

    peptide_data = []
    positions_cache = {}
    columns = ['Run', 'Sequence', 'Intensity', 'Charge', p_value_column, 'Proteotypic']
    for run, peptide_sequence, intensity, charge, p_value, proteotypic in (
        protein_report_df[columns].itertuples(index=False, name=None)
    ):
        start_positions = positions_cache.get(peptide_sequence)
        if start_positions is None:
            start_positions = []
            start_idx = protein_sequence.find(peptide_sequence)
            while start_idx != -1:
                start_positions.append(start_idx)
                start_idx = protein_sequence.find(peptide_sequence, start_idx + 1)
            positions_cache[peptide_sequence] = start_positions
        if not start_positions:
            continue
        peptide_length = len(peptide_sequence)
        for start_position in start_positions:
            peptide_data.append({
                'Run': run,
                'Peptide': peptide_sequence,
                'Start': start_position + 1,
                'End': start_position + peptide_length,
                'Intensity': intensity,
                'Charge': charge,
                p_value_column: p_value,
                'Proteotypic': proteotypic
            })

    result = pd.DataFrame(peptide_data)
    return result
