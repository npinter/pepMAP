import pandas as pd

from .utils import read_report_table

REQUIRED_COLUMNS = [
    'Spectrum',
    'Protein ID',
    'Intensity',
    'Peptide',
    'Charge',
    'Is Unique',
    ['Q.value', 'PeptideProphet Probability']
]


def parse_fragpipe(tsv_stream, filename=None):
    report_df = read_report_table(tsv_stream, filename, delimiter='\t')
    if 'Q.value' in report_df.columns:
        # FragPipe v23+ includes Q.value directly.
        ep_source = 'Q.value'
        ep_transform = None
    elif 'PeptideProphet Probability' in report_df.columns:
        # FragPipe <= v23: use PeptideProphet Probability as a PEP-like proxy
        ep_source = 'PeptideProphet Probability'
        ep_transform = lambda series: 1 - pd.to_numeric(series, errors='coerce')
    else:
        ep_source = 'Q.value'
        ep_transform = None
    ep_column = f"EP ({ep_source})"

    report_df = report_df.rename(columns={
        'Spectrum': 'Run',
        'Protein ID': 'Protein.Ids',
        'Intensity': 'Intensity',
        'Peptide': 'Sequence',
        'Charge': 'Charge',
        'Is Unique': 'Proteotypic',
        ep_source: ep_column
    })
    if ep_transform is not None:
        report_df[ep_column] = ep_transform(report_df[ep_column])

    report_df['Run'] = report_df['Run'].str.rsplit('.', n=3).str[0]

    report_df = report_df.groupby(['Run', 'Sequence', 'Charge'], as_index=False).agg({
        'Intensity': 'sum',
        ep_column: 'min',
        'Proteotypic': 'first',
        'Protein.Ids': 'first',
        'Mapped Proteins': 'first'
    })

    def extract_protein_ids(mapped_proteins):
        if not isinstance(mapped_proteins, str):
            return ''
        return ';'.join([p.split('|')[1] for p in mapped_proteins.split(',') if '|' in p])

    report_df['Protein.Ids'] = report_df.apply(
        lambda row: (
            f"{row['Protein.Ids']};{extract_protein_ids(row['Mapped Proteins'])}"
            if pd.notna(row['Mapped Proteins']) else row['Protein.Ids']
        ),
        axis=1
    )

    report_df = report_df[[
        'Run',
        'Protein.Ids',
        'Intensity',
        'Sequence',
        'Charge',
        'Proteotypic',
        ep_column,
    ]]

    report_df = report_df[report_df['Intensity'] != 0]
    return report_df, 'FragPipe.png'
