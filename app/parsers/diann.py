import pandas as pd

from .utils import read_report_table

REQUIRED_COLUMNS = [
    'Run',
    'Protein.Ids',
    'Precursor.Normalised',
    'Stripped.Sequence',
    'Precursor.Charge',
    'Q.Value',
    'Proteotypic'
]


def parse_diann(tsv_stream, filename=None):
    report_df = read_report_table(tsv_stream, filename, delimiter='\t')
    ep_source = 'Q.Value'

    ep_column = f"EP ({ep_source})"
    report_df = report_df.rename(columns={
        'Precursor.Normalised': 'Intensity',
        'Stripped.Sequence': 'Sequence',
        'Precursor.Charge': 'Charge',
        ep_source: ep_column
    })
    report_df = report_df[[
        'Run',
        'Protein.Ids',
        'Intensity',
        'Sequence',
        'Charge',
        'Proteotypic',
        ep_column
    ]]
    report_df = report_df[report_df['Intensity'] != 0]
    return report_df, 'DIA-NN.png'
