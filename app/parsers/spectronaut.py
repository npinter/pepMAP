import pandas as pd

from .utils import read_report_table

REQUIRED_COLUMNS = [
    'R.FileName',
    'PG.ProteinGroups',
    'EG.TotalQuantity (Settings)',
    'PEP.StrippedSequence',
    'FG.Charge',
    'EG.Qvalue',
    'PEP.IsProteotypic'
]


def parse_spectronaut(tsv_stream, filename=None):
    report_df = read_report_table(tsv_stream, filename, delimiter='\t')
    ep_source = 'EG.Qvalue'
    ep_column = f"EP ({ep_source})"

    report_df = report_df.rename(columns={
        'R.FileName': 'Run',
        'PG.ProteinGroups': 'Protein.Ids',
        'EG.TotalQuantity (Settings)': 'Intensity',
        'PEP.StrippedSequence': 'Sequence',
        'FG.Charge': 'Charge',
        'PEP.IsProteotypic': 'Proteotypic',
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
    return report_df, 'Spectronaut.png'
