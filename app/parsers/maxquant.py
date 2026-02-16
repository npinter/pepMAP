import pandas as pd

from .utils import read_report_table

REQUIRED_COLUMNS = [
    'Experiment',
    'Proteins',
    'Intensity',
    'Sequence',
    'Charge',
    'PIF'
]


def parse_maxquant(tsv_stream, filename=None):
    report_df = read_report_table(tsv_stream, filename, delimiter='\t')

    report_df['Proteotypic'] = ~report_df['Proteins'].str.contains(';', na=False)

    report_df = report_df.rename(columns={
        'Experiment': 'Run',
        'Proteins': 'Protein.Ids',
        'Intensity': 'Intensity',
        'Charge': 'Charge',
        'PIF': 'EP (PIF)'
    })

    report_df = report_df[[
        'Run',
        'Protein.Ids',
        'Intensity',
        'Sequence',
        'Charge',
        'Proteotypic',
        'EP (PIF)'
    ]]

    report_df = report_df[report_df['Intensity'] != 0]
    return report_df, 'MaxQuant.png'
