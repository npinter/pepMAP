import os

import pandas as pd


def is_parquet_filename(filename):
    if not filename:
        return False
    _, ext = os.path.splitext(str(filename))
    return ext.lower() in {'.parquet', '.pq'}


def read_report_table(report_stream, filename=None, delimiter='\t'):
    if is_parquet_filename(filename):
        try:
            report_stream.seek(0)
            return pd.read_parquet(report_stream)
        except Exception as exc:
            raise ValueError(
                "Unable to read Parquet report. Install a Parquet engine such as pyarrow."
            ) from exc
    report_stream.seek(0)
    return pd.read_csv(report_stream, delimiter=delimiter)


def read_report_columns(report_stream, filename=None, delimiter='\t'):
    if is_parquet_filename(filename):
        try:
            import pyarrow.parquet as pq
            report_stream.seek(0)
            parquet_file = pq.ParquetFile(report_stream)
            columns = parquet_file.schema.names
            report_stream.seek(0)
            return columns
        except Exception:
            report_stream.seek(0)
            df = pd.read_parquet(report_stream)
            report_stream.seek(0)
            return list(df.columns)

    header = report_stream.readline().decode('utf-8').strip().split(delimiter)
    report_stream.seek(0)
    return header
