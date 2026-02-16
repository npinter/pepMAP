from .diann import REQUIRED_COLUMNS as DIANN_REQUIRED, parse_diann
from .fragpipe import REQUIRED_COLUMNS as FRAGPIPE_REQUIRED, parse_fragpipe
from .maxquant import REQUIRED_COLUMNS as MAXQUANT_REQUIRED, parse_maxquant
from .spectronaut import REQUIRED_COLUMNS as SPECTRONAUT_REQUIRED, parse_spectronaut
from .utils import read_report_columns

PARSERS = [
    ('DIA-NN', DIANN_REQUIRED, parse_diann),
    ('FragPipe', FRAGPIPE_REQUIRED, parse_fragpipe),
    ('MaxQuant', MAXQUANT_REQUIRED, parse_maxquant),
    ('Spectronaut', SPECTRONAUT_REQUIRED, parse_spectronaut)
]


def get_missing_columns(header, required_columns):
    header_set = set(header)
    missing = []
    for entry in required_columns:
        if isinstance(entry, (list, tuple, set)):
            if not any(candidate in header_set for candidate in entry):
                missing.append(' or '.join(entry))
        elif entry not in header_set:
            missing.append(entry)
    return missing

def detect_parser(tsv_stream, filename=None):
    header = read_report_columns(tsv_stream, filename)

    missing_by_parser = {}
    for name, required_columns, parser in PARSERS:
        missing = get_missing_columns(header, required_columns)
        missing_by_parser[name] = missing
        if not missing:
            return name, parser

    details = '; '.join(
        f"{name}: missing {', '.join(missing) if missing else 'none'}"
        for name, missing in missing_by_parser.items()
    )
    raise ValueError(
        "Unable to determine file type. Provide a supported DIA-NN, FragPipe, MaxQuant, or Spectronaut report. "
        f"Missing columns -> {details}."
    )
