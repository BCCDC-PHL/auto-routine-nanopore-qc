import csv
import glob
import logging

from pathlib import Path

from auto_routine_nanopore_qc.model import Config, CustomJSONEncoder, Run, Pipeline

log = logging.getLogger(__name__)

def find_samplesheet(run: Run) -> Path|None:
    """
    Given a sequencing run, find its sample_sheet.csv file.
    Returns None if not found.
    """
    samplesheet_path = None
    samplesheet_glob = str(run.path / "sample_sheet*.csv")
    samplesheets = glob.glob(samplesheet_glob)
    if len(samplesheets) > 0:
        samplesheet_path = Path(samplesheets[0])

    return samplesheet_path


def parse_samplesheet(samplesheet_path: Path) -> list[dict]:
    """
    Parse a sample_sheet.csv file, breaking the 'alias' field
    into sample_id and project_id.

    Return list of dicts with keys 'sample_id', 'project_id', 'barcode'
    """
    parsed_samplesheet = []
    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            parsed_row = {}
            parsed_row['barcode'] = row['barcode']
            alias_split = row['alias'].split('_', 1)
            sample_id = alias_split[0]
            parsed_row['sample_id'] = sample_id
            project_id = None
            if len(alias_split) > 1:
                project_id = alias_split[1]
            parsed_row['project_id'] = project_id
            parsed_samplesheet.append(parsed_row)

    return parsed_samplesheet
