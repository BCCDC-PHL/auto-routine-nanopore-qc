import csv
import json
import logging

from pathlib import Path

log = logging.getLogger(__name__)

def parse_samplesheet(samplesheet_path: Path) -> list[dict]:
    """
    """
    parsed_samplesheet = []
    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            parsed_samplesheet.append(row)
            
    return parsed_samplesheet
