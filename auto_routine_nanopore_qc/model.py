import json
from datetime import datetime

from dataclasses import dataclass, field, fields
from enum import StrEnum
from pathlib import Path
from typing import Optional


class CustomJSONEncoder(json.JSONEncoder):
    """
    Helper class to assist with JSON serialization. Pass to json.dumps as cls=CustomJSONEncoder
    """
    def default(self, o):
        if isinstance(o, Path):
            return str(o)
        if isinstance(o, datetime):
            return o.isoformat()
        return super().default(o)


class InstrumentType(StrEnum):
    """
    Instrument type ('gridion', 'promethion' or 'unknown')
    """
    gridion   = "gridion"
    promethion = "promethion"
    unknown = "unknown"


@dataclass
class Config:
    """
    Main application config.
    """
    analysis_output_dir: Path
    analysis_work_dir: Path
    run_parent_dirs: list[Path] = field(default_factory=list)
    notification: dict = field(default_factory=dict)
    scan_interval_seconds: float = 3600
    excluded_runs_list: Optional[Path] = None
    excluded_runs: list[str] = field(default_factory=list)
    known_species_list: Optional[Path] = None
    known_species: list[dict] = field(default_factory=list)
    projects_definition_file: Optional[Path] = None
    projects: list[dict] = field(default_factory=list)

    def __post_init__(self):
        self.scan_interval_seconds = float(self.scan_interval_seconds)
        self.run_parent_dirs = [Path(p) for p in self.run_parent_dirs]
        self.analysis_output_dir = Path(self.analysis_output_dir)
        self.analysis_work_dir = Path(self.analysis_work_dir)

    @classmethod
    def from_dict(cls, data: dict):
        # Get all valid field names for this dataclass
        valid_fields = {f.name for f in fields(cls)}
        # Filter the input dictionary
        filtered_data = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered_data)


@dataclass
class Run:
    """
    A sequencing run directory that is ready to be QC checked.
    """
    sequencing_run_id: str
    path: Path
    instrument_type: InstrumentType
