import datetime
import json
import logging
import os
import re
import shutil
import subprocess
import uuid

from typing import Iterator, Optional
from dataclasses import asdict
from pathlib import Path

import auto_routine_nanopore_qc.instrument as instrument
import auto_routine_nanopore_qc.pre_analysis as pre_analysis
import auto_routine_nanopore_qc.analysis as analysis
import auto_routine_nanopore_qc.post_analysis as post_analysis

from auto_routine_nanopore_qc.model import Config, CustomJSONEncoder, InstrumentType, Run


log = logging.getLogger(__name__)


def is_readable(run_dir: Path) -> bool:
    """
    Check whether the files needed for the QC analysis are readable.
    Newly-uploaded runs are only readable by the uploading user until permissions are updated,
    and the run directory may be readable before the files inside it are.

    :param run_dir: Path to the run directory.
    :return: True if the run directory and other required files are readable
    """
    if not os.access(run_dir, os.R_OK | os.X_OK):
        return False
    # TODO: determine exactly what we need to be able to read, to proceed.
    # SampleSheet? combine_fastq_complete.json file? 
    paths = []

    return all(os.access(path, os.R_OK) for path in paths)


def find_run_dirs(config: Config):
    """
    Find sequencing run directories under the 'run_parent_dirs' listed in the config
    that are ready to be QC checked.

    :param config: Application config.
    :return: Runs that are ready to be QC checked
    """
    for run_parent_dir in config.run_parent_dirs:
        subdirs = []
        try:
            subdirs = list(os.scandir(run_parent_dir))
        except OSError as e:
            # Parent dirs may be on network storage that's temporarily unavailable.
            # Skip it for this scan, and try again on the next one.
            log.error({
                "event_type": "failed_to_scan_run_parent_dir",
                "run_parent_dir": str(run_parent_dir),
                "error": str(e)
            })
            continue

        for subdir in subdirs:
            run_id = subdir.name
            run_dir = Path(subdir.path).resolve()
            instrument_type = instrument.determine_instrument_type(run_id)
            analysis_not_already_initiated = not (config.analysis_output_dir / run_id).exists()
            combine_fastq_complete_json_file_exists = (run_dir / 'combine_fastq_complete.json').exists()
            combine_fastq_complete_output_dir_exists = (run_dir / 'fastq_pass_combined').exists()

            conditions_checked = {
                "is_directory": subdir.is_dir(),
                "supported_run_id_format": instrument_type != InstrumentType.unknown,
                "combine_fastq_complete": combine_fastq_complete_json_file_exists and combine_fastq_complete_output_dir_exists,
                "readable": is_readable(run_dir),
                "analysis_not_already_initiated": analysis_not_already_initiated,
                "not_excluded": run_id not in config.excluded_runs
            }

            if all(conditions_checked.values()):
                log.info({
                    "event_type": "run_directory_found",
                    "sequencing_run_id": run_id,
                    "run_directory_path": str(run_dir),
                })
                log.debug({
                    "event_type": "run_directory_checked",
                    "run_directory_path": str(run_dir),
                    "conditions_checked": conditions_checked
                })
                run = Run(
                    sequencing_run_id=run_id,
                    path=run_dir,
                    instrument_type=instrument_type,
                    fastq_directory=Path(run_dir / 'fastq_pass_combined'),
                )
                yield run
            else:
                log.debug({
                    "event_type": "run_directory_skipped",
                    "run_directory_path": str(run_dir),
                    "conditions_checked": conditions_checked
                })
    

def scan(config: Config) -> Iterator[Run]:
    """
    Scanning involves looking for all existing runs.

    :param config: Application config.
    :return: A run directory to analyze
    """
    log.info({"event_type": "scan_start"})
    yield from find_run_dirs(config)


def analyze_run(config: Config, run: Run):
    """
    Initiate an analysis on one directory of fastq files.
    
    :param config: Application config.
    :type config: dict[str, object]
    :param run: Sequencing run. Keys: ['run_dir', 'sequencing_run_id', 'instrument_type']
    :type run: dict[str, str]
    :return: None
    :rtype: None
    """
    base_analysis_outdir = config.analysis_output_dir
    base_analysis_work_dir = config.analysis_work_dir
    
    for pipeline in config.pipelines:
        log.debug({
            "event_type": "pre_analysis_starting",
            "sequencing_run_id": run.sequencing_run_id,
            "pipeline_name": pipeline.name,
        })

        pipeline = pre_analysis.prepare_analysis(config, pipeline, run)

        log.debug({
            "event_type": "prepare_analysis_complete",
            "sequencing_run_id": run.sequencing_run_id,
            "pipeline_name": pipeline.name
        })

        analysis_dependencies_complete = pre_analysis.check_analysis_dependencies_complete(config, pipeline, run)

        analysis_not_already_started = not os.path.exists(pipeline.parameters['outdir'])
        conditions_checked = {
            'pipeline_dependencies_met': analysis_dependencies_complete,
            'analysis_not_already_started': analysis_not_already_started,
        }

        if not all(conditions_checked.values()):
            log.warning({
                "event_type": "analysis_skipped",
                "pipeline_name": pipeline.name,
                "pipeline_version": pipeline.version,
                "pipeline_dependencies": pipeline.dependencies,
                "sequencing_run_id": run.sequencing_run_id,
                "conditions_checked": conditions_checked,
            })
            continue

        analysis_result = analysis.run_pipeline(config, pipeline, run)

        post_analysis.post_analysis(config, pipeline, run, analysis_result)

        
