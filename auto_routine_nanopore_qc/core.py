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

            conditions_checked = {
                "is_directory": subdir.is_dir(),
                "supported_run_id_format": instrument_type != InstrumentType.unknown,
                "combine_fastq_complete": (run_dir / 'combine_fastq_complete.json').exists(),
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
                run = Run(sequencing_run_id=run_id, path=run_dir, instrument_type=instrument_type)
                print(json.dumps(asdict(run), indent=2, cls=CustomJSONEncoder))
                exit()
                yield run
            else:
                log.debug({
                    "event_type": "directory_skipped",
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
        pipeline_parameters = pipeline['pipeline_parameters']
        pipeline_short_name = pipeline['pipeline_name'].split('/')[1].replace('_', '-')
        pipeline_minor_version = '.'.join(pipeline['pipeline_version'].split('.')[0:2])
        analysis_timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        analysis_run_id = run.sequencing_run_id
        analysis_output_dir = os.path.join(config.analysis_output_dir, analysis_run_id, pipeline_short_name + '-' + pipeline_minor_version + '-output')
        pipeline_parameters['fastq_input'] = os.path.join(run.path, 'fastq_pass_combined')
        pipeline_parameters['outdir'] = analysis_output_dir
        analysis_work_dir = os.path.abspath(os.path.join(base_analysis_work_dir, 'work-' + analysis_run_id + '-' + analysis_timestamp))
        analysis_trace_path = os.path.abspath(os.path.join(base_analysis_outdir, analysis_run_id, pipeline_short_name + '-' + pipeline_minor_version + '-output', analysis_run_id + '_trace.tsv'))
        analysis_report_path = os.path.abspath(os.path.join(base_analysis_outdir, analysis_run_id, pipeline_short_name + '-' + pipeline_minor_version + '-output', analysis_run_id + '_nextflow_report.html'))
        pipeline_command = [
            'nextflow',
            'run',
            pipeline['pipeline_name'],
            '-r', pipeline['pipeline_version'],
            '-profile', 'conda',
            '--cache', os.path.join(os.path.expanduser('~'), '.conda/envs'),
            '-work-dir', analysis_work_dir,
            '-with-trace', analysis_trace_path,
        ]
        if 'send_notification_emails' in config.notification and config.notification['send_notification_emails']:
            if 'recipient_email_addresses' in config.notification:
                pipeline_command += ['-with-notification', ','.join(config.notification['recipient_email_addresses'])]

        for flag, config_value in pipeline_parameters.items():
            if config_value is None:
                value = run[flag]
            else:
                value = config_value
            pipeline_command += ['--' + flag, value]
            pipeline_command = list(map(str, pipeline_command))
        log.info({"event_type": "analysis_started", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)})

        try:
            os.makedirs(analysis_work_dir, exist_ok=True)
            timestamp_analysis_start = datetime.datetime.now().isoformat()
            subprocess.run(pipeline_command, capture_output=True, check=True, cwd=analysis_work_dir)
            timestamp_analysis_complete = datetime.datetime.now().isoformat()
            analysis_complete_path = os.path.join(analysis_output_dir, 'analysis_complete.json')
            analysis_complete = {
                'timestamp_analysis_start': timestamp_analysis_start,
                'timestamp_analysis_complete': timestamp_analysis_complete,
            }
            with open(analysis_complete_path, 'w') as f:
                json.dump(analysis_complete, f, indent=2)
            log.info({"event_type": "analysis_completed", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)})
            shutil.rmtree(analysis_work_dir, ignore_errors=True)
            log.info({"event_type": "analysis_work_dir_deleted", "sequencing_run_id": analysis_run_id, "analysis_work_dir_path": analysis_work_dir})
        except subprocess.CalledProcessError as e:
            log.error({"event_type": "analysis_failed", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)})
        except OSError as e:
            log.error({"event_type": "delete_analysis_work_dir_failed", "sequencing_run_id": analysis_run_id, "analysis_work_dir_path": analysis_work_dir})
