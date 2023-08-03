import datetime
import json
import logging
import os
import re
import shutil
import subprocess
import uuid

from typing import Iterator, Optional


def find_run_dirs(config, check_symlinks_complete=True):
    """
    Find sequencing run directories under the 'run_parent_dirs' listed in the config.

    :param config: Application config.
    :type config: dict[str, object]
    :param check_upload_complete: Check for presence of 'upload_complete.json' file.
    :type check_upload_complete: bool
    :return: Run directory. Keys: ['sequencing_run_id', 'run_dir', 'instrument_type']
    :rtype: Iterator[Optional[dict[str, str]]]
    """
    fastq_by_run_dir = config['fastq_by_run_dir']
    gridion_run_id_regex    = "\d{8}_\d{4}_X\d_[A-Z0-9]{8}_[a-z0-9]{8}$"
    promethion_run_id_regex = "\d{8}_\d{4}_P2S_\d+-\d_[A-Z0-9]{8}_[a-z0-9]{8}$"

    subdirs = os.scandir(fastq_by_run_dir)

    for subdir in subdirs:
        run_id = subdir.name
        matches_gridion_regex = re.match(gridion_run_id_regex, run_id)
        matches_promethion_regex = re.match(promethion_run_id_regex, run_id)
        instrument_type = 'unknown'
        if matches_gridion_regex:
            instrument_type = 'gridion'
        elif matches_promethion_regex:
            instrument_type = 'promethion'
        symlinks_complete = os.path.exists(os.path.join(subdir, 'symlinks_complete.json'))
        analysis_not_already_initiated = not os.path.exists(os.path.join(config['analysis_output_dir'], run_id))
        not_excluded = True
        if 'excluded_runs' in config:
            not_excluded = not run_id in config['excluded_runs']

        conditions_checked = {
            "is_directory": subdir.is_dir(),
            "matches_nanopore_run_id_format": ((matches_gridion_regex is not None) or
                                               (matches_promethion_regex is not None)),
            "analysis_not_already_initiated": analysis_not_already_initiated,
            "not_excluded": not_excluded,
        }

        if check_symlinks_complete:
            conditions_checked["symlinks_complete"] = symlinks_complete

        conditions_met = list(conditions_checked.values())
        if all(conditions_met):
            logging.info(json.dumps({"event_type": "run_directory_found", "sequencing_run_id": run_id, "run_directory_path": os.path.abspath(subdir.path)}))
            run = {}
            run['run_dir'] = os.path.abspath(subdir.path)
            run['sequencing_run_id'] = run_id
            run['instrument_type'] = instrument_type
            yield run
        else:
            logging.debug(json.dumps({"event_type": "directory_skipped", "run_directory_path": os.path.abspath(subdir.path), "conditions_checked": conditions_checked}))
            yield None
    

def scan(config: dict[str, object]) -> Iterator[Optional[dict[str, object]]]:
    """
    Scanning involves looking for all existing runs.

    :param config: Application config.
    :type config: dict[str, object]
    :return: A run directory to analyze, or None. Keys: ['run_dir', 'sequencing_run_id', 'instrument_type']
    :rtype: Iterator[Optional[dict[str, object]]]
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    for run_dir in find_run_dirs(config):    
        yield run_dir


def analyze_run(config, run):
    """
    Initiate an analysis on one directory of fastq files.
    
    :param config: Application config.
    :type config: dict[str, object]
    :param run: Sequencing run. Keys: ['run_dir', 'sequencing_run_id', 'instrument_type']
    :type run: dict[str, str]
    :return: None
    :rtype: None
    """
    base_analysis_outdir = config['analysis_output_dir']
    base_analysis_work_dir = config['analysis_work_dir']
    if 'notification_email_addresses' in config:
        notification_email_addresses = config['notification_email_addresses']
    else:
        notification_email_addresses = []
    for pipeline in config['pipelines']:
        pipeline_parameters = pipeline['pipeline_parameters']
        pipeline_short_name = pipeline['pipeline_name'].split('/')[1].replace('_', '-')
        pipeline_minor_version = '.'.join(pipeline['pipeline_version'].split('.')[0:2])
        analysis_timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        analysis_run_id = run['sequencing_run_id']
        analysis_output_dir = os.path.join(base_analysis_outdir, analysis_run_id, pipeline_short_name + '-' + pipeline_minor_version + '-output')
        run['fastq_input'] = run['run_dir']
        run['outdir'] = analysis_output_dir
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
        if 'send_notification_emails' in config and config['send_notification_emails']:
            pipeline_command += ['-with-notification', ','.join(notification_email_addresses)]
        for flag, config_value in pipeline_parameters.items():
            if config_value is None:
                value = run[flag]
            else:
                value = config_value
            pipeline_command += ['--' + flag, value]
            pipeline_command = list(map(str, pipeline_command))
        logging.info(json.dumps({"event_type": "analysis_started", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)}))

        try:
            timestamp_analysis_start = datetime.datetime.now().isoformat()
            subprocess.run(pipeline_command, capture_output=True, check=True)
            timestamp_analysis_complete = datetime.datetime.now().isoformat()
            analysis_complete_path = os.path.join(analysis_output_dir, 'analysis_complete.json')
            analysis_complete = {
                'timestamp_analysis_start': timestamp_analysis_start,
                'timestamp_analysis_complete': timestamp_analysis_complete,
            }
            with open(analysis_complete_path, 'w') as f:
                json.dump(analysis_complete, f, indent=2)
            logging.info(json.dumps({"event_type": "analysis_completed", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)}))
            shutil.rmtree(analysis_work_dir, ignore_errors=True)
            logging.info(json.dumps({"event_type": "analysis_work_dir_deleted", "sequencing_run_id": analysis_run_id, "analysis_work_dir_path": analysis_work_dir}))
        except subprocess.CalledProcessError as e:
            logging.error(json.dumps({"event_type": "analysis_failed", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)}))
        except OSError as e:
            logging.error(json.dumps({"event_type": "delete_analysis_work_dir_failed", "sequencing_run_id": analysis_run_id, "analysis_work_dir_path": analysis_work_dir}))
