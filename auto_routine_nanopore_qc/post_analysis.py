import csv
import datetime
import glob
import json
import logging
import os
import shutil
import subprocess

from . import parsers

from auto_routine_nanopore_qc.model import Config, CustomJSONEncoder, Pipeline, Run

log = logging.getLogger(__name__)

def post_analysis_routine_nanopore_qc(config: Config, pipeline: Pipeline, run: Run):
    """
    Perform post-analysis tasks for the basic-nanopore-qc pipeline.

    :param config: The app config
    :param pipeline: The pipeline
    :param run: The sequencing run
    :return: None
    :rtype: None
    """
    log.info({
        "event_type": "post_analysis_started",
        "sequencing_run_id": run.sequencing_run_id,
        "pipeline_name": pipeline.name,
        "run": run,
    })

    analysis_run_output_dir = os.path.join(config.analysis_output_dir, run.sequencing_run_id)

    pipeline_short_name = pipeline.name.split('/')[1]
    pipeline_minor_version = ''.join(pipeline.version.rsplit('.', 1)[0])
    analysis_pipeline_output_dir = pipeline.parameters.get('outdir')
    log.debug({
        "event_type": "analysis_pipeline_output_dir",
        "sequencing_run_id": run.sequencing_run_id,
        "analysis_pipeline_output_dir": analysis_pipeline_output_dir
    })

    
    return None


def post_analysis_autocycler_nf(config: Config, pipeline: Pipeline, run: Run):
    """
    Perform post-analysis tasks for the BCCDC-PHL/autocycler-nf assembly pipeline.

    :param config: The config
    :param pipeline: The pipeline
    :param run: The sequencing run
    :return: None
    :rtype: None
    """

    log.info({
        "event_type": "post_analysis_started",
        "sequencing_run_id": run.sequencing_run_id,
        "pipeline_name": pipeline.name,
    })

    analysis_run_output_dir = os.path.join(config.analysis_output_dir, run.sequencing_run_id)

    pipeline_short_name = pipeline.name.split('/')[1]
    pipeline_minor_version = ''.join(pipeline.version.rsplit('.', 1)[0])
    analysis_pipeline_output_dir = pipeline.parameters.get('outdir')
    log.debug({
        "event_type": "analysis_pipeline_output_dir",
        "sequencing_run_id": run.sequencing_run_id,
        "analysis_pipeline_output_dir": analysis_pipeline_output_dir
    })

    return None


def post_analysis(config: Config, pipeline: Pipeline, run: Run, analysis_result: subprocess.CompletedProcess):
    """
    Perform post-analysis tasks for a pipeline.

    :param config: The config
    :param pipeline: The pipeline
    :param run: The sequencing run
    :return: None
    """
    pipeline_short_name = pipeline.name.split('/')[1]
    delete_pipeline_work_dir = analysis_result and analysis_result.returncode == 0
    base_analysis_work_dir = config.analysis_work_dir
    work_dir = pipeline.parameters.get('work_dir')

    if work_dir and delete_pipeline_work_dir:
        shutil.rmtree(work_dir, ignore_errors=True)
        log.info({
            "event_type": "analysis_work_dir_deleted",
            "sequencing_run_id": run.sequencing_run_id,
            "analysis_work_dir_path": work_dir
        })
    elif not work_dir or not os.path.exists(work_dir):        
        log.warning({
            "event_type": "analysis_work_dir_not_found",
            "sequencing_run_id": run.sequencing_run_id,
        })

    if pipeline.name == 'BCCDC-PHL/routine-nanopore-qc':
        return post_analysis_routine_nanopore_qc(config, pipeline, run)
    elif pipeline.name == 'BCCDC-PHL/autocycler-nf':
        return post_analysis_autocycler_nf(config, pipeline, run)
    else:
        log.warning({
            "event_type": "post_analysis_not_implemented",
            "sequencing_run_id": run.sequencing_run_id,
            "pipeline_name": pipeline.name
        })
        return None
