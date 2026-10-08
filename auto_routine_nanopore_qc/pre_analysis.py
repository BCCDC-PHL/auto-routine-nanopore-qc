import csv
import datetime
import glob
import json
import logging
import os
import shutil
import subprocess

from auto_routine_nanopore_qc.model import Config, CustomJSONEncoder, Run, Pipeline

log = logging.getLogger(__name__)

def check_analysis_dependencies_complete(config: Config, pipeline: Pipeline, run: Run) -> bool:
    """
    Check that all of the entries in the pipeline's `dependencies` config have completed. If so, return True. Return False otherwise.

    Pipeline completion is determined by the presence of an `analysis_complete.json` file in the analysis output directory.

    :param config: The config dictionary
    :type config: dict
    :param pipeline:
    :type pipeline: dict[str, object]
    :param run: The run dictionary
    :type run: dict
    :return: Whether or not all of the pipelines listed in `dependencies` have completed.
    :rtype: bool
    """
    all_dependencies_complete = False
    dependencies = pipeline.dependencies
    if dependencies is None:
        return True
    dependencies_complete = []
    dependency_infos = []
    base_analysis_output_dir = config.analysis_output_dir
    analysis_run_output_dir = os.path.join(base_analysis_output_dir, run.sequencing_run_id)
    for dependency in dependencies:
        dependency_pipeline_short_name = dependency['pipeline_name'].split('/')[1]
        dependency_pipeline_minor_version = ''.join(dependency['pipeline_version'].rsplit('.', 1)[0])
        dependency_analysis_output_dir_name = '-'.join([dependency_pipeline_short_name, dependency_pipeline_minor_version, 'output'])
        dependency_analysis_complete_path = os.path.join(analysis_run_output_dir, dependency_analysis_output_dir_name, 'analysis_complete.json')
        dependency_analysis_complete = os.path.exists(dependency_analysis_complete_path)
        dependency_info = {
            'pipeline_name': dependency['pipeline_name'],
            'pipeline_version': dependency['pipeline_version'],
            'analysis_complete_path': dependency_analysis_complete_path,
            'analysis_complete': dependency_analysis_complete
        }
        dependency_infos.append(dependency_info)
    dependencies_complete = [dep['analysis_complete'] for dep in dependency_infos]
    log.info({
        "event_type": "checked_analysis_dependencies",
        "sequencing_run_id": run.sequencing_run_id,
        "pipeline_name": pipeline.name,
        "all_analysis_dependencies_complete": all(dependencies_complete),
        "analysis_dependencies": dependency_infos
    })
    if all(dependencies_complete):
        all_dependencies_complete = True

    return all_dependencies_complete


def pre_analysis_routine_nanopore_qc(config: Config, pipeline: Pipeline, run: Run):
    """
    Prepare the BCCDC-PHL/routnine-nanopore-qc analysis pipeline for execution.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :param run: The run dictionary
    :type run: dict
    :return: The prepared pipeline dictionary
    :rtype: dict
    """
    pipeline_short_name = pipeline.name.split('/')[1]
    pipeline_minor_version = ''.join(pipeline.version.rsplit('.', 1)[0])
    
    base_analysis_outdir = config.analysis_output_dir
    pipeline_output_dirname = '-'.join([pipeline_short_name, pipeline_minor_version, 'output'])
    outdir = os.path.abspath(os.path.join(
        base_analysis_outdir,
        run.sequencing_run_id,
        pipeline_output_dirname,
    ))
    pipeline.parameters['fastq_input'] = run.fastq_directory
    pipeline.parameters['outdir'] = outdir

    return pipeline


def pre_analysis_autocycler_nf(config: Config, pipeline: Pipeline, run: Run):
    """
    Prepare the BCCDC-PHL/routine-assembly analysis pipeline for execution.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary
    :type pipeline: dict
    :param run: The run dictionary
    :type run: dict
    :return: The prepared pipeline dictionary
    :rtype: dict
    """
    pipeline_short_name = pipeline.name.split('/')[1]
    pipeline_minor_version = ''.join(pipeline.version.rsplit('.', 1)[0])

    base_analysis_outdir = config.analysis_output_dir
    run_analysis_outdir = os.path.join(base_analysis_outdir, run.sequencing_run_id)
    pipeline_output_dirname = '-'.join([pipeline_short_name, pipeline_minor_version, 'output'])
    pipeline_analysis_output_dir = os.path.join(run_analysis_outdir, pipeline_output_dirname)
    samplesheet_path = os.path.join(
        run_analysis_outdir,
        'samplesheets',
        run.sequencing_run_id + '_autocycler-nf_samplesheet.csv'
    )

    pipeline.parameters['samplesheet_input'] = samplesheet_path
    pipeline.parameters['outdir'] = pipeline_analysis_output_dir

    return pipeline


def prepare_analysis(config: Config, pipeline: Pipeline, run: Run):
    """
    Prepare the pipeline for execution.

    :param config: The config dictionary
    :type config: dict
    :param pipeline: The pipeline dictionary. Expected keys: ['name', 'version', 'parameters']
    :type pipeline: dict
    :param run: The run dictionary. Expected keys: ['sequencing_run_id', 'analysis_parameters']
    :type run: dict
    :return: The prepared pipeline dictionary
    :rtype: dict
    """
    pipeline_short_name = pipeline.name.split('/')[1]
    pipeline_minor_version = ''.join(pipeline.version.rsplit('.', 1)[0])
    pipeline_output_dirname = '-'.join([pipeline_short_name, pipeline_minor_version, 'output'])

    base_analysis_work_dir = config.analysis_work_dir
    analysis_timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
    work_dir = os.path.abspath(os.path.join(base_analysis_work_dir, 'work-' + run.sequencing_run_id + '_' + pipeline_short_name + '_' + analysis_timestamp))
    pipeline.parameters['work_dir'] = work_dir

    base_analysis_outdir = config.analysis_output_dir
    run_analysis_outdir = os.path.join(base_analysis_outdir, run.sequencing_run_id)
    pipeline_output_dir = os.path.abspath(os.path.join(run_analysis_outdir, pipeline_output_dirname))
    pipeline.parameters['outdir'] = pipeline_output_dir

    report_path = os.path.abspath(os.path.join(pipeline_output_dir, run.sequencing_run_id + '_' + pipeline_short_name + '_report.html'))
    pipeline.parameters['report_path'] = report_path

    trace_path = os.path.abspath(os.path.join(pipeline_output_dir, run.sequencing_run_id + '_' + pipeline_short_name + '_trace.tsv'))
    pipeline.parameters['trace_path'] = trace_path

    timeline_path = os.path.abspath(os.path.join(pipeline_output_dir, run.sequencing_run_id + '_' + pipeline_short_name + '_timeline.html'))
    pipeline.parameters['timeline_path'] = timeline_path

    log_path = os.path.abspath(os.path.join(pipeline_output_dir, run.sequencing_run_id + '_' + pipeline_short_name + '_nextflow.log'))
    pipeline.parameters['log_path'] = log_path

    analysis_dependencies_complete = check_analysis_dependencies_complete(config, pipeline, run)

    if not analysis_dependencies_complete:
        log.info({
            "event_type": "analysis_dependencies_incomplete",
            "pipeline_name": pipeline.name,
            "sequencing_run_id": run.sequencing_run_id
        })
        return None

    if pipeline.name == 'BCCDC-PHL/routine-nanopore-qc':
        return pre_analysis_routine_nanopore_qc(config, pipeline, run)
    elif pipeline.name == 'BCCDC-PHL/autocycler-nf':
        return pre_analysis_autocycler_nf(config, pipeline, run)
    else:
        logging.error(json.dumps({
            "event_type": "pipeline_not_supported",
            "pipeline_name": pipeline.name,
            "sequencing_run_id": run.sequencing_run_id
        }))
        return None
