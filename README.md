# auto-routine-nanopore-qc
Automated analysis of nanopore sequence data using the [BCCDC-PHL/routine-nanopore-qc](https://github.com/BCCDC-PHL/routine-nanopore-qc) pipeline.

# Installation

# Usage
Start the tool as follows:

```bash
auto-routine-nanopore-qc --config config.json
```

See the Configuration section of this document for details on preparing a configuration file.

More detailed logs can be produced by controlling the log level using the `--log-level` flag:

```bash
auto-routine-nanopore-qc --config config.json --log-level debug
```

# Configuration
This tool takes a single config file, in JSON format, with the following structure:

```json
{
  "fastq_by_run_dir": "/path/to/fastq_symlinks_by_run",
  "analysis_output_dir": "/path/to/analysis_by_run",
  "analysis_work_dir": "/path/to/auto-routine-nanopore-qc-work",
  "notification_email_addresses": [
	"someone@example.org",
	"someone_else@example.org"
  ],
  "send_notification_emails": true,
  "scan_interval_seconds": 3600,
  "pipelines": [
    {
      "pipeline_name": "BCCDC-PHL/routine-nanopore-qc",
      "pipeline_version": "v0.1.0",
      "pipeline_parameters": {
	    "fastq_input": null,
		"kraken2_db": "/path/to/kraken2_db",
		"outdir": null
      }
    }
  ]
}
```

# Logging
This tool outputs [structured logs](https://www.honeycomb.io/blog/structured-logging-and-your-team/) in [JSON Lines](https://jsonlines.org/) format:

Every log line should include the fields:

- `timestamp`
- `level`
- `module`
- `function_name`
- `line_num`
- `message`

...and the contents of the `message` key will be a JSON object that includes at `event_type`. The remaining keys inside the `message` will vary by event type.

```json
{"timestamp": "2022-09-22T11:32:52.287", "level": "INFO", "module", "core", "function_name": "scan", "line_num", 56, "message": {"event_type": "scan_start"}}
```
