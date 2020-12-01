#!/usr/bin/env python
"""
Create a job submission script using values from a workflow
configuration file and a template job submission script.

Usage::

    python -m riboviz.tools.create_job_script [-h]
        -i [INPUT_FILE] [-o [OUTPUT_FILE]] [-t [TOKEN_TAG]]
        --r-libs R_LIBS --config-file [CONFIG_FILE]
        [--job-name [JOB_NAME]]
        [--job-runtime [JOB_RUNTIME]]
        [--job-memory [JOB_MEMORY]]
        [--job-num-cpus [JOB_NUM_CPUS]]
        [--job-email [JOB_EMAIL]]
        [--job-email-events [JOB_EMAIL_EVENTS]]
        [--validate-only]
        [--nextflow-work-dir [NEXTFLOW_WORK_DIR]]
        [--nextflow-report-file [NEXTFLOW_REPORT_FILE]]

      -h, --help            show this help message and exit
      -i INPUT_FILE, --input_file INPUT_FILE
                            Job submission script template
      -o [OUTPUT_FILE], --output [OUTPUT_FILE]
                            Job submission script (if not provided
                            then the job submission script is printed
                            to standard output)
      -t [TOKEN_TAG], --token-tag [TOKEN_TAG]
                            Tag marking up tokens for replacement in
                            job submission script template

See :py:mod:`riboviz.create_job_script` for details on the expected
format of the template and how it is customised.
"""
import argparse
from riboviz import create_job_script
from riboviz import params


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: recognised command-line options and extra arguments supplied
    :rtype: argparse.Namespace, list(str or unicode)
    """
    parser = argparse.ArgumentParser(
        description="Create a job submission script using values from a workflow configuration file and a template job submission script")
    parser.add_argument("-i",
                        "--input_file",
                        dest="input_file",
                        required=True,
                        help="Job submission script template")
    parser.add_argument("-o",
                        "--output",
                        dest="output_file",
                        nargs='?',
                        help="Job submission script (if not provided then the job submission script is printed to standard output)")
    parser.add_argument("-t",
                        "--token-tag",
                        dest="token_tag",
                        nargs='?',
                        default=create_job_script.TOKEN_TAG,
                        help="Tag marking up tokens for replacement in job submission script template")
    for param in [params.R_LIBS, params.CONFIG_FILE]:
        parser.add_argument("--" + param.replace("_", "-"),
                            dest=param,
                            required=True)
    return parser.parse_known_args()


def invoke_create_job_script():
    """
    Parse command-line options then invoke
    :py:mod:`riboviz.create_job_script.create_job_script`.

    :raises AssertionError: if any additional parameter provided
    by the user is not in
    :py:const:`riboviz.params.DEFAULT_JOB_CONFIG`.
    """
    options, extras = parse_command_line_options()
    config_file = options.config_file
    input_file = options.input_file
    output_file = options.output_file
    token_tag = options.token_tag
    config = {}
    for param in [params.R_LIBS, params.CONFIG_FILE]:
        config[param] = vars(options)[param]
    # Validate additional parameters and cast to relevant type.
    # Any parameter with value "None" is given value None.
    while len(extras) > 0:
        key = extras.pop(0)
        key = key.strip("-").replace("-", "_")
        assert key in params.DEFAULT_JOB_CONFIG.keys(), \
            "Unknown configuration parameter: {}".format(key)
        value_type = params.JOB_CONFIG_TYPE[key]
        if value_type is bool:
            # Booleans are solo flags without associated values
            config[key] = True
        else:
            value = extras.pop(0)
            if value == "None":
                value = None
            else:
                value = value_type(value)
            config[key] = value
    create_job_script.create_job_script(config_file,
                                        config,
                                        input_file,
                                        output_file,
                                        token_tag)


if __name__ == "__main__":
    invoke_create_job_script()
