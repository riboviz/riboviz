#!/usr/bin/env python
"""
Create a job submission script using values from a workflow
configuration file and a template job submission script.

Usage::

    python -m riboviz.tools.create_job_script [-h]
        -c CONFIG_FILE -i [INPUT_FILE]
         [-o [OUTPUT_FILE]] [-t [TOKEN_TAG]]

      -h, --help            show this help message and exit
      -c CONFIG_FILE, --configuration CONFIG_FILE
                            Workflow configuration file
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
import os.path
import sys
from riboviz import create_job_script


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options and extra configuration values
    provided at the command-line
    :rtype: (argparse.Namespace, dict(str or unicode => str or unicode))
    """
    parser = argparse.ArgumentParser(
        description="Create a job submission script using values from a workflow configuration file and a template job submission script")
    parser.add_argument("-c",
                        "--configuration",
                        dest="config_file",
                        required=True,
                        help="Workflow configuration file")
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
    options, extras = parser.parse_known_args()
    if len(extras) % 2:
        parser.print_usage(sys.stderr)
        error_msg = "{}: error: expected every configuration parameters to have a value\n".format(os.path.basename(__file__))
        sys.stderr.write(error_msg)
        sys.exit(2)
    extra_config = dict(zip(extras[::2], extras[1::2]))
    extra_config = {k.strip("-"): v for k, v in extra_config.items()}
    return options, extra_config


def invoke_create_job_script():
    """
    Parse command-line options then invoke
    :py:mod:`riboviz.create_job_script.create_job_script`.
    """
    options, extra_config = parse_command_line_options()
    config_file = options.config_file
    input_file = options.input_file
    output_file = options.output_file
    token_tag = options.token_tag
    create_job_script.create_job_script(config_file,
                                        extra_config,
                                        input_file,
                                        output_file,
                                        token_tag)


if __name__ == "__main__":
    try:
        invoke_create_job_script()
    except Exception as e:
        print(e)
