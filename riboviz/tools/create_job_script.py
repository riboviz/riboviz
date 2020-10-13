#!/usr/bin/env python
"""
Create a job submission script using values from a workflow
configuration file and a template job submission script.

Usage::

    python -m riboviz.tools.create_job_script [-h]
        -c CONFIG_FILE -i [INPUT_FILE]
         [-o [OUTPUT_FILE]]

  -h, --help            show this help message and exit
  -c CONFIG_FILE, --configuration CONFIG_FILE
                        Workflow configuration file
  -i [INPUT_FILE], --input_file [INPUT_FILE]
                        Job submission script template
  -o [OUTPUT_FILE], --output [OUTPUT_FILE]
                        Job submission script (if not provided then
                        the job submission script is printed to
                        standard output)

See :py:mod:`riboviz.create_job_script` for details on the expected
format of the template and how it is customised.
"""
import argparse
from riboviz import create_job_script


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
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
    options = parser.parse_args()
    return options


def invoke_create_job_script():
    """
    Parse command-line options then invoke
    :py:mod:`riboviz.create_job_script.create_job_script`.
    """
    options = parse_command_line_options()
    config_file = options.config_file
    input_file = options.input_file
    output_file = options.output_file
    create_job_script.create_job_script(config_file,
                                        input_file,
                                        output_file)


if __name__ == "__main__":
    invoke_create_job_script()
