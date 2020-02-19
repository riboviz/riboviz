#!/usr/bin/env python
"""
Upgrade workflow configuration file to be compatible with current
configuration.

Usage::

    python -m riboviz.tools.upgrade_config_file [-h]
        -i INPUT_FILE [-o [OUTPUT_FILE]]

    -h, --help            show this help message and exit
    -i INPUT_FILE, --input INPUT_FILE
                          Input file
    -o [OUTPUT_FILE], --output [OUTPUT_FILE]
                          Output file (if not provided then the
                          upgraded configuration is printed)

See :py:mod:`riboviz.upgrade_config.upgrade_config_file`.
"""
import argparse
from riboviz import upgrade_config


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Upgrade workflow configuration file to be compatible with current configuration")
    parser.add_argument("-i",
                        "--input",
                        dest="input_file",
                        required=True,
                        help="Input file")
    parser.add_argument("-o",
                        "--output",
                        dest="output_file",
                        nargs='?',
                        help="Output file (if not provided then the upgraded configuration is printed)")
    options = parser.parse_args()
    return options


def invoke_upgrade_config_file():
    """
    Parse command-line options then invoke
    :py:mod:`riboviz.upgrade_config.upgrade_config_file`.
    """
    options = parse_command_line_options()
    input_file = options.input_file
    output_file = options.output_file
    upgrade_config.upgrade_config_file(input_file, output_file)


if __name__ == "__main__":
    invoke_upgrade_config_file()
