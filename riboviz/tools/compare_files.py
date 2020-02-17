#!/usr/bin/env python
"""
Compare two files for equality.

Usage:

    python -m riboviz.tools.compare_files [-h] \
        -i INPUT_FILE -o OUTPUT_FILE [-n]

Arguments:

* '-h', '--help': show this help message and exit
* '-i INPUT_FILE', '--input-file OINPUT_FILE': Input file
* '-o OUTPUT_FILE', '--output-file OUTPUT_FILE': Output file
* '-n', '--names': Compare file names
"""
import argparse
from riboviz import compare_files


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Compare two files for equality")
    parser.add_argument("-i",
                        "--input-file",
                        dest="input_file",
                        required=True,
                        help="Input file")
    parser.add_argument("-o",
                        "--output-file",
                        dest="output_file",
                        required=True,
                        help="Output file")
    parser.add_argument("-n",
                        "--names",
                        dest='names',
                        action='store_true',
                        help="Compare file names")
    options = parser.parse_args()
    return options


def invoke_compare_files():
    """
    Parse command-line options then invoke "compare_files".
    """
    options = parse_command_line_options()
    input_file = options.input_file
    output_file = options.output_file
    names = options.names
    compare_files.compare_files(input_file, output_file, names)


if __name__ == "__main__":
    invoke_compare_files()
