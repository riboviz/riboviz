#!/usr/bin/env python
"""
Compare two files for equality.

Usage:

    python -m riboviz.tools.compare_files [-h] \
        -1 FILE1 -2 FILE2 [-n]

Arguments:

* '-h', '--help': show this help message and exit
* '-1 FILE1', '--file1 FILE1': File1
* '-2 FILE2', '--file2 FILE2': File2
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
    parser.add_argument("-1",
                        "--file1",
                        dest="file1",
                        required=True,
                        help="File1")
    parser.add_argument("-2",
                        "--file2",
                        dest="file2",
                        required=True,
                        help="File2")
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
    file1 = options.file1
    file2 = options.file2
    names = options.names
    compare_files.compare_files(file1, file2, names)


if __name__ == "__main__":
    invoke_compare_files()
