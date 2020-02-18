#!/usr/bin/env python
"""
Create simulated FASTQ files to test UMI/deduplication, adaptor
trimming, and demultiplexing.

Usage::

    python -m riboviz.tools.create_fastq_simdata [-h]
        -o OUTPUT_DIR

    -h, --help            show this help message and exit
    -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                          Output directory

See :py:func:`riboviz.create_fastq_simdata.create_fastq_simdata`.
"""
import argparse
from riboviz import create_fastq_simdata


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Create simulated FASTQ files to test UMI/deduplication, adaptor trimming, and demultiplexing")
    parser.add_argument("-o",
                        "--output-dir",
                        dest="output_dir",
                        required=True,
                        help="Output directory")
    options = parser.parse_args()
    return options


def invoke_create_fastq_simdata():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.create_fastq_simdata.create_fastq_simdata`.
    """
    options = parse_command_line_options()
    output_dir = options.output_dir
    create_fastq_simdata.create_fastq_simdata(output_dir)


if __name__ == "__main__":
    invoke_create_fastq_simdata()
