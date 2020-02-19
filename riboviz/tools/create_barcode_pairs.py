#!/usr/bin/env python
"""
Create barcode pairs and write each pair plus the Hamming distance
between then to a file of tab-separated values.

Usage::

    python -m riboviz.tools.create_barcode_pairs [-h]
        -o OUTPUT_FILE -l LENGTH

    -h, --help            show this help message and exit
    -o OUTPUT_FILE, --output-file OUTPUT_FILE
                          Output file
    -l LENGTH, --length LENGTH
                          Barcode length (default 3)

See :py:func:`riboviz.barcodes_umis.create_barcode_pairs`.
"""
import argparse
from riboviz import barcodes_umis


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Create barcode pairs and write each pair plus the Hamming distance between then to a file of tab-separated values")
    parser.add_argument("-o",
                        "--output-file",
                        dest="output_file",
                        required=True,
                        help="Output file")
    parser.add_argument("-l",
                        "--length",
                        dest="length",
                        default=3,
                        type=int,
                        help="Barcode length (default 3)")
    options = parser.parse_args()
    return options


def invoke_create_barcode_pairs():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.barcodes_umis.create_barcode_pairs`.
    """
    options = parse_command_line_options()
    output_file = options.output_file
    length = options.length
    barcodes_umis.create_barcode_pairs(output_file, length)


if __name__ == "__main__":
    invoke_create_barcode_pairs()
