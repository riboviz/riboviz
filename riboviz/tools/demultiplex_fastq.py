#!/usr/bin/env python
"""
Demultiplex FASTQ files using UMI-tools-compliant barcodes present
within the FASTQ headers and a sample sheet file.

Usage::

    python -m riboviz.tools.demultiplex_fastq [-h]
        -s SAMPLE_SHEET_FILE -1 READ1_FILE
        [-2 [READ2_FILE]] [-m MISMATCHES] [-o [OUT_DIR]]
        [-d [DELIMITER]]

    -h, --help            show this help message and exit
    -s SAMPLE_SHEET_FILE, --sample-sheet SAMPLE_SHEET_FILE
                          Sample sheet file
    -1 READ1_FILE, --read1 READ1_FILE
                          FASTQ file
    -2 [READ2_FILE], --read2 [READ2_FILE]
                          FASTQ file, for paired reads
    -m MISMATCHES, --mismatches MISMATCHES
                          Number of mismatches permitted in barcode
                          (default 1)
    -o [OUT_DIR], --outdir [OUT_DIR]
                          Output directory
    -d [DELIMITER], --delimiter [DELIMITER]
                          Barcode delimiter (default _)

See :py:func:`riboviz.demultiplex_fastq.demultiplex`.
"""

import argparse
from riboviz import barcodes_umis
from riboviz import demultiplex_fastq
from riboviz import provenance


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Demultiplex FASTQ files using UMI-tools-compliant barcodes present within the FASTQ headers and a sample sheet file")
    parser.add_argument("-s",
                        "--sample-sheet",
                        dest="sample_sheet_file",
                        required=True,
                        help="Sample sheet file")
    parser.add_argument("-1",
                        "--read1",
                        dest="read1_file",
                        required=True,
                        help="FASTQ file")
    parser.add_argument("-2",
                        "--read2",
                        dest="read2_file",
                        default=None,
                        nargs='?',
                        help="FASTQ file, for paired reads")
    parser.add_argument("-m",
                        "--mismatches",
                        dest="mismatches",
                        default=1,
                        type=int,
                        help="Number of mismatches permitted in barcode (default 1)")
    parser.add_argument("-o",
                        "--outdir",
                        dest="out_dir",
                        nargs='?',
                        default="output",
                        help="Output directory")
    parser.add_argument("-d",
                        "--delimiter",
                        dest="delimiter",
                        nargs='?',
                        default=barcodes_umis.BARCODE_DELIMITER,
                        help="Barcode delimiter (default " +
                        barcodes_umis.BARCODE_DELIMITER + ")")
    options = parser.parse_args()
    return options


def invoke_demultiplex_fastq():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.demultiplex_fastq.demultiplex`.
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    sample_sheet_file = options.sample_sheet_file
    read1_file = options.read1_file
    read2_file = options.read2_file
    mismatches = options.mismatches
    out_dir = options.out_dir
    delimiter = options.delimiter
    demultiplex_fastq.demultiplex(sample_sheet_file,
                                  read1_file,
                                  read2_file,
                                  mismatches,
                                  out_dir,
                                  delimiter)


if __name__ == "__main__":
    invoke_demultiplex_fastq()
