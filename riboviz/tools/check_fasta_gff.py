#!/usr/bin/env python
"""
Check FASTA and GFF files for compatibility.

Usage::

    python -m riboviz.tools.check_fasta_gff [-h] -f FASTA -g GFF

    -h, --help            show this help message and exit
    -f FASTA, --fasta FASTA
                          fasta file input
    -g GFF, --gff GFF     gff3 file input

See :py:func:`riboviz.fasta_gff.check_fasta_gff`.
"""
import argparse
from riboviz import fasta_gff
from riboviz import provenance


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Check FASTA and GFF files for compatibility")
    parser.add_argument("-f",
                        "--fasta",
                        dest="fasta",
                        required=True,
                        help="fasta file input")
    parser.add_argument("-g",
                        "--gff",
                        dest="gff",
                        required=True,
                        help="gff3 file input")
    options = parser.parse_args()
    return options


def invoke_check_fasta_gff():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.fasta_gff.check_fasta_gff`.
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    fasta = options.fasta
    gff = options.gff
    fasta_gff.check_fasta_gff(fasta, gff)


if __name__ == "__main__":
    invoke_check_fasta_gff()
