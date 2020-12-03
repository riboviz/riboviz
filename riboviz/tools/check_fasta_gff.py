#!/usr/bin/env python
"""
Check FASTA and GFF files for compatibility.

Usage::

    python -m riboviz.tools.check_fasta_gff [-h] \
        -f FASTA -g GFF [-o FEATURES_ISSUES]

    -h, --help            show this help message and exit
    -f FASTA, --fasta FASTA
                          fasta file input
    -g GFF, --gff GFF     gff3 file input
    -o FEATURES_ISSUES, --features-issues FEATURES_ISSUES
                          Features issues file output
                          (default features_issues.tsv)

See :py:func:`riboviz.check_fasta_gff.check_fasta_gff`.
"""
import argparse
from riboviz import check_fasta_gff
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
    parser.add_argument("-o",
                        "--features-issues",
                        dest="features_issues",
                        default="features_issues.tsv",
                        help="Features issues file output")
    options = parser.parse_args()
    return options


def invoke_check_fasta_gff():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.check_fasta_gff.check_fasta_gff`.
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    fasta = options.fasta
    gff = options.gff
    features_issues = options.features_issues
    check_fasta_gff.check_fasta_gff(fasta, gff, features_issues)


if __name__ == "__main__":
    invoke_check_fasta_gff()
