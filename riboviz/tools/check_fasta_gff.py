#!/usr/bin/env python
"""
Check FASTA and GFF files for compatibility.

Usage::

    python -m riboviz.tools.check_fasta_gff [-h] \
        -f FASTA -g GFF [-o FEATURES_ISSUES] \
        [--feature-format FEATURE_FORMAT]

    -h, --help            show this help message and exit
    -f FASTA, --fasta FASTA
                          FASTA file input
    -g GFF, --gff GFF     GFF3 file input
    -o FEATURES_ISSUES, --features-issues FEATURES_ISSUES
                          Issues file output
    --feature-format FEATURE_FORMAT
                          Feature name format for features which do
                          not define ``ID`` or ``Name``
                          attributes. This format is applied to the
                          sequence ID to create a feature name.

See :py:func:`riboviz.check_fasta_gff.check_fasta_gff`.
"""
import argparse
from riboviz import check_fasta_gff
from riboviz.fasta_gff import CDS_FEATURE_FORMAT
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
                        help="FASTA file input")
    parser.add_argument("-g",
                        "--gff",
                        dest="gff",
                        required=True,
                        help="GFF3 file input")
    parser.add_argument("-o",
                        "--features-issues",
                        dest="features_issues",
                        default="features_issues.tsv",
                        help="Issues file output")
    parser.add_argument("--feature-format",
                        dest="feature_format",
                        default=CDS_FEATURE_FORMAT,
                        help="Feature name format for features which do not define 'ID'  or 'Name' attributes. This format is applied to the sequence ID to create a feature name.")
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
    feature_format = options.feature_format
    check_fasta_gff.check_fasta_gff(fasta, gff, features_issues,
                                    feature_format)


if __name__ == "__main__":
    invoke_check_fasta_gff()
