#!/usr/bin/env python
"""
Check FASTA and GFF files for compatibility.

Usage::

    python -m riboviz.tools.check_fasta_gff [-h] \
        -f FASTA -g GFF [-o FEATURES_ISSUES] \
        [--use-feature-name] \
        [--feature-format FEATURE_FORMAT]
        [--start-codon START_CODON [START_CODON ...]]

    -h, --help            show this help message and exit
    -f FASTA, --fasta FASTA
                          FASTA file input
    -g GFF, --gff GFF     GFF3 file input
    -o FEATURES_ISSUES, --features-issues FEATURES_ISSUES
                          Issues file output
    --use-feature-name    If a CDS feature defines both 'ID' and 'Name'
                          attributes then use 'Name' in reporting,
                          otherwise use 'ID' (default 'false')
    --feature-format FEATURE_FORMAT
                          Feature name format for features which do
                          not define ``ID`` or ``Name``
                          attributes. This format is applied to the
                          sequence ID to create a feature name.
    --start-codon START_CODON [START_CODON ...]
                          Allowable start codons (default 'ATG')

See :py:func:`riboviz.check_fasta_gff.check_fasta_gff`.
"""
import argparse
from riboviz import check_fasta_gff
from riboviz.fasta_gff import CDS_FEATURE_FORMAT
from riboviz.fasta_gff import START_CODON
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
    parser.add_argument("--use-feature-name",
                        dest="use_feature_name",
                        action='store_true',
                        default=False,
                        help="If a CDS feature defines both 'ID' and 'Name' attributes then use 'Name' in reporting, otherwise use 'ID' (default 'false')")
    parser.add_argument("--feature-format",
                        dest="feature_format",
                        default=CDS_FEATURE_FORMAT,
                        help="Feature name format for features which do not define 'ID'  or 'Name' attributes. This format is applied to the sequence ID to create a feature name.")
    parser.add_argument("--start-codon",
                        dest="start_codon",
                        nargs="+",
                        default=[START_CODON],
                        help="Allowable start codons (default '{}')".format(START_CODON))
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
    use_feature_name = options.use_feature_name
    start_codons = options.start_codon
    check_fasta_gff.check_fasta_gff(fasta,
                                    gff,
                                    features_issues,
                                    feature_format=feature_format,
                                    use_feature_name=use_feature_name,
                                    start_codons=start_codons)


if __name__ == "__main__":
    invoke_check_fasta_gff()
