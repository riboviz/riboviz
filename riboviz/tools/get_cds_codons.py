#!/usr/bin/env python
"""
Extract coding sequence codons and export as a tab-separated values
file.

Usage::

    python -m riboviz.tools.get_cds_codons [-h] \
        -f FASTA -g GFF [-c CDS_CODONS] [-e] \
        [--use-feature-name] \
        [--cds-feature-format CDS_FEATURE_FORMAT]

    -h, --help            show this help message and exit
    -f FASTA, --fasta FASTA
                          FASTA file input
    -g GFF, --gff GFF     GFF3 file input
    -c CDS_CODONS, --cds-codons CDS_CODONS
                          Coding sequence codons file output
    -e, --exclude-stop-codons
                          Exclude stop codons (default false)
    --use-feature-name    If a CDS feature defines both 'ID' and 'Name'
                          attributes then use 'Name' in reporting,
                          otherwise use 'ID' (default 'false')
    --cds-feature-format CDS_FEATURE_FORMAT
                          CDS feature name format for CDS features
                          which do not define ``ID`` or
                          ``Name`` attributes. This format is applied
                          to the sequence ID to create a feature
                          name.

See :py:func:`riboviz.get_cds_codons.get_cds_codons_file` for
information on the tab-separated values file format.
"""
import argparse
from pyfaidx import FastaIndexingError
from riboviz import get_cds_codons
from riboviz.fasta_gff import CDS_FEATURE_FORMAT
from riboviz import provenance


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Extract coding sequence codons and export as a tab-separated values file.")
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
    parser.add_argument("-c",
                        "--cds-codons",
                        dest="cds_codons",
                        default="cds_codons.tsv",
                        help="Coding sequence codons file output")
    parser.add_argument("-e",
                        "--exclude-stop-codons",
                        dest="exclude_stop_codons",
                        action='store_true',
                        help="Exclude stop codons (default false)")
    parser.add_argument("--use-feature-name",
                        dest="use_feature_name",
                        action='store_true',
                        default=False,
                        help="If a CDS feature defines both 'ID' and 'Name' attributes then use 'Name` in reporting, otherwise use 'ID' (default false)")
    parser.add_argument("--cds-feature-format",
                        dest="cds_feature_format",
                        default=CDS_FEATURE_FORMAT,
                        help="CDS feature name format for CDS features which do not define 'ID'  or 'Name' attributes. This format is applied to the sequence ID to create a feature name.")
    options = parser.parse_args()
    return options


def invoke_get_cds_codons():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.get_cds_codons.get_cds_codons_file`.
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    fasta = options.fasta
    gff = options.gff
    cds_codons = options.cds_codons
    exclude_stop_codons = options.exclude_stop_codons
    cds_feature_format = options.cds_feature_format
    use_feature_name = options.use_feature_name
    try:
        get_cds_codons.get_cds_codons_file(fasta,
                                           gff,
                                           cds_codons,
                                           exclude_stop_codons,
                                           cds_feature_format,
                                           use_feature_name)
    except FastaIndexingError as e:
        print("{}: {}".format(type(e).__name__, e))
    except FileNotFoundError as e:
        print("{}: {}".format(type(e).__name__, e))
    except ValueError as e:
        print("{}: {}".format(type(e).__name__, e))


if __name__ == "__main__":
    invoke_get_cds_codons()
