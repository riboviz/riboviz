#! python
"""
Remove a single 5' mismatched nt AND filter reads with more than
specified mismatches from a SAM file.

Usage:

    trim_5p_mismatch.py [-h] -i SAMFILEIN -o SAMFILEOUT [-m [MISMATCHES]]
                               [-5 | -k]

Arguments:

*  '-h', --help': show this help message and exit
*  '-i SAMFILEIN', '--input SAMFILEIN': SAM file input
*  '-o SAMFILEOUT', '--output SAMFILEOUT': SAM file output
*  '-m [MISMATCHES]', '--mismatches [MISMATCHES]': Number of
   mismatches to allow (default 1)
*  '-5', '--5p-remove': Remove 5p mismatches (default True)
*  '-k', '--5p-keep': Keep 5p mismatches (default False)

Examples:

    python -m riboviz.tools.trim_5p_mismatch \
        -i testdata_trim_5p_mismatch.sam \
        -o testdata_trim_5p_mismatch_clean.sam
    python -m riboviz.tools.trim_5p_mismatch \
        -i testdata_trim_5pos5neg.sam \
        -o testdata_trim_5pos5neg_clean.sam
    python -m riboviz.tools.trim_5p_mismatch \
        -i data_map1.sam -o data_map1_clean.sam
"""
import argparse
from riboviz import trim_5p_mismatch
from riboviz import provenance


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Remove a single 5' mismatched nt AND filter reads with more than specified mismatches from a SAM file")
    parser.add_argument("-i",
                        "--input",
                        dest="samfilein",
                        required=True,
                        help="SAM file input")
    parser.add_argument("-o",
                        "--output",
                        dest="samfileout",
                        required=True,
                        help="SAM file output")
    parser.add_argument("-m",
                        "--mismatches",
                        dest="mismatches",
                        nargs='?',
                        default=1,
                        type=int,
                        help="Number of mismatches to allow")
    fivep_parser = parser.add_mutually_exclusive_group(required=False)
    fivep_parser.add_argument("-5",
                              "--5p-remove",
                              dest='fivepremove',
                              action='store_true',
                              help="Remove 5p mismatches")
    fivep_parser.add_argument("-k",
                              "--5p-keep",
                              dest='fivepremove',
                              action='store_false',
                              help="Keep 5p mismatches")
    parser.set_defaults(fivepremove=True)
    options = parser.parse_args()
    return options


def invoke_trim_5p_mismatch():
    """
    Parse command-line options then invoke "trim_5p_mismatches".
    """
    print(provenance.get_version(__file__))
    options = parse_command_line_options()
    sam_file_in = options.samfilein
    sam_file_out = options.samfileout
    fivep_remove = options.fivepremove
    max_mismatches = options.mismatches
    print("trim_5p_mismatch.py running")
    trim_5p_mismatch.trim_5p_mismatch(sam_file_in,
                                      sam_file_out,
                                      fivep_remove,
                                      max_mismatches)


if __name__ == "__main__":
    invoke_trim_5p_mismatch()
