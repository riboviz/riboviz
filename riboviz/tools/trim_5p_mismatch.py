#! python
"""
Remove a single 5' mismatched nt AND filter reads with more than
specified mismatches from a SAM file.

Usage:

    python -m riboviz.tools.trim_5p_mismatch [-h] \
        -i SAM_FILE_IN -o SAM_FILE_OUT \
        [-m [MAX_MISMATCHES]] [-5 | -k] [-s SUMMARY_FILE]

Arguments:

* '-h', --help': show this help message and exit
* '-i SAM_FILE_IN', '--input SAM_FILE_IN': SAM file input
* '-o SAM_FILE_OUT', '--output SAM_FILE_OUT': SAM file output
* '-m [MAX_MISMATCHES]', '--max-ismatches [MAX_MISMATCHES]': Number of
  mismatches to allow (default 1)
* '-5', '--5p-remove': Remove 5p mismatches (default True)
* '-k', '--5p-keep': Keep 5p mismatches (default False)
* '-s SUMMARY_FILE', '--summary-file SUMMARY_FILE': Summary file for
  TSV summary of reads processed, discarded, trimmed and written
  (default "trim_5p_mismatch.tsv")
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
                        dest="sam_file_in",
                        required=True,
                        help="SAM file input")
    parser.add_argument("-o",
                        "--output",
                        dest="sam_file_out",
                        required=True,
                        help="SAM file output")
    parser.add_argument("-m",
                        "--max-mismatches",
                        dest="max_mismatches",
                        nargs='?',
                        default=1,
                        type=int,
                        help="Number of mismatches to allow")
    fivep_parser = parser.add_mutually_exclusive_group(required=False)
    fivep_parser.add_argument("-5",
                              "--5p-remove",
                              dest='fivep_remove',
                              action='store_true',
                              help="Remove 5p mismatches")
    fivep_parser.add_argument("-k",
                              "--5p-keep",
                              dest='fivep_remove',
                              action='store_false',
                              help="Keep 5p mismatches")
    parser.set_defaults(fivepremove=True)
    parser.add_argument("-s",
                        "--summary-file",
                        dest="summary_file",
                        default=trim_5p_mismatch.TRIM_5P_MISMATCH_FILE,
                        help="Summary file output")
    options = parser.parse_args()
    return options


def invoke_trim_5p_mismatch():
    """
    Parse command-line options then invoke "trim_5p_mismatches".
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    sam_file_in = options.sam_file_in
    sam_file_out = options.sam_file_out
    fivep_remove = options.fivep_remove
    max_mismatches = options.max_mismatches
    summary_file = options.summary_file
    trim_5p_mismatch.trim_5p_mismatch_file(sam_file_in,
                                           sam_file_out,
                                           fivep_remove,
                                           max_mismatches,
                                           summary_file)


if __name__ == "__main__":
    invoke_trim_5p_mismatch()
