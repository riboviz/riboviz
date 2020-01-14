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
import re
import pysam
from riboviz import provenance


def trim_5p_mismatches(sam_file_in,
                       sam_file_out,
                       fivep_remove=True,
                       max_mismatches=1):
    """
    Remove a single 5' mismatched nt AND filter reads with more than
    specified mismatches from a SAM file.

    :param sam_file_in: SAM file input
    :type sam_file_in: str or unicode
    :param sam_file_out: SAM file output
    :type sam_file_out: str or unicode
    :param fivep_remove: Remove mismatched 5' nt?
    :type fivep_remove: bool
    :param max_mismatches: Number of mismatches to allow
    :type max_mismatches: int
    """
    # TODO with...open
    num_processed = 0
    num_discarded = 0
    num_trimmed = 0
    num_written = 0
    with pysam.AlignmentFile(sam_file_in, "r") as sam_in,\
        pysam.AlignmentFile(sam_file_out, "wh", template=sam_in) as sam_out:
        for read in sam_in.fetch():
            num_processed += 1
            if (num_processed % 1000000) == 1:
                print(("processed " + str(num_processed - 1) + " reads"))
            try:
                # Get MD tag for read, encoding mismatches.
                md_tag = read.get_tag('MD')
            except KeyError:
                # MD tag not present, assume read not aligned, discard.
                num_discarded += 1
                continue
            # Count mismatches in read.
            num_mismatch = read.get_tag('NM')

            if num_mismatch > 0 and fivep_remove:
                # If there are any mismatches...

                if md_tag[0] == "0" and read.flag == 0:
                    # If the 5' nt is mismatched on a plus-strand read...
                    if md_tag[2] in ["A", "T", "C", "G", "0"]:
                        # 2nd nt is also mismatched; discard.
                        num_discarded += 1
                        continue
                    # ... soft-clip 5' nt
                    # Increment position of alignment.
                    read.pos += 1
                    # Edit MD tag to remove leading mismatch.
                    read.set_tag('MD', re.sub("^0[ATCG]", "", md_tag))
                    # Lower the number of mismatches by 1.
                    num_mismatch -= 1
                    read.set_tag('NM', num_mismatch)
                    # Edit CIGAR string to increase soft clip.
                    cigar_string = read.cigarstring
                    if "S" not in cigar_string[:3]:
                        # Read is not soft-clipped on left.
                        # Find number of initial matches.
                        num_init_match = int(re.findall(r"^([0-9]+)M",
                                                        cigar_string)[0])
                        # Add initial "1S" and reduce initial match by 1.
                        new_init_bit = "1S" + str(num_init_match - 1) + "M"
                        read.cigarstring = re.sub(r"^([0-9]+)M",
                                                  new_init_bit,
                                                  cigar_string)
                    else:
                        # Read is soft-clipped on left.
                        num_soft_clip = int(re.findall(r"^([0-9]+)S",
                                                       cigar_string)[0])
                        num_init_match = int(re.findall(r"([0-9]+)M",
                                                        cigar_string)[0])
                        # Add 1 to left soft-clip and reduce initial match by 1.    
                        new_init_bit = str(num_soft_clip + 1) + \
                            "S" + str(num_init_match - 1) + "M"
                        read.cigarstring = re.sub(r"^([0-9]+)S([0-9]+)M",
                                                  new_init_bit,
                                                  cigar_string)
                    # Count the read as trimmed.
                    num_trimmed += 1

                if bool(re.search("[ATCG]0$", md_tag)) and read.flag == 16:
                    # If the 5' nt is mismatched on a minus strand read...
                    # Positive sense is with template.
                    # Read is reverse-complement.
                    if md_tag[-3] in ["A", "T", "C", "G"]:
                        # 2nd nt is also mismatched; discard.
                        num_discarded += 1
                        continue
                    # ... soft-clip 5' nt
                    # Don't increment position of alignment!
                    # Edit MD tag to remove trailing mismatch.
                    read.set_tag('MD', re.sub("[ATCG]0$", "", md_tag))
                    # Lower the number of mismatches by 1.
                    num_mismatch -= 1
                    read.set_tag('NM', num_mismatch)
                    # Edit CIGAR string to increase soft clip.
                    cigar_string = read.cigarstring
                    if cigar_string[-1] != "S":
                        # Read is not soft-clipped on left (= right along template)
                        # Find number of terminal matches.
                        num_term_match = int(re.findall(r"([0-9]+)M$",
                                                        cigar_string)[0])
                        # Add initial "1S" and reduce initial match by 1.
                        new_term_bit = str(num_term_match - 1) + "M1S"
                        read.cigarstring = re.sub(r"([0-9]+)M$",
                                                  new_term_bit,
                                                  cigar_string)
                    else:
                        # Read is soft-clipped on left (= right along template).
                        # Find number of terminal matches.
                        num_soft_clip = int(re.findall(r"([0-9]+)S$",
                                                       cigar_string)[0])
                        num_term_match = int(re.findall(r"([0-9]+)M",
                                                        cigar_string)[-1])
                        # Add initial "1S" and reduce initial match by 1.
                        new_term_bit = str(num_term_match - 1) + \
                                       "M" + str(num_soft_clip + 1) + "S"
                        read.cigarstring = re.sub(r"([0-9]+)M([0-9]+)S$",
                                                  new_term_bit,
                                                  cigar_string)
                    # Count the read as trimmed.
                    num_trimmed += 1

            if num_mismatch <= max_mismatches:
                num_written += 1
                sam_out.write(read)
            else:
                num_discarded += 1

    print("trim_5p_mismatch.py finished, number of reads")
    print(("processed:\t" + str(num_processed)))
    print(("discarded:\t" + str(num_discarded)))
    print(("trimmed:\t" + str(num_trimmed)))
    print(("written:\t" + str(num_written)))


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


def main():
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
    trim_5p_mismatches(sam_file_in,
                       sam_file_out,
                       fivep_remove,
                       max_mismatches)


if __name__ == "__main__":
    main()
