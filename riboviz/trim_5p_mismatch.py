#! python
"""
Remove a single 5' mismatched nt AND filter reads with more than
specified mismatches from a SAM file.
"""
import re
import pysam


NUM_PROCESSED = "num_processed"
""" Key for number of reads processed """
NUM_DISCARDED = "num_discarded"
""" Key for number of reads discarded """
NUM_TRIMMED = "num_trimmed"
""" Key for number of reads trimmed """
NUM_WRITTEN = "num_written"
""" Key for number of reads written """


def increase_soft_clip_init(read):
    """
    Edit CIGAR string of a read to increase soft clip, reducing
    number of initial matches.

    :param read: Read in SAM file.
    :type read: pysam.libcalignedsegment.AlignedSegment
    """
    cigar_string = read.cigarstring
    if "S" not in cigar_string[:3]:
        # Read is not soft-clipped on left.
        # Find number of initial matches.
        num_init_match = int(re.findall(r"^([0-9]+)M", cigar_string)[0])
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


def increase_soft_clip_term(read):
    """
    Edit CIGAR string of a read to increase soft clip, reducing
    number of terminal matches.

    :param read: Read in SAM file.
    :type read: pysam.libcalignedsegment.AlignedSegment
    """
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


def trim_5p_mismatch(sam_file_in,
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
    :return dict with keys "processed", "discarded", "trimmed" and
    "written" and numbers of reads corresponding to each
    :rtype: dict
    """
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
            num_mismatches = read.get_tag('NM')
            if num_mismatches > 0 and fivep_remove:
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
                    num_mismatches -= 1
                    read.set_tag('NM', num_mismatches)
                    increase_soft_clip_init(read)
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
                    num_mismatches -= 1
                    read.set_tag('NM', num_mismatches)
                    increase_soft_clip_term(read)
                    num_trimmed += 1

            if num_mismatches <= max_mismatches:
                num_written += 1
                sam_out.write(read)
            else:
                num_discarded += 1
    print("Number of reads:")
    print(("{}:\t{}".format(NUM_PROCESSED, num_processed)))
    print(("{}:\t{}".format(NUM_DISCARDED, num_discarded)))
    print(("{}:\t{}".format(NUM_TRIMMED, num_trimmed)))
    print(("{}:\t{}".format(NUM_WRITTEN, num_written)))
    return {NUM_PROCESSED: num_processed,
            NUM_DISCARDED: num_discarded,
            NUM_TRIMMED: num_trimmed,
            NUM_WRITTEN: num_written}
