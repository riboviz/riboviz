"""
SAM and BAM-related constants and functions.
"""
import os
import pysam


PG_TAG = "PG"
""" SAM file PG (program) tag """
BAM_FORMAT = "{}.bam"
""" Format string for BAM files """
BAM_BAI_FORMAT = "{}.bai"
"""
Format string for BAM index file based on default behaviour of
"samtools index"
"""


def is_bam(file_name):
    """
    Does the given file end with .bam or .BAM?

    :param file_name: File name
    :type file_name: str or unicode
    :return: True if file_name ends with .bam or .BAM, False otherwise
    :rtype: bool
    """
    _, ext = os.path.splitext(os.path.basename(file_name))
    return ext.lower() == ".bam"


def is_sam(file_name):
    """
    Does the given file end with .sam or .SAM?

    :param file_name: File name
    :type file_name: str or unicode
    :return: True if file_name ends with .sam or .SAM, False otherwise
    :rtype: bool
    """
    _, ext = os.path.splitext(os.path.basename(file_name))
    return ext.lower() == ".sam"


def count_sequences(file_name):
    """
    Count number of sequences and mapped (primary aligned) sequences
    in a SAM or BAM file.

    :param file_name: SAM/BAM file name
    :type file_name: str or unicode
    :return: (number of sequences, number of mapped sequences)
    :rtype: tuple(int, int)
    """
    if is_bam(file_name):
        mode = "rb"
    else:
        mode = "r"
    num_sequences = 0
    num_mapped_sequences = 0
    with pysam.AlignmentFile(file_name, mode=mode) as f:
        for sequence in f:
            num_sequences = num_sequences + 1
            if sequence.is_unmapped or sequence.is_secondary:
                continue
            num_mapped_sequences = num_mapped_sequences + 1
    return (num_sequences, num_mapped_sequences)
