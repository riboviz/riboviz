"""
Subsample an input FASTQ (or other sequencing) file, to produce a
smaller file whose reads are randomly sampled from of the input with a
fixed probability.
"""
import gzip
import os
import os.path
from random import random
from Bio import SeqIO


def subsample_bioseqfile(input_file, prob, output_file, file_type,
                         verbose=False):
    """
    Subsample an input FASTQ (or other sequencing) file, to produce a
    smaller file whose reads are randomly sampled from of the input
    with a fixed probability.

    See https://biopython.org/wiki/SeqIO for description of valid
    filetypes (``fastq``, etc).

    :param input_file: Input file
    :type input_file: str or unicode
    :param prob: Proportion to sample
    :type prob: float
    :param output_file: Output file
    :type output_file: str or unicode
    :param file_type: `Bio.SeqIO` file type
    :type file_type: str or unicode
    :param verbose: Print progress statements?
    :type verbose: bool
    """
    ext = os.path.splitext(input_file)[1]
    is_gz = ext in [".gz", ".gzip"]
    if is_gz:
        open_file = gzip.open
        open_r = "rt"
        open_w = "wt"
    else:
        open_file = open
        open_r = "r"
        open_w = "w"
    with open_file(input_file, open_r) as in_handle, \
        open_file(output_file, open_w) as out_handle:
        for record in SeqIO.parse(in_handle, file_type):
            if random() < prob:
                if verbose:
                    print(record.id)
                SeqIO.write(record, out_handle, file_type)
    if verbose:
        print("subsampling complete")
