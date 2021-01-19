"""
Subsample .fastq, .fastq.gz or other sequence file.
"""
import argparse, gzip, os
from Bio import SeqIO
from random import random


def subsample_bioseqfile(seqfilein, seqfileout, filetype, prob, overwrite, verbose):
    """
    Subsample a *gzipped* biological sequence file using Bio.SeqIO
    See https://biopython.org/wiki/SeqIO for description of valid filetypes (fastq, etc)

    :param seqfilein: File name of input sequence file
    :type seqfilein: str or unicode
    :param seqfileout: File name of input sequence file
    :type seqfileout: str or unicode
    :param filetype: SeqIO file type (default 'fastq')
    :type filetype: str or unicode
    :param prob: probability / proportion to sample (default 0.01)
    :type prob: float
    :param overwrite: overwrite output file if it already exists? (default store_true)
    :type overwrite: bool
    :param verbose: print progress statements (default False)
    :type verbose: bool
    :return: Bio.SeqIO
    :rtype:
    :raise FileNotFoundError: If the file cannot be found or is \
    not a file
    """
    filein_ext = os.path.splitext(seqfilein)[1]

    if filein_ext == ".gz" or filein_ext == ".gzip":
        in_handle = gzip.open(seqfilein, "rt")
        out_handle = gzip.open(seqfileout, "wt")
    else:
        in_handle = open(seqfilein, "r")
        out_handle = open(seqfileout, "w")

    row_count = 0
    row_count_out = 0
    for record in SeqIO.parse(in_handle, filetype):
        row_count += 1
        if row_count % 100000 == 0:
            print("read {rowcount}".format(rowcount=row_count))
        if random() < prob:
            row_count_out += 1
            if verbose:
                print(record.id)
            SeqIO.write(record, out_handle, filetype)
    print("subsampling complete; read {} records from {}, wrote {} records to {}".format(
        row_count, seqfilein, row_count_out, seqfileout
    ))
