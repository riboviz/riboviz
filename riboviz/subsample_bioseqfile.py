"""
Subsample .fastq, .fastq.gz or other sequence file.
"""
import gzip
import os
from random import random
from Bio import SeqIO


def subsample_bioseqfile(
        seqfilein, seqfileout, filetype, prob, overwrite, verbose
):
    """
    Subsample a *gzipped* biological sequence file using Bio.SeqIO
    See https://biopython.org/wiki/SeqIO for description of valid filetypes

    :param seqfilein: File name of input sequence file
    :type seqfilein: str or unicode
    :param seqfileout: File name of input sequence file
    :type seqfileout: str or unicode
    :param filetype: SeqIO file type (default 'fastq')
    :type filetype: str or unicode
    :param prob: probability / proportion to sample (default 0.01)
    :type prob: float
    :param overwrite: overwrite if output file exists? (default False)
    :type overwrite: bool
    :param verbose: print progress statements (default False)
    :type verbose: bool
    :return: Bio.SeqIO
    :rtype:
    :raise FileNotFoundError: If the file cannot be found or is \
    not a file
    """

    ext = os.path.splitext(seqfilein)[1].lower()

    # files exist, overwrite output?
    if os.path.exists(seqfileout) and not overwrite:
        raise ValueError(
            "output file {} already exists, use '-overwrite' to replace"
            .format(seqfileout))
    if not os.path.exists(seqfilein):
        raise ValueError(
            "input file {} doesn't exist".format(seqfilein))

    is_gz = ext in [".gz", ".gzip"]

    if is_gz:
        open_file = gzip.open
        open_r = "rt"
        open_w = "wt"
    else:
        open_file = open
        open_r = "r"
        open_w = "w"

    row_count = 0
    row_count_out = 0

    with open_file(seqfilein, open_r) as in_handle, \
            open_file(seqfileout, open_w) as out_handle:
        for record in SeqIO.parse(in_handle, filetype):
            row_count += 1
            if row_count % 100000 == 0:
                print(("read {rowcount}".format(rowcount=row_count)))
                if random() < prob:
                    row_count_out += 1
                    if verbose:
                        print((record.id))
                    SeqIO.write(record, out_handle, filetype)
    print(("subsampling complete; read {} records from {}, \
        wrote {} records to {}".format(
            row_count, seqfilein, row_count_out, seqfileout
        )))
