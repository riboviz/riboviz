#! python

## subsample_bioseqfile.py
## subsamples an input fastq (or other sequencing) file, to produce a smaller file  
## whose reads are randomly sampled from of the input with a fixed probability
## 
## example:
##   python riboviz/tools/subsample_bioseqfile.py -in vignette/input/SRR1042855_s1mi.fastq.gz -prob 0.001 -out /tmp/SRR1042855_s1000.fastq.gz -ftype fastq

import argparse, gzip, os
from Bio import SeqIO
from random import random


def subsample_bioseqfile(seqfilein, seqfileout, filetype, prob, verbose=False):
    """
    Subsample a *gzipped* biological sequence file using Bio.SeqIO
    See https://biopython.org/wiki/SeqIO for description of valid filetypes (fastq, etc)
    
    :param seqfilein: string
    :param seqfileout: string
    :param filetype: string
    :param prob: float
    :param verbose: bool

    :return: Bio.SeqIO
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


if __name__ == "__main__":
    # take input options
    parser = argparse.ArgumentParser(description="Subsample reads from an input fastq file")
    parser.add_argument("-in", dest="filein", help="SeqIO file input")
    parser.add_argument("-out", dest="fileout", help="SeqIO file output")
    parser.add_argument("-ftype", dest="filetype", default="fastq", help="SeqIO filetype (default is fastq)",
                        choices=['fastq'])
    parser.add_argument("-prob", dest="prob", type=float, default=0.01, help="proportion to sample (default 0.01)")
    parser.add_argument("-overwrite", dest="overwrite", action='store_true', help="overwrite output file if it exists")
    parser.add_argument("-verb", dest="verbose", action='store_true', help="print progress statements")
    options = parser.parse_args()

    # files exist, overwrite output?
    if not options.overwrite and os.path.exists(options.fileout):
        raise ValueError("output file {} already exists, use '-overwrite' to overwrite".format(options.fileout))
    if not os.path.exists(options.filein):
        raise ValueError("input file {} doesn't exist".format(options.filein))

    # run the function
    subsample_bioseqfile(seqfilein=options.filein,
                         seqfileout=options.fileout,
                         filetype=options.filetype,
                         prob=options.prob,
                         verbose=options.verbose)