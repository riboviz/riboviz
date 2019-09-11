#! python

## subsample_bioseqfile.py
## subsamples an input fastq (or other sequencing) file, to produce a smaller file  
## whose reads are randomly sampled from of the input with a fixed probability
## 
## example:
##   python pyscripts/subsample_bioseqfile.py -in vignette/input/SRR1042855_s1mi.fastq.gz -prob 0.001 -out vignette/tmp/SRR1042855_s1000.fastq.gz -ftype fastq

import argparse, gzip, os
from Bio import SeqIO
from random import random 
from os.path import splitext


def subsample_bioseqfile(seqfilein, prob, seqfileout, filetype, verbose=False):
    """
    Subsample a biological sequence file using Bio.SeqIO
    See https://biopython.org/wiki/SeqIO for description of valid filetypes (fastq, etc)

    :param seqfilein: string
    :param prob: float
    :param seqfileout: string
    :param filetype: string
    :return: Bio.SeqIO
    """
    with open(seqfilein, "r") as in_handle, open(seqfileout,"w") as out_handle:
        for record in SeqIO.parse(in_handle, filetype) :
            if random() < prob :
                if verbose : 
                    print(record.id)
                SeqIO.write(record, out_handle, filetype)
    if verbose : 
        print("subsampling complete")

def subsample_bioseqfile_gz(seqfilein, prob, seqfileout, filetype, verbose=False):
    """
    Subsample a *gzipped* biological sequence file using Bio.SeqIO
    See https://biopython.org/wiki/SeqIO for description of valid filetypes (fastq, etc)
    
    :param seqfilein: string
    :param prob: float
    :param seqfileout: string
    :return: Bio.SeqIO
    """
    with gzip.open(seqfilein, "r") as in_handle, gzip.open(seqfileout,"w") as out_handle:
        for record in SeqIO.parse(in_handle, filetype) :
            if random() < prob :
                if verbose : 
                    print(record.id)
                SeqIO.write(record, out_handle, filetype)
    if verbose : 
        print("subsampling complete")

# These are the lines we used to test the functions
""""
subsample_bioseqfile(seqfilein = "vignette/input/SRR1042855_s1mi.fastq", 
prob = 0.00001, 
seqfileout= "vignette/tmp/SRR1042855_s10.fastq",
filetype="fastq",
verbose=True)

subsample_bioseqfile_gz(seqfilein = "vignette/input/SRR1042855_s1mi.fastq.gz", 
prob = 0.00001, 
seqfileout= "vignette/tmp/SRR1042855_s10.fastq.gz",
filetype="fastq",
verbose=True)
""""

if __name__=="__main__" :
    # take input options
    parser = argparse.ArgumentParser(description="Subsample reads from an input fastq file")
    parser.add_argument("-in", dest="filein", nargs='?', help="SeqIO file input")
    parser.add_argument("-out", dest="fileout", nargs='?', help="SeqIO file output")
    parser.add_argument("-ftype", dest="filetype", nargs='?', default="fastq", help="SeqIO filetype (default is fastq)")
    parser.add_argument("-prob", dest="prob", type=float, nargs='?', default=0.01, help="proportion to sample (default 0.01)")
    parser.add_argument("-verb", dest="verbose", type=bool, nargs='?', default=False, help="print progress statements")
    options = parser.parse_args()
    
    # Tests we should have: files exist, all arguments present, overwrite output?
    
    filein = options.filein
    
    filein_ext = os.path.splitext(filein)[1] 
    
    if filein_ext == ".gz" or filein_ext == ".gzip" :
        subsample_bioseqfile_gz(seqfilein = options.filein, 
                                prob = options.prob, 
                                seqfileout= options.fileout,
                                filetype=options.filetype,
                                verbose=options.verbose)
    else : 
        subsample_bioseqfile(seqfilein = options.filein, 
                             prob = options.prob, 
                             seqfileout= options.fileout,
                             filetype=options.filetype,
                             verbose=options.verbose)
