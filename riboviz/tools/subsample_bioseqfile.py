#! python
"""
Subsample an input fastq (or other sequencing) file, to produce a
smaller file whose reads are randomly sampled from of the input with a
fixed probability.

Usage:

   python -m riboviz.tools.subsample_bioseqfile [-h] \
        [-i [FILEIN]] [-o [FILEOUT]] \
        [-t [FILETYPE]] [-p [PROB]] [-v [VERBOSE]]

Example:

    python -m riboviz.tools.subsample_bioseqfile \
        -i vignette/input/SRR1042855_s1mi.fastq \
        -p 0.00001 \
        -o vignette/tmp/SRR1042855_s10.fastq \
        -t fastq
        -v

    python -m riboviz.tools.subsample_bioseqfile \
        -i vignette/input/SRR1042855_s1mi.fastq.gz \
        -p 0.00001 \
        -o vignette/tmp/SRR1042855_s10.fastq.gz \
        -t fastq
"""
import argparse
import gzip
import os
import os.path
from random import random
from Bio import SeqIO
from riboviz import provenance


def subsample_bioseq_file(input_file, prob, output_file, file_type,
                          verbose=False):
    """
    Subsample a biological sequence file using Bio.SeqIO.

    See https://biopython.org/wiki/SeqIO for description of valid
    filetypes (fastq, etc).

    :param input_file: Input file
    :type input_file: str or unicode
    :param prob: Proportion to sample
    :type prob: float
    :param output_file: Output file
    :type output_file: str or unicode
    :param file_type: SeqIO file type
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


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Subsample reads from an input fastq file")
    parser.add_argument("-i",
                        "--input",
                        dest="input_file",
                        required=True,
                        help="SeqIO file input")
    parser.add_argument("-o",
                        "--output",
                        dest="output_file",
                        required=True,
                        help="SeqIO file output")
    parser.add_argument("-t",
                        "--type",
                        dest="file_type",
                        default="fastq",
                        help="SeqIO file type (default 'fastq')")
    parser.add_argument("-p",
                        "--probability",
                        dest="prob",
                        type=float,
                        default=0.01,
                        help="proportion to sample (default 0.01)")
    parser.add_argument("-v",
                        "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="print progress statements")
    options = parser.parse_args()
    return options


def invoke_subsample_bioseq_file():
    """
    Parse command-line options then invoke "subsample_bioseq_file".
    """
    print(provenance.get_provenance_str(__file__))
    options = parse_command_line_options()
    input_file = options.input_file
    output_file = options.output_file
    file_type = options.file_type
    prob = options.prob
    verbose = options.verbose
    subsample_bioseq_file(input_file,
                          prob,
                          output_file,
                          file_type,
                          verbose)


if __name__ == "__main__":
    invoke_subsample_bioseq_file()
