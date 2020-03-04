#!/usr/bin/env python
"""
Subsample an input FASTQ (or other sequencing) file, to produce a
smaller file whose reads are randomly sampled from of the input with a
fixed probability.

Usage::

    python -m riboviz.tools.subsample_bioseqfile [-h]
        -i INPUT_FILE -o OUTPUT_FILE
        [-t FILE_TYPE] [-p PROB] [-v]

    -h, --help            show this help message and exit
    -i INPUT_FILE, --input INPUT_FILE
                          SeqIO file input
    -o OUTPUT_FILE, --output OUTPUT_FILE
                          SeqIO file output
    -t FILE_TYPE, --type FILE_TYPE
                          SeqIO file type (default 'fastq')
    -p PROB, --probability PROB
                          proportion to sample (default 0.01)
    -v, --verbose         print progress statements

Examples::

    python -m riboviz.tools.subsample_bioseqfile
        -i vignette/input/SRR1042855_s1mi.fastq
        -p 0.00001
        -o vignette/tmp/SRR1042855_s10.fastq
        -t fastq
        -v

    python -m riboviz.tools.subsample_bioseqfile
        -i vignette/input/SRR1042855_s1mi.fastq
        -p 0.00001
        -o vignette/tmp/SRR1042855_s10.fastq.gz
        -t fastq

See :py:func:`riboviz.subsample_bioseqfile.subsample_bioseqfile`.
"""
import argparse
from riboviz import provenance
from riboviz import subsample_bioseqfile


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Subsample an input FASTQ (or other sequencing) file, to produce a smaller file whose reads are randomly sampled from of the input with a fixed probability")
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


def invoke_subsample_bioseqfile():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.subsample_bioseqfile.subsample_bioseqfile`.
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    input_file = options.input_file
    output_file = options.output_file
    file_type = options.file_type
    prob = options.prob
    verbose = options.verbose
    subsample_bioseqfile.subsample_bioseqfile(input_file,
                                              prob,
                                              output_file,
                                              file_type,
                                              verbose)


if __name__ == "__main__":
    invoke_subsample_bioseqfile()
