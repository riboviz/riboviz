#!/usr/bin/env python
"""
Generate barcode pairs and write each pair plus the Hamming distance
between then to a file of tab-separated values.

Usage:

    python -m riboviz.tools.generate_barcode_pairs \
        <OUTPUT_FILE> <BARCODE_LENGTH>

Example:

    python -m riboviz.tools.generate_barcode_pairs pairs.tsv 3
"""
import csv
import itertools
import sys
from riboviz.tools.demultiplex_fastq import hamming_distance


NUCLEOTIDES = "ACGT"
""" Nucleotide letters """


def generate_barcode_pairs(filename, length=1):
    """
    Generate barcode pairs and write each pair plus the Hamming
    distance between then to a file of tab-separated values.

    :param filename: Filename
    :type filename: str or unicode
    :param length: Barcode length
    :type length: int
    """
    barcodes = [''.join(i) for i in itertools.product(NUCLEOTIDES,
                                                      repeat=length)]
    with open(filename, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for (a, b) in itertools.product(barcodes, repeat=2):
            distance = hamming_distance(a, b)
            writer.writerow([a, b, distance])


if __name__ == "__main__":
    filename = sys.argv[1]
    length = int(sys.argv[2])
    generate_barcode_pairs(filename, length)
