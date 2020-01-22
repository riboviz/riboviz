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
import sys
from riboviz.barcodes_umis import generate_barcode_pairs


if __name__ == "__main__":
    filename = sys.argv[1]
    length = int(sys.argv[2])
    generate_barcode_pairs(filename, length)
