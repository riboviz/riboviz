"""
Barcode and UMI-related utilities.
"""
import csv
import itertools


NUCLEOTIDES = "ACGT"
""" Nucleotide letters """
BARCODE_DELIMITER = "_"
""" Default barcode delmiter in fastq header """
UMI_DELIMITER = "_"
""" Default UMI delmiter in fastq header """


def hamming_distance(str1, str2):
    """
    Returns the hamming distance between two strings.

    :param str1: String
    :type str1: str or unicode
    :param str2: String
    :type str2: str or unicode
    :return: hamming distance
    :rtype: int
    """
    return sum(1 for (a, b) in zip(str1, str2) if a != b)


def barcode_matches(record,
                    barcode,
                    mismatches=0,
                    delimiter=BARCODE_DELIMITER):
    """
    Returns True if fastq record header includes barcode, up to a given
    number of mismatches. The header is assumed to be of form:

        @...<DELIMITER><BARCODE><DELIMITER>...

    If <BARCODE> differs from barcode by greater than mismatches, o
    there is no <BARCODE> or <BARCODE> has a different length than
    barcode then False is returned.

    :param record: Record
    :type record: str or unicode
    :param barcode: Barcode
    :type barcode: str or unicode
    :param mismatches: Number of mismatches
    :type mismatches: int
    :param delimiter: Barcode delimiter
    :type delimiter: str or unicode
    :returns: True or False
    :rtype: bool
    """
    chunks = record.split(delimiter)
    if len(chunks) == 1:
        return False
    candidate = chunks[1]
    if len(candidate) != len(barcode):
        return False
    return hamming_distance(candidate, barcode) <= mismatches


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
