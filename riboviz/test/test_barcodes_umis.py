"""
riboviz.barcodes_umis test suite.
"""
import csv
import os
import tempfile
import pytest
from riboviz.barcodes_umis import NUCLEOTIDES
from riboviz.barcodes_umis import barcode_matches
from riboviz.barcodes_umis import hamming_distance
from riboviz.barcodes_umis import generate_barcode_pairs


@pytest.fixture(scope="function")
def temporary_file():
    """
    Create a temporary file with a ".dat" suffix.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp", suffix=".dat")
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


def test_hamming_distance_empty():
    """
    Test hamming_distance with empty strings.
    """
    assert hamming_distance("", "") == 0


def test_hamming_distance_equal_characters():
    """
    Test hamming_distance with equal characters.
    """
    assert hamming_distance("A", "A") == 0


def test_hamming_distance_nonequal_characters():
    """
    Test hamming_distance with non-equal characters.
    """
    assert hamming_distance("A", "T") == 1


def test_hamming_distance_equal_strings():
    """
    Test hamming_distance with equal strings.
    """
    assert hamming_distance("GATTACCA", "GATTACCA") == 0


def test_hamming_distance_one():
    """
    Test hamming_distance with strings 1 apart.
    """
    assert hamming_distance("GATTACCA", "GATTGCCA") == 1


def test_hamming_distance_eight():
    """
    Test hamming_distance with strings 8 apart.
    """
    assert hamming_distance("GATTACCA", "CTAATGGT") == 8


def test_barcode_matches():
    """
    Test barcode_matches with default mismatches and delimiter.
    """
    record = "@X1:Tag_AAA_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert barcode_matches(record, barcode)


def test_barcode_matches_no_match():
    """
    Test barcode_matches with a non-matching record.
    """
    record = "@X1:Tag_AAA_ 1:N:0:XXXXXXXX"
    barcode = "AAC"
    assert not barcode_matches(record, barcode)


def test_barcode_matches_different_length_barcode():
    """
    Test barcode_matches with a record with a different length barcode.
    """
    record = "@X1:Tag_AAAA_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert not barcode_matches(record, barcode)


def test_barcode_matches_delimiter():
    """
    Test barcode_matches with a non-default delimiter.
    """
    record = "@X1:Tag.AAA. 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert barcode_matches(record, barcode, delimiter=".")


def test_barcode_matches_no_barcode():
    """
    Test barcode_matches with a record with no barcode.
    """
    record = "@X1:Tag 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert not barcode_matches(record, barcode)


def test_barcode_matches_one_mismatch():
    """
    Test barcode_matches with 1 allowed mismatch.
    """
    record = "@X1:Tag_AAC_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert barcode_matches(record, barcode, 1)


def test_barcode_matches_one_mismatch_false():
    """
    Test barcode_matches with 1 allowed mismatch and a non-matching
    record.
    """
    record = "@X1:Tag_ACC_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert not barcode_matches(record, barcode, 1)


def test_barcode_matches_two_mismatch():
    """
    Test barcode_matches with 1 allowed mismatch and a non-matching
    record.
    """
    record = "@X1:Tag_ACC_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert barcode_matches(record, barcode, 2)


def test_generate_barcode_pairs_0(temporary_file):
    """
    Test generating barcode pairs of length 0 creates an empty file.

    :param temporary_file: Temporary file
    :type temporary_file: str or unicode
    """
    generate_barcode_pairs(temporary_file, length=0)
    with open(temporary_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter="\t")
        rows = [row for row in csv_reader]
    assert len(rows) == 0, "Expected zero rows"


@pytest.mark.parametrize("delimiter", ["\t", ","])
def test_generate_barcode_pairs_1(temporary_file, delimiter):
    """
    Test generating barcode pairs of length 1.

    :param temporary_file: Temporary file
    :type temporary_file: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    generate_barcode_pairs(temporary_file, length=1, delimiter=delimiter)
    with open(temporary_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=delimiter)
        rows = [row for row in csv_reader]
    assert len(rows) == 2 ** len(NUCLEOTIDES)
    for row in rows:
        for nucleotide in row[0:2]:
            assert nucleotide in NUCLEOTIDES, \
                ("{} is not in {}".format(nucleotide, str(NUCLEOTIDES)))
        if row[0] == row[1]:
            assert int(row[2]) == 0,\
                "Hamming distance of {} and {} is not 0".format(row[0], row[1])
        else:
            assert int(row[2]) == 1,\
                "Hamming distance of {} and {} is not 1".format(row[0], row[1])
