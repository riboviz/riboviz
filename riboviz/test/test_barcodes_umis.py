"""
:py:mod:`riboviz.barcodes_umis` tests.
"""
import csv
import os
import tempfile
import pytest
from riboviz import barcodes_umis


@pytest.fixture(scope="function")
def tmp_file():
    """
    Create a temporary file with a ``dat`` suffix.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp", suffix=".dat")
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


def test_hamming_distance_empty():
    """
    Test :py:func:`riboviz.barcodes_umis.hamming_distance` with empty
    strings.
    """
    assert barcodes_umis.hamming_distance("", "") == 0


def test_hamming_distance_equal_characters():
    """
    Test :py:func:`riboviz.barcodes_umis.hamming_distance` with equal
    characters.
    """
    assert barcodes_umis.hamming_distance("A", "A") == 0


def test_hamming_distance_nonequal_characters():
    """
    Test :py:func:`riboviz.barcodes_umis.hamming_distance` with
    non-equal characters.
    """
    assert barcodes_umis.hamming_distance("A", "T") == 1


def test_hamming_distance_equal_strings():
    """
    Test :py:func:`riboviz.barcodes_umis.hamming_distance` with equal
    strings.
    """
    assert barcodes_umis.hamming_distance("GATTACCA", "GATTACCA") == 0


def test_hamming_distance_one():
    """
    Test :py:func:`riboviz.barcodes_umis.hamming_distance` with
    strings distance 1 apart.
    """
    assert barcodes_umis.hamming_distance("GATTACCA", "GATTGCCA") == 1


def test_hamming_distance_eight():
    """
    Test :py:func:`riboviz.barcodes_umis.hamming_distance` with
    strings distance 8 apart.
    """
    assert barcodes_umis.hamming_distance("GATTACCA", "CTAATGGT") == 8


def test_barcode_matches():
    """
    Test :py:func:`riboviz.barcodes_umis.barcode_matches` with default
    mismatches and delimiter.
    """
    record = "@X1:Tag_AAA_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert barcodes_umis.barcode_matches(record, barcode)


def test_barcode_matches_no_match():
    """
    Test :py:func:`riboviz.barcodes_umis.barcode_matches` with a
    non-matching record.
    """
    record = "@X1:Tag_AAA_ 1:N:0:XXXXXXXX"
    barcode = "AAC"
    assert not barcodes_umis.barcode_matches(record, barcode)


def test_barcode_matches_different_length_barcode():
    """
    Test :py:func:`riboviz.barcodes_umis.barcode_matches` with a
    record with that has a barcode a different length from that being
    matched.
    """
    record = "@X1:Tag_AAAA_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert not barcodes_umis.barcode_matches(record, barcode)


def test_barcode_matches_delimiter():
    """
    Test :py:func:`riboviz.barcodes_umis.barcode_matches` with a
    non-default delimiter.
    """
    record = "@X1:Tag.AAA. 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert barcodes_umis.barcode_matches(record, barcode, delimiter=".")


def test_barcode_matches_no_barcode():
    """
    Test :py:func:`riboviz.barcodes_umis.barcode_matches` with a
    record with no barcode.
    """
    record = "@X1:Tag 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert not barcodes_umis.barcode_matches(record, barcode)


def test_barcode_matches_one_mismatch():
    """
    Test :py:func:`riboviz.barcodes_umis.barcode_matches` with 1
    allowed mismatch.
    """
    record = "@X1:Tag_AAC_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert barcodes_umis.barcode_matches(record, barcode, 1)


def test_barcode_matches_one_mismatch_false():
    """
    Test :py:func:`riboviz.barcodes_umis.barcode_matches` with 1
    allowed mismatch and a non-matching record.
    """
    record = "@X1:Tag_ACC_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert not barcodes_umis.barcode_matches(record, barcode, 1)


def test_barcode_matches_two_mismatch():
    """
    Test :py:func:`riboviz.barcodes_umis.barcode_matches` with 2
    allowed mismatches and a non-matching record.
    """
    record = "@X1:Tag_ACC_ 1:N:0:XXXXXXXX"
    barcode = "AAA"
    assert barcodes_umis.barcode_matches(record, barcode, 2)


def test_create_barcode_pairs_0(tmp_file):
    """
    Test :py:func:`riboviz.barcodes_umis.create_barcode_pairs` of
    length 0 creates an empty file.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    barcodes_umis.create_barcode_pairs(tmp_file, length=0)
    with open(tmp_file, 'r', newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter="\t")
        rows = [row for row in csv_reader]
    assert len(rows) == 0, "Expected zero rows"


@pytest.mark.parametrize("delimiter", ["\t", ","])
def test_create_barcode_pairs_1(tmp_file, delimiter):
    """
    Test :py:func:`riboviz.barcodes_umis.create_barcode_pairs` of
    length 1.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    barcodes_umis.create_barcode_pairs(tmp_file,
                                       length=1,
                                       delimiter=delimiter)
    with open(tmp_file, 'r', newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=delimiter)
        rows = [row for row in csv_reader]
    assert len(rows) == 2 ** len(barcodes_umis.NUCLEOTIDES)
    for row in rows:
        for nucleotide in row[0:2]:
            assert nucleotide in barcodes_umis.NUCLEOTIDES, \
                ("{} is not in {}".format(nucleotide,
                                          str(barcodes_umis.NUCLEOTIDES)))
        if row[0] == row[1]:
            assert int(row[2]) == 0,\
                "Hamming distance of {} and {} is not 0".format(row[0],
                                                                row[1])
        else:
            assert int(row[2]) == 1,\
                "Hamming distance of {} and {} is not 1".format(row[0],
                                                                row[1])
