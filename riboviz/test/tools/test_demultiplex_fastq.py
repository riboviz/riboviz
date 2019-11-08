"""
demultiplex_fastq.py test suite.
"""
from riboviz.tools.demultiplex_fastq import barcode_matches
from riboviz.tools.demultiplex_fastq import hamming_distance


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


def test_barcode_matches_false():
    """
    Test barcode_matches with a non-matching record.
    """
    record = "@X1:Tag_AAA_ 1:N:0:XXXXXXXX"
    barcode = "AAC"
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
