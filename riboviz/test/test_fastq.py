"""
riboviz.fastq tests.
"""
import gzip
import itertools
import os
import tempfile
import pytest
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from riboviz import fastq
from riboviz.barcodes_umis import NUCLEOTIDES


@pytest.fixture(scope="function")
def fastq_file():
    """
    Create a temporary FASTQ file to write data to.

    :return: path to temporary FASTQ file
    :rtype: str or unicode
    """
    _, fastq_file = tempfile.mkstemp(prefix="tmp", suffix=".fastq")
    yield fastq_file
    if os.path.exists(fastq_file):
        os.remove(fastq_file)


@pytest.fixture(scope="function")
def fastq_gz_file():
    """
    Create a temporary FASTQ.GZ file to write data to.

    :return: path to temporary FASTQ.GZ file
    :rtype: str or unicode
    """
    _, fastq_file = tempfile.mkstemp(prefix="tmp", suffix=".fastq.gz")
    yield fastq_file
    if os.path.exists(fastq_file):
        os.remove(fastq_file)


@pytest.mark.parametrize("extension", [".gz", ".gzip", ".GZ", ".GZIP"])
def test_is_fastq_gz(extension):
    """
    Test is_fastq_gz with GZIP extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert fastq.is_fastq_gz("sample.fastq{}".format(extension))


@pytest.mark.parametrize("extension", [".txt", ""])
def test_not_is_fastq_gz(extension):
    """
    Test is_fastq_gz with non-GZIP extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert not fastq.is_fastq_gz("sample.fastq.{}".format(extension))


@pytest.mark.parametrize("extension", [".gz", ".gzip", ".GZ", ".GZIP"])
def test_strip_fastq_gz(extension):
    """
    Test strip_fastq_gz with GZIP extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert fastq.strip_fastq_gz("sample.fastq{}".format(extension)) == \
        "sample.fastq"


@pytest.mark.parametrize("extension", [".txt", ""])
def test_not_strip_fastq_gz(extension):
    """
    Test that strip_fastq_gz with non-GZIP extensions returns the
    original file_name.

    :param extension: Extension
    :type extension: str or unicode
    """
    file_name = "sample.fastq{}".format(extension)
    assert file_name == fastq.strip_fastq_gz(file_name)


def test_get_fastq_filename_gz_false():
    """
    Test get_fastq_filename with is_gz=False.
    """
    assert fastq.get_fastq_filename("sample") == "sample.fastq"


def test_get_fastq_filename_gz_true():
    """
    Test get_fastq_filename with is_gz=True.
    """
    assert fastq.get_fastq_filename("sample", True) == "sample.fastq.gz"


def test_get_fastq_filenames_gz_false():
    """
    Test get_fastq_filenames with is_gz=False.
    """
    file_names = ["sample{}".format(i) for i in range(0, 3)]
    expected_names = ["sample{}.fastq".format(i) for i in range(0, 3)]
    actual_names = fastq.get_fastq_filenames(file_names)
    assert expected_names == actual_names


def test_get_fastq_filenames_gz_true():
    """
    Test get_fastq_filenames with is_gz=True.
    """
    file_names = ["sample{}".format(i) for i in range(0, 3)]
    expected_names = ["sample{}.fastq.gz".format(i) for i in range(0, 3)]
    actual_names = fastq.get_fastq_filenames(file_names, True)
    assert expected_names == actual_names


def get_test_fastq_records(read_length, count):
    """
    Get FASTQ records for test FASTQ files.

    :param read_length: Read length
    :type read_length: int
    :param count: Number of records
    :type count: int
    :return: List of records, consisting of the first count reads of
    read_length found by enumerating combinations of
    riboviz.barcodes_umis.NUCLEOTIDES, with quality scores each
    [0,...,read_kength - 1]
    :rtype: list(Bio.SeqRecord.SeqRecord)
    """
    # Create a list of "count" reads AAAA, AAAC etc.
    reads = [''.join(i)
             for i in itertools.product(NUCLEOTIDES,
                                        repeat=read_length)][0:count]
    # Create a list of SeqRecords
    records = [SeqRecord(Seq(read, IUPAC.ambiguous_dna),
                         id="read{}".format(i),
                         name="read{}".format(i),
                         description="read{}".format(i))
               for read, i in zip(reads, range(0, len(reads)))]
    quality = list(range(0, read_length))
    for record in records:
        record.letter_annotations["phred_quality"] = quality
    return records


@pytest.mark.parametrize("count", [0, 1, 10])
def test_count_records(fastq_file, count):
    """
    Test count_records.

    :param fastq_file: path to temporary FASTQ file
    :type fastq_file: str or unicode
    :param count: Number of records
    :type count: int
    """
    records = get_test_fastq_records(4, count)
    with open(fastq_file, "wt") as f:
        SeqIO.write(records, f, "fastq")
    assert fastq.count_records(fastq_file) == count


@pytest.mark.parametrize("count", [0, 1, 10])
def test_count_records_gz(fastq_gz_file, count):
    """
    Test count_records with a .FASTQ.GZ file.

    :param fastq_gz_file: path to temporary FASTQ.GZ file
    :type fastq_gz_file: str or unicode
    :param count: Number of records
    :type count: int
    """
    records = get_test_fastq_records(4, count)
    with gzip.open(fastq_gz_file, "wt") as f:
        SeqIO.write(records, f, "fastq")
    assert fastq.count_records(fastq_gz_file) == count
