"""
riboviz.fastq tests.
"""
import gzip
import itertools
import os
import tempfile
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from riboviz import fastq
from riboviz.barcodes_umis import NUCLEOTIDES


@pytest.fixture(scope="function", params=fastq.FASTQ_EXTS)
def fastq_file(request):
    """
    Create a temporary FASTQ file to write data to.

    :param request: pytest SubRequest with param member which has \
    FASTQ file extension
    :type request: _pytest.fixtures.SubRequest
    :return: path to temporary FASTQ file
    :rtype: str or unicode
    """
    _, fastq_file = tempfile.mkstemp(prefix="tmp",
                                     suffix="." + request.param)
    yield fastq_file
    if os.path.exists(fastq_file):
        os.remove(fastq_file)


@pytest.fixture(scope="function", params=fastq.FASTQ_GZ_EXTS)
def fastq_gz_file(request):
    """
    Create a temporary FASTQ.GZ file to write data to.

    :param request: pytest SubRequest with param member which has \
    FASTQ file extension
    :type request: _pytest.fixtures.SubRequest
    :return: path to temporary FASTQ.GZ file
    :rtype: str or unicode
    """
    _, fastq_file = tempfile.mkstemp(prefix="tmp",
                                     suffix="." + request.param)
    yield fastq_file
    if os.path.exists(fastq_file):
        os.remove(fastq_file)


@pytest.mark.parametrize("extension", [fastq.FASTQ_GZ_FORMAT,
                                       fastq.FQ_GZ_FORMAT,
                                       fastq.FASTQ_GZ_FORMAT.upper(),
                                       fastq.FQ_GZ_FORMAT.upper()])
def test_is_fastq_gz(extension):
    """
    Test is_fastq_gz with GZIP extensions.

    :param extension: extension file format
    :type extension: str or unicode
    """
    assert fastq.is_fastq_gz(extension.format("sample"))


@pytest.mark.parametrize("extension", [".txt", ""])
def test_not_is_fastq_gz(extension):
    """
    Test is_fastq_gz with non-GZIP extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert not fastq.is_fastq_gz("sample{}".format(extension))


@pytest.mark.parametrize("extension",
                         [(fastq.FASTQ_GZ_FORMAT,
                           fastq.FASTQ_FORMAT),
                          (fastq.FQ_GZ_FORMAT,
                           fastq.FQ_FORMAT),
                          (fastq.FASTQ_GZ_FORMAT.upper(),
                           fastq.FASTQ_FORMAT.upper()),
                          (fastq.FQ_GZ_FORMAT.upper(),
                           fastq.FQ_FORMAT.upper())])
def test_strip_fastq_gz(extension):
    """
    Test strip_fastq_gz with GZIP extensions.

    :param extension: (GZIP extension file format, \
    corresponding non-GZIP extension file format)
    :type extension: tuple(str or unicode, str or unicode)
    """
    gz_ext, ext = extension
    assert fastq.strip_fastq_gz(gz_ext.format("sample")) == \
        ext.format("sample")


@pytest.mark.parametrize("extension", [".txt", ""])
def test_not_strip_fastq_gz(extension):
    """
    Test that strip_fastq_gz with non-GZIP extensions returns the
    original file_name.

    :param extension: Extension
    :type extension: str or unicode
    """
    file_name = "sample{}".format(extension)
    assert file_name == fastq.strip_fastq_gz(file_name)


def get_test_fastq_sequences(read_length, count):
    """
    Get FASTQ sequences for test FASTQ files.

    :param read_length: Read length
    :type read_length: int
    :param count: Number of sequences
    :type count: int
    :return: List of sequences, consisting of the first count reads \
    of read_length found by enumerating combinations of \
    riboviz.barcodes_umis.NUCLEOTIDES, with quality scores each \
    [0,...,read_kength - 1]
    :rtype: list(Bio.SeqRecord.SeqRecord)
    """
    # Create a list of "count" reads AAAA, AAAC etc.
    reads = [''.join(i)
             for i in itertools.product(NUCLEOTIDES,
                                        repeat=read_length)][0:count]
    # Create a list of SeqRecords
    sequences = [SeqRecord(Seq(read),
                           id="read{}".format(i),
                           name="read{}".format(i),
                           description="read{}".format(i))
                 for read, i in zip(reads, range(0, len(reads)))]
    quality = list(range(0, read_length))
    for sequence in sequences:
        sequence.letter_annotations["phred_quality"] = quality
    return sequences


@pytest.mark.parametrize("count", [0, 1, 10])
def test_count_sequences(fastq_file, count):
    """
    Test count_sequences.

    :param fastq_file: path to temporary FASTQ file
    :type fastq_file: str or unicode
    :param count: Number of sequences
    :type count: int
    """
    sequences = get_test_fastq_sequences(4, count)
    with open(fastq_file, "wt") as f:
        SeqIO.write(sequences, f, "fastq")
    assert fastq.count_sequences(fastq_file) == count


@pytest.mark.parametrize("count", [0, 1, 10])
def test_count_sequences_gz(fastq_gz_file, count):
    """
    Test count_sequences with a .FASTQ.GZ file.

    :param fastq_gz_file: path to temporary FASTQ.GZ file
    :type fastq_gz_file: str or unicode
    :param count: Number of sequences
    :type count: int
    """
    sequences = get_test_fastq_sequences(4, count)
    with gzip.open(fastq_gz_file, "wt") as f:
        SeqIO.write(sequences, f, "fastq")
    assert fastq.count_sequences(fastq_gz_file) == count
