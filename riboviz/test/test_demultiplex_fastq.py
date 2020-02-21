"""
:py:mod:`riboviz.demultiplex_fastq` tests.
"""
from io import StringIO
from contextlib import ExitStack
import gzip
import os
import shutil
import tempfile
import pytest
from riboviz import demultiplex_fastq
from riboviz import fastq
from riboviz import utils
import riboviz.test


FASTQ_RECORD1 = ["@X1:Tag_AAC_ 1:N:0:XXXXXXXX\n",
                 "GATTACCA\n",
                 "+\n",
                 "IIIIIIII\n"]
""" Sample FASTQ record. """
FASTQ_RECORD2 = ["@X1:Tag_AAC_ 1:N:0:XXXXXXXX\n",
                 "AAAAAAAA\n",
                 "+\n",
                 "IIIIIIII\n"]
""" Sample FASTQ record """


@pytest.fixture(scope="function")
def tmp_dir():
    """
    Create a temporary directory.

    :return: path to temporary directory
    :rtype: str or unicode
    """
    tmp_dir = tempfile.mkdtemp(__name__)
    yield tmp_dir
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def test_assign_sample():
    """
    Test :py:func:`riboviz.demultiplex_fastq.assign_sample`
    with a record with a matching barcode.
    """
    with StringIO() as read1_fh:
        barcode = "AAA"
        is_assigned = demultiplex_fastq.assign_sample(
            FASTQ_RECORD1, None,
            barcode,
            read1_fh, None,
            False, 1, "_")
        assert is_assigned
        assert "".join(FASTQ_RECORD1) == read1_fh.getvalue()


def test_assign_sample_no_match():
    """
    Test :py:func:`riboviz.demultiplex_fastq.assign_sample`
    with a record with a non-matching barcode.
    """
    with StringIO() as read1_fh:
        barcode = "GGG"
        is_assigned = demultiplex_fastq.assign_sample(
            FASTQ_RECORD1, None,
            barcode,
            read1_fh, None,
            False, 1, "_")
        assert not is_assigned
        assert read1_fh.getvalue() == ""


def test_assign_sample_paired_end():
    """
    Test :py:func:`riboviz.demultiplex_fastq.assign_sample`
    with a record and a paired end record with a matching barcode.
    """
    with StringIO() as read1_fh, StringIO() as read2_fh:
        barcode = "AAA"
        is_assigned = demultiplex_fastq.assign_sample(
            FASTQ_RECORD1, FASTQ_RECORD2,
            barcode,
            read1_fh, read2_fh,
            True, 1, "_")
        assert is_assigned
        assert "".join(FASTQ_RECORD1) == read1_fh.getvalue()
        assert "".join(FASTQ_RECORD2) == read2_fh.getvalue()


def test_assign_sample_paired_end_no_match():
    """
    Test :py:func:`riboviz.demultiplex_fastq.assign_sample`
    with a record and a paired end record with a non-matching
    barcode.
    """
    with StringIO() as read1_fh, StringIO() as read2_fh:
        barcode = "GGG"
        is_assigned = demultiplex_fastq.assign_sample(
            FASTQ_RECORD1, FASTQ_RECORD2,
            barcode,
            read1_fh, read2_fh,
            True, 1, "_")
        assert not is_assigned
        assert read1_fh.getvalue() == ""
        assert read2_fh.getvalue() == ""


def test_assign_samples():
    """
    Test :py:func:`riboviz.demultiplex_fastq.assign_samples` with
    paired ends records and matching barcodes.
    """
    with ExitStack() as stack:
        read1_fhs = [stack.enter_context(StringIO()) for f in range(2)]
        read2_fhs = [stack.enter_context(StringIO()) for f in range(2)]
        barcodes = ["CCC", "AAA"]
        num_reads = [0] * len(barcodes)
        is_assigned = demultiplex_fastq.assign_samples(
            FASTQ_RECORD1, FASTQ_RECORD2,
            barcodes,
            read1_fhs, read2_fhs,
            True,
            num_reads,
            1, "_")
        assert is_assigned
        assert num_reads[0] == 0
        assert num_reads[1] == 1
        assert read1_fhs[0].getvalue() == ""
        assert read2_fhs[0].getvalue() == ""
        assert "".join(FASTQ_RECORD1) == read1_fhs[1].getvalue()
        assert "".join(FASTQ_RECORD2) == read2_fhs[1].getvalue()


def test_assign_samples_no_match():
    """
    Test :py:func:`riboviz.demultiplex_fastq.assign_samples` with
    paired end records and non-matching barcodes.
    """
    with ExitStack() as stack:
        read1_fhs = [stack.enter_context(StringIO()) for f in range(2)]
        read2_fhs = [stack.enter_context(StringIO()) for f in range(2)]
        barcodes = ["GGG", "TTT"]
        num_reads = [0] * len(barcodes)
        is_assigned = demultiplex_fastq.assign_samples(
            FASTQ_RECORD1, FASTQ_RECORD2,
            barcodes,
            read1_fhs, read2_fhs,
            True,
            num_reads,
            1, "_")
        assert not is_assigned
        assert num_reads[0] == 0
        assert num_reads[1] == 0
        assert read1_fhs[0].getvalue() == ""
        assert read2_fhs[0].getvalue() == ""
        assert read1_fhs[1].getvalue() == ""
        assert read2_fhs[1].getvalue() == ""


def test_demultiplex_no_sample_sheet(tmp_dir):
    """
    Test :py:func:`riboviz.demultiplex_fastq.demultiplex` raises
    ``FileNotFoundError`` if the sample sheet is not found.

    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    """
    with pytest.raises(FileNotFoundError):
        demultiplex_fastq.demultiplex(
            "nosuchfile.tsv",
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex.fastq"),
            out_dir=tmp_dir)


def test_demultiplex_no_read1_file(tmp_dir):
    """
    Test :py:func:`riboviz.demultiplex_fastq.demultiplex` raises
    ``FileNotFoundError`` if the FASTQ file is not found.

    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    """
    with pytest.raises(FileNotFoundError):
        demultiplex_fastq.demultiplex(
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex_barcodes.tsv"),
            "nosuchfile.fastq",
            out_dir=tmp_dir)


def test_demultiplex_no_read2_file(tmp_dir):
    """
    Test :py:func:`riboviz.demultiplex_fastq.demultiplex` raises
    ``FileNotFoundError`` if the paired FASTQ file is not found.

    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    """
    with pytest.raises(FileNotFoundError):
        demultiplex_fastq.demultiplex(
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex_barcodes.tsv"),
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex.fastq"),
            "nosuchfile.fastq",
            out_dir=tmp_dir)


def test_demultiplex_output_error():
    """
    Test :py:func:`riboviz.demultiplex_fastq.demultiplex` raises
    ``IOError`` if the output directory cannot be created.
    """
    with pytest.raises(IOError):
        demultiplex_fastq.demultiplex(
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex_barcodes.tsv"),
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex.fastq"),
            # Pass existing file as out_dir value.
            out_dir=os.path.join(riboviz.test.SIMDATA_DIR,
                                 "multiplex_barcodes.tsv"))


@pytest.mark.parametrize("file_format",
                         [fastq.FASTQ_FORMAT,
                          fastq.FQ_FORMAT,
                          fastq.FASTQ_FORMAT.upper(),
                          fastq.FQ_FORMAT.upper()])
def test_demultiplex(tmp_dir, file_format):
    """
    Test :py:func:`riboviz.demultiplex_fastq.demultiplex`.

    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    :param file_format: FASTQ file format
    :type file_format: str or unicode
    """
    tmp_fastq_file = os.path.join(tmp_dir,
                                  file_format.format("test_multiplex"))
    shutil.copyfile(os.path.join(riboviz.test.SIMDATA_DIR,
                                 "multiplex.fastq"),
                    tmp_fastq_file)
    demultiplex_fastq.demultiplex(
        os.path.join(riboviz.test.SIMDATA_DIR,
                     "multiplex_barcodes.tsv"),
        tmp_fastq_file,
        mismatches=2,
        out_dir=tmp_dir)

    actual_num_reads = os.path.join(
        tmp_dir,
        demultiplex_fastq.NUM_READS_FILE)
    expected_num_reads = os.path.join(
        riboviz.test.SIMDATA_DIR,
        "deplex",
        demultiplex_fastq.NUM_READS_FILE)
    utils.equal_tsv(expected_num_reads, actual_num_reads)
    for tag in ["Tag0", "Tag1", "Tag2", "Unassigned"]:
        # Actual data has extension matching lower-case version
        # of multiplexed file's extension.
        actual_fq = os.path.join(tmp_dir,
                                 file_format.lower().format(tag))
        # Simulated data always has a .fastq extension.
        expected_fq = os.path.join(riboviz.test.SIMDATA_DIR,
                                   "deplex",
                                   fastq.FASTQ_FORMAT.format(tag))
        fastq.equal_fastq(expected_fq, actual_fq)


@pytest.mark.parametrize("file_format",
                         [(fastq.FASTQ_GZ_FORMAT,
                           fastq.FASTQ_FORMAT),
                          (fastq.FQ_GZ_FORMAT,
                           fastq.FQ_FORMAT),
                          (fastq.FASTQ_GZ_FORMAT.upper(),
                           fastq.FASTQ_FORMAT),
                          (fastq.FQ_GZ_FORMAT.upper(),
                           fastq.FQ_FORMAT)])
def test_demultiplex_gz(tmp_dir, file_format):
    """
    Test :py:func:`riboviz.demultiplex_fastq.demultiplex` using
    GZIPped FASTQ files.

    Each ``file_format`` consists of a FASTQ GZIP file name format and
    the corresponding non-GZIP FASTQ file name format.

    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    :param file_format: File name format
    :type file_format: tuple(str or unicode, str or unicode)
    """
    gz_fmt, fmt = file_format
    tmp_fastq_file = os.path.join(tmp_dir,
                                  gz_fmt.format("test_multiplex"))
    with open(os.path.join(riboviz.test.SIMDATA_DIR,
                           "multiplex.fastq"), "rb") as fr:
        with gzip.open(tmp_fastq_file, "wb") as fw:
            shutil.copyfileobj(fr, fw)
    demultiplex_fastq.demultiplex(
        os.path.join(riboviz.test.SIMDATA_DIR,
                     "multiplex_barcodes.tsv"),
        tmp_fastq_file,
        mismatches=2,
        out_dir=tmp_dir)

    actual_num_reads = os.path.join(
        tmp_dir,
        demultiplex_fastq.NUM_READS_FILE)
    expected_num_reads = os.path.join(
        riboviz.test.SIMDATA_DIR,
        "deplex",
        demultiplex_fastq.NUM_READS_FILE)
    utils.equal_tsv(expected_num_reads, actual_num_reads)
    for tag in ["Tag0", "Tag1", "Tag2", "Unassigned"]:
        # Actual data has extension matching lower-case version
        # of multiplexed file's extension.
        actual_fq_gz = os.path.join(tmp_dir,
                                    gz_fmt.lower().format(tag))
        actual_fq = os.path.join(tmp_dir,
                                 fmt.format(tag))
        # Simulated data always has a .fastq extension.
        expected_fq = os.path.join(riboviz.test.SIMDATA_DIR,
                                   "deplex",
                                   fastq.FASTQ_FORMAT.format(tag))
        # Decompress actual_fq_gz
        with gzip.open(actual_fq_gz, "rb") as fr:
            with open(actual_fq, "wb") as fw:
                shutil.copyfileobj(fr, fw)
        fastq.equal_fastq(expected_fq, actual_fq)
