"""
riboviz.demultiplex_fastq test suite.
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
""" Sample fastq record """
FASTQ_RECORD2 = ["@X1:Tag_AAC_ 1:N:0:XXXXXXXX\n",
                 "AAAAAAAA\n",
                 "+\n",
                 "IIIIIIII\n"]
""" Sample fastq record """


@pytest.fixture(scope="function")
def temporary_dir():
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
    Test assign_sample with matching barcode.
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
    Test assign_sample with non-matching barcode.
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
    Test assign_sample with paired end.
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
    Test assign_sample with paired end and non-matching barcode.
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
    Test assign_samples with paired ends.
    """
    with ExitStack() as stack:
        read1_fhs = [stack.enter_context(StringIO()) for f in range(2)]
        read2_fhs = [stack.enter_context(StringIO()) for f in range(2)]
        barcodes = ["CCC", "AAA"]
        num_samples = 2
        num_reads = [0] * num_samples
        is_assigned = demultiplex_fastq.assign_samples(
            FASTQ_RECORD1, FASTQ_RECORD2,
            barcodes,
            read1_fhs, read2_fhs,
            True,
            num_samples,
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
    Test assign_samples with paired ends and non-matching barcode.
    """
    with ExitStack() as stack:
        read1_fhs = [stack.enter_context(StringIO()) for f in range(2)]
        read2_fhs = [stack.enter_context(StringIO()) for f in range(2)]
        barcodes = ["GGG", "TTT"]
        num_samples = 2
        num_reads = [0] * num_samples
        is_assigned = demultiplex_fastq.assign_samples(
            FASTQ_RECORD1, FASTQ_RECORD2,
            barcodes,
            read1_fhs, read2_fhs,
            True,
            num_samples,
            num_reads,
            1, "_")
        assert not is_assigned
        assert num_reads[0] == 0
        assert num_reads[1] == 0
        assert read1_fhs[0].getvalue() == ""
        assert read2_fhs[0].getvalue() == ""
        assert read1_fhs[1].getvalue() == ""
        assert read2_fhs[1].getvalue() == ""


def test_demultiplex_no_sample_sheet(temporary_dir):
    """
    Test demultiplex raises FileNotFoundError if the sample sheet is
    not found.

    :param temporary_dir: Temporary directory
    :type temporary_dir: str or unicode
    """
    with pytest.raises(FileNotFoundError):
        demultiplex_fastq.demultiplex(
            "nosuchfile.tsv",
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex.fastq"),
            out_dir=temporary_dir)


def test_demultiplex_no_read1_file(temporary_dir):
    """
    Test demultiplex raises FileNotFoundError if the
    read1 file is not found.

    :param temporary_dir: Temporary directory
    :type temporary_dir: str or unicode
    """
    with pytest.raises(FileNotFoundError):
        demultiplex_fastq.demultiplex(
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex_barcodes.tsv"),
            "nosuchfile.fastq",
            out_dir=temporary_dir)


def test_demultiplex_no_read2_file(temporary_dir):
    """
    Test demultiplex raises FileNotFoundError if the
    read2 file is not found.

    :param temporary_dir: Temporary directory
    :type temporary_dir: str or unicode
    """
    with pytest.raises(FileNotFoundError):
        demultiplex_fastq.demultiplex(
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex_barcodes.tsv"),
            os.path.join(riboviz.test.SIMDATA_DIR,
                         "multiplex.fastq"),
            "nosuchfile.fastq",
            out_dir=temporary_dir)


def test_demultiplex_output_error():
    """
    Test demultiplex raises FileNotFoundError if the
    output directory cannot be created.
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


@pytest.mark.parametrize("fmt",
                         [fastq.FASTQ_FORMAT,
                          fastq.FQ_FORMAT,
                          fastq.FASTQ_FORMAT.upper(),
                          fastq.FQ_FORMAT.upper()])
def test_demultiplex(temporary_dir, fmt):
    """
    Validate that TSV and FASTQ files produced by
    riboviz.demultiplex_fastq have the expected
    content.

    :param temporary_dir: Temporary directory
    :type temporary_dir: str or unicode
    :param fmt: fastq file format
    :type fmt: str or unicode
    """
    tmp_fastq_file = os.path.join(temporary_dir,
                                  fmt.format("test_multiplex"))
    shutil.copyfile(os.path.join(riboviz.test.SIMDATA_DIR,
                                 "multiplex.fastq"),
                    tmp_fastq_file)
    demultiplex_fastq.demultiplex(
        os.path.join(riboviz.test.SIMDATA_DIR,
                     "multiplex_barcodes.tsv"),
        tmp_fastq_file,
        mismatches=2,
        out_dir=temporary_dir)

    actual_num_reads = os.path.join(
        temporary_dir,
        demultiplex_fastq.NUM_READS_FILE)
    expected_num_reads = os.path.join(
        riboviz.test.SIMDATA_DIR,
        "deplex",
        demultiplex_fastq.NUM_READS_FILE)
    utils.equal_tsv_files(expected_num_reads, actual_num_reads)
    for tag in ["Tag0", "Tag1", "Tag2", "Unassigned"]:
        # Actual data has extension matching lower-case version
        # of multiplexed file's extension.
        actual_fq = os.path.join(temporary_dir,
                                 fmt.lower().format(tag))
        # Simulated data always has a .fastq extension.
        expected_fq = os.path.join(riboviz.test.SIMDATA_DIR,
                                   "deplex",
                                   fastq.FASTQ_FORMAT.format(tag))
        fastq.equal_fastq(expected_fq, actual_fq)


@pytest.mark.parametrize("extension",
                         [(fastq.FASTQ_GZ_FORMAT,
                           fastq.FASTQ_FORMAT),
                          (fastq.FQ_GZ_FORMAT,
                           fastq.FQ_FORMAT),
                          (fastq.FASTQ_GZ_FORMAT.upper(),
                           fastq.FASTQ_FORMAT),
                          (fastq.FQ_GZ_FORMAT.upper(),
                           fastq.FQ_FORMAT)])
def test_demultiplex_gz(temporary_dir, extension):
    """
    Validate that TSV and FASTQ files produced by
    riboviz.demultiplex_fastq on multiplexed GZ files have the
    expected content.

    :param temporary_dir: Temporary directory
    :type temporary_dir: str or unicode
    :param extension: (gz format,  non-gz format)
    :type extension: tuple(str or unicode, str or unicode)
    """
    gz_fmt, fmt = extension
    tmp_fastq_file = os.path.join(temporary_dir,
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
        out_dir=temporary_dir)

    actual_num_reads = os.path.join(
        temporary_dir,
        demultiplex_fastq.NUM_READS_FILE)
    expected_num_reads = os.path.join(
        riboviz.test.SIMDATA_DIR,
        "deplex",
        demultiplex_fastq.NUM_READS_FILE)
    utils.equal_tsv_files(expected_num_reads, actual_num_reads)
    for tag in ["Tag0", "Tag1", "Tag2", "Unassigned"]:
        # Actual data has extension matching lower-case version
        # of multiplexed file's extension.
        actual_fq_gz = os.path.join(temporary_dir,
                                    gz_fmt.lower().format(tag))
        actual_fq = os.path.join(temporary_dir,
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
