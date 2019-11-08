"""
demultiplex_fastq.py test suite.
"""
from io import StringIO
from contextlib import ExitStack
from riboviz.tools.demultiplex_fastq import assign_sample
from riboviz.tools.demultiplex_fastq import assign_samples


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


def test_assign_sample():
    """
    Test assign_sample with matching barcode.
    """
    with StringIO() as read1_fh:
        barcode = "AAA"
        is_assigned = assign_sample(FASTQ_RECORD1, None,
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
        is_assigned = assign_sample(FASTQ_RECORD1, None,
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
        is_assigned = assign_sample(FASTQ_RECORD1, FASTQ_RECORD2,
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
        is_assigned = assign_sample(FASTQ_RECORD1, FASTQ_RECORD2,
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
        is_assigned = assign_samples(FASTQ_RECORD1, FASTQ_RECORD2,
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
        is_assigned = assign_samples(FASTQ_RECORD1, FASTQ_RECORD2,
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
