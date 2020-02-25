"""
:py:mod:`riboviz.sam_bam` tests.
"""
import os
import pytest
from riboviz import sam_bam
from riboviz.test import data


@pytest.mark.parametrize("extension",
                         [sam_bam.SAM_EXT, sam_bam.SAM_EXT.upper()])
def test_is_sam(extension):
    """
    Test :py:func:`riboviz.sam_bam.is_sam` with SAM extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert sam_bam.is_sam("sample.{}".format(extension))


@pytest.mark.parametrize("extension",
                         [sam_bam.BAM_EXT,
                          sam_bam.BAM_EXT.upper(),
                          "txt"])
def test_not_is_sam(extension):
    """
    Test :py:func:`riboviz.sam_bam.is_sam` with non-SAM extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert not sam_bam.is_sam("sample.{}".format(extension))


@pytest.mark.parametrize("extension",
                         [sam_bam.BAM_EXT, sam_bam.BAM_EXT.upper()])
def test_is_bam(extension):
    """
    Test :py:func:`riboviz.sam_bam.is_bam` with BAM extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert sam_bam.is_bam("sample.{}".format(extension))


@pytest.mark.parametrize("extension",
                         [sam_bam.SAM_EXT, sam_bam.SAM_EXT.upper(), "txt"])
def test_not_is_bam(extension):
    """
    Test :py:func:`riboviz.sam_bam.is_bam` with non-BAM extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert not sam_bam.is_bam("sample.{}".format(extension))


@pytest.mark.parametrize("test_case",
                         [("WTnone_rRNA_map_20", (20, 6)),
                          ("WTnone_rRNA_map_6_primary", (6, 6)),
                          ("WTnone_rRNA_map_14_secondary", (14, 0))])
@pytest.mark.parametrize("file_format",
                         [sam_bam.SAM_FORMAT, sam_bam.BAM_FORMAT])
def test_count_sequences(test_case, file_format):
    """
    Test :py:func:`riboviz.sam_bam.count_sequences`.

    Each ``test_case`` includes a SAM/BAM file name prefix (which
    will be converted into a file name using ``file_format``), and a
    tuple with the expected number of sequences and expected number of
    mapped (primary) sequences).

    :param test_case: Test case
    :type test_case: tuple(str or unicode, tuple(int, int))
    :param file_format: File name format
    :type file_format: str or unicode
    """
    file_name, expected_counts = test_case
    sam_bam_file = os.path.join(os.path.dirname(data.__file__),
                                file_format.format(file_name))
    actual_counts = sam_bam.count_sequences(sam_bam_file)
    assert expected_counts == actual_counts
