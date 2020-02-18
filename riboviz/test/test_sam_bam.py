"""
riboviz.sam_bam tests.
"""
import os
import pytest
from riboviz import sam_bam
from riboviz.test import data


@pytest.mark.parametrize("extension",
                         [sam_bam.SAM_EXT, sam_bam.SAM_EXT.upper()])
def test_is_sam(extension):
    """
    Test is_sam with SAM extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert sam_bam.is_sam("sample.{}".format(extension))


@pytest.mark.parametrize("extension",
                         [sam_bam.BAM_EXT, sam_bam.BAM_EXT.upper(), "txt"])
def test_not_is_sam(extension):
    """
    Test is_sam with non-SAM extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert not sam_bam.is_sam("sample.{}".format(extension))


@pytest.mark.parametrize("extension",
                         [sam_bam.BAM_EXT, sam_bam.BAM_EXT.upper()])
def test_is_bam(extension):
    """
    Test is_bam with BAM extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert sam_bam.is_bam("sample.{}".format(extension))


@pytest.mark.parametrize("extension",
                         [sam_bam.SAM_EXT, sam_bam.SAM_EXT.upper(), "txt"])
def test_not_is_bam(extension):
    """
    Test is_sam with non-BAM extensions.

    :param extension: Extension
    :type extension: str or unicode
    """
    assert not sam_bam.is_bam("sample.{}".format(extension))


@pytest.mark.parametrize("test_case",
                         [("WTnone_rRNA_map_20", (20, 6)),
                          ("WTnone_rRNA_map_6_primary", (6, 6)),
                          ("WTnone_rRNA_map_14_secondary", (14, 0))])
@pytest.mark.parametrize("file_type",
                         [sam_bam.SAM_FORMAT, sam_bam.BAM_FORMAT])
def test_count_sequences(test_case, file_type):
    """
    Test test_count_sequences.

    :param test_case: Test case with test SAM/BAM file name prefix \
    (to which ".file_type" will be appended, and a tuple with the \
    expected number of sequences and expected number of mapped \
    (primary) sequences).
    :type test_case: tuple(str or unicode, tuple(int, int))
    :param file_type: SAM or BAM
    :type file_type: str or unicode
    """
    file_name, expected_counts = test_case
    sam_bam_file = os.path.join(os.path.dirname(data.__file__),
                                file_type.format(file_name))
    actual_counts = sam_bam.count_sequences(sam_bam_file)
    assert expected_counts == actual_counts
