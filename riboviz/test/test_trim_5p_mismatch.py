"""
:py:mod:`riboviz.trim_5p_mismatch` tests.
"""
import os
import tempfile
import pytest
import pandas as pd
from riboviz.test import data
from riboviz import sam_bam
from riboviz import trim_5p_mismatch


TEST_5P_FILE = "trim_5p_mismatch.sam"
""" Test file in :py:mod:`riboviz.test.data`. """
TEST_5P_EXPECTED = [(0, {trim_5p_mismatch.NUM_PROCESSED: 13,
                         trim_5p_mismatch.NUM_DISCARDED: 7,
                         trim_5p_mismatch.NUM_TRIMMED: 4,
                         trim_5p_mismatch.NUM_WRITTEN: 6}),
                    (1, {trim_5p_mismatch.NUM_PROCESSED: 13,
                         trim_5p_mismatch.NUM_DISCARDED: 4,
                         trim_5p_mismatch.NUM_TRIMMED: 4,
                         trim_5p_mismatch.NUM_WRITTEN: 9}),
                    (2, {trim_5p_mismatch.NUM_PROCESSED: 13,
                         trim_5p_mismatch.NUM_DISCARDED: 2,
                         trim_5p_mismatch.NUM_TRIMMED: 4,
                         trim_5p_mismatch.NUM_WRITTEN: 11})]
"""
Expected results for testing :py:const:`TEST_5P_FILE` with 0, 1 and 2
mismatches. Format is of a list of tuples. Each tuple holds:

* Number of mismatches for :py:mod:`riboviz.trim_5p_mismatch` function
* calls.
* Dictionary with expected trimming summary resulting from
  :py:mod:`riboviz.trim_5p_mismatch` function calls.
"""
TEST_5P_CASES = list(zip([TEST_5P_FILE] * len(TEST_5P_EXPECTED),
                         TEST_5P_EXPECTED))
"""
List of tuples of form (:py:const:`TEST_5P_FILE`, tuple from
:py:const:`TEST_5P_CASES`.
"""

TEST_5POS_5NEG_FILE = "trim_5pos5neg.sam"
""" Test file in :py:mod:`riboviz.test.data`. """
TEST_5POS_5NEG_EXPECTED = [(0, {trim_5p_mismatch.NUM_PROCESSED: 10,
                                trim_5p_mismatch.NUM_DISCARDED: 0,
                                trim_5p_mismatch.NUM_TRIMMED: 10,
                                trim_5p_mismatch.NUM_WRITTEN: 10}),
                           (1, {trim_5p_mismatch.NUM_PROCESSED: 10,
                                trim_5p_mismatch.NUM_DISCARDED: 0,
                                trim_5p_mismatch.NUM_TRIMMED: 10,
                                trim_5p_mismatch.NUM_WRITTEN: 10}),
                           (2, {trim_5p_mismatch.NUM_PROCESSED: 10,
                                trim_5p_mismatch.NUM_DISCARDED: 0,
                                trim_5p_mismatch.NUM_TRIMMED: 10,
                                trim_5p_mismatch.NUM_WRITTEN: 10})]
"""
Expected results for testing :py:const:`TEST_5POS_5NEG_FILE` with 0, 1
and 2 mismatches. Format is of a list of tuples. Each tuple holds:

* Number of mismatches for :py:mod:`riboviz.trim_5p_mismatch` function
* calls.
* Dictionary with expected trimming summary resulting from
  :py:mod:`riboviz.trim_5p_mismatch` function calls.
"""
TEST_5POS_5NEG_CASES = list(zip(
    [TEST_5POS_5NEG_FILE] * len(TEST_5POS_5NEG_EXPECTED),
    TEST_5POS_5NEG_EXPECTED))
"""
List of tuples of form (:py:const:`TEST_5POS_5NEG_FILE`, tuple from
:py:const:`TEST_5POS_5NEG_CASES`.
"""


@pytest.fixture(scope="function")
def tmp_sam_file():
    """
    Create a temporary file with a ``sam`` extension.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_sam_file = tempfile.mkstemp(
        prefix="tmp", suffix="." + sam_bam.SAM_EXT)
    yield tmp_sam_file
    if os.path.exists(tmp_sam_file):
        os.remove(tmp_sam_file)


@pytest.fixture(scope="function")
def tmp_tsv_file():
    """
    Create a temporary file with a ``tsv`` extension.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_tsv_file = tempfile.mkstemp(prefix="tmp", suffix=".tsv")
    yield tmp_tsv_file
    if os.path.exists(tmp_tsv_file):
        os.remove(tmp_tsv_file)


@pytest.mark.parametrize("test_case", TEST_5P_CASES + TEST_5POS_5NEG_CASES)
def test_trim_5p_mismatch(test_case, tmp_sam_file):
    """
    Run :py:func:`riboviz.trim_5p_mismatch.trim_5p_mismatch`.

    Each test case is a tuple with a SAM file, and another tuple. The
    other tuple holds:

    * Number of mismatches for
      :py:func:`riboviz.trim_5p_mismatch.trim_5p_mismatch` call.
    * Dictionary with expected trimming summary resulting from
      :py:func:`riboviz.trim_5p_mismatch.trim_5p_mismatch` call.

    :param test_case: Test case
    :type test_case: tuple(str or unicode, tuple(int, dict))
    :param tmp_sam_file: path to temporary file
    :type tmp_sam_file: str or unicode
    """
    sam_file_name, (max_mismatches, expected_summary) = test_case
    sam_file = os.path.join(os.path.dirname(data.__file__), sam_file_name)
    summary = trim_5p_mismatch.trim_5p_mismatch(sam_file,
                                                tmp_sam_file,
                                                True,
                                                max_mismatches)
    assert summary == expected_summary, "Unexpeted summary"


@pytest.mark.parametrize("test_case", TEST_5P_CASES + TEST_5POS_5NEG_CASES)
def test_trim_5p_mismatch_file(test_case, tmp_sam_file, tmp_tsv_file):
    """
    Run :py:func:`riboviz.trim_5p_mismatch.trim_5p_mismatch_file`.

    Each test case is a tuple with a SAM file, and another tuple. The
    other tuple holds:

    * Number of mismatches for
      :py:func:`riboviz.trim_5p_mismatch.trim_5p_mismatch_file` call.
    * Dictionary with expected trimming summary resulting from
      :py:func:`riboviz.trim_5p_mismatch.trim_5p_mismatch_file` call.

    :param test_case: Test case
    :type test_case: tuple(str or unicode, tuple(int, dict))
    :param tmp_sam_file: path to temporary file
    :type tmp_sam_file: str or unicode
    :param tmp_tsv_file: path to temporary file
    :type tmp_tsv_file: str or unicode
    """
    sam_file_name, (max_mismatches, expected_summary) = test_case
    sam_file = os.path.join(os.path.dirname(data.__file__), sam_file_name)
    trim_5p_mismatch.trim_5p_mismatch_file(sam_file,
                                           tmp_sam_file,
                                           True,
                                           max_mismatches,
                                           tmp_tsv_file)
    summary_df = pd.read_csv(tmp_tsv_file, sep="\t", comment="#")
    summary = summary_df.to_dict('records')
    assert len(summary_df) == 1, "Expected 1 summary row only"
    assert summary[0] == expected_summary, "Unexpeted summary"
