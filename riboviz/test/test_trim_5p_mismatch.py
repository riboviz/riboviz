"""
riboviz.trim_5p_mismatch test suite.
"""
import os
import tempfile
import pytest
import pandas as pd
import riboviz
from riboviz.test import data
from riboviz import trim_5p_mismatch


TEST_5P_FILE = "trim_5p_mismatch.sam"
""" 5p test file name, assumed to be in riboviz.test.data """
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
Expected results for testing TEST_5P_FILE with 0, 1 and 2
mismatches. Format is of list(tuple(max_mismatches, dict with expected
results from trim_5p_mismatch(TEST_5P_FILE, _, max_mismatches)
"""
TEST_5P_CASES = list(zip([TEST_5P_FILE] * len(TEST_5P_EXPECTED),
                         TEST_5P_EXPECTED))
""" Test cases for TEST_5P_FILE. """

TEST_5POS5NEG_FILE = "trim_5pos5neg.sam"
""" 5p test file name, assumed to be in riboviz.test.data """
TEST_5POS5NEG_EXPECTED = [(0, {trim_5p_mismatch.NUM_PROCESSED: 10,
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
Expected results for testing TEST_5POS_NEG_FILE with 0, 1 and 2
mismatches. Format is of list(tuple(max_mismatches, dict with expected
results from trim_5p_mismatch(TEST_5POS_NEG_FILE, _, max_mismatches)
"""
TEST_5POS5NEG_CASES = list(zip(
    [TEST_5POS5NEG_FILE] * len(TEST_5POS5NEG_EXPECTED),
    TEST_5POS5NEG_EXPECTED))
""" Test cases for TEST_5POS_NEG_FILE. """


@pytest.fixture(scope="function")
def trimmed_sam_file():
    """
    Create a temporary SAM file to write trimmed SAM file to.

    :return: path to SAM file
    :rtype: str or unicode
    """
    _, sam_file = tempfile.mkstemp(prefix="tmp", suffix=".sam")
    yield sam_file
    if os.path.exists(sam_file):
        os.remove(sam_file)


@pytest.fixture(scope="function")
def summary_file():
    """
    Create a temporary TSV summary file to write summary data to.

    :return: path to summary file
    :rtype: str or unicode
    """
    _, summary_file = tempfile.mkstemp(prefix="tmp", suffix=".tsv")
    yield summary_file
    if os.path.exists(summary_file):
        os.remove(summary_file)


@pytest.mark.parametrize("test_case", TEST_5P_CASES + TEST_5POS5NEG_CASES)
def test_trim_5p_mismatch(test_case, trimmed_sam_file):
    """
    Run trim_5p_mismatch and validate summary returned.

    :param test_case: Test case with test SAM file name, maximum
    number of mismatches to test and dictionary with expected results
    from trim_5p_mismatch.
    :type test_case: tuple(str or unicode, tuple(int, dict))
    :param trimmed_sam_file: path to trimmed SAM file
    :type trimmed_sam_file: str or unicode
    """
    sam_file_name, (max_mismatches, expected_summary) = test_case
    sam_file = os.path.join(os.path.dirname(data.__file__), sam_file_name)
    summary = trim_5p_mismatch.trim_5p_mismatch(sam_file,
                                                trimmed_sam_file,
                                                True,
                                                max_mismatches)
    assert summary == expected_summary, "Unexpeted summary"


@pytest.mark.parametrize("test_case", TEST_5P_CASES + TEST_5POS5NEG_CASES)
def test_trim_5p_mismatch_file(test_case, trimmed_sam_file, summary_file):
    """
    Run trim_5p_mismatch_file and validate summary file created.

    :param test_case: Test case with test SAM file name, maximum
    number of mismatches to test and dictionary with expected results
    from trim_5p_mismatch.
    :type test_case: tuple(str or unicode, tuple(int, dict))
    :param trimmed_sam_file: path to trimmed SAM file
    :type trimmed_sam_file: str or unicode
    :param summary_file: path to TSV summary file
    :type summary_file: str or unicode
    """
    sam_file_name, (max_mismatches, expected_summary) = test_case
    sam_file = os.path.join(os.path.dirname(data.__file__), sam_file_name)
    trim_5p_mismatch.trim_5p_mismatch_file(sam_file,
                                           trimmed_sam_file,
                                           True,
                                           max_mismatches,
                                           summary_file)
    summary_df = pd.read_csv(summary_file, sep="\t", comment="#")
    summary = summary_df.to_dict('records')
    assert len(summary_df) == 1, "Expected 1 summary row only"
    assert summary[0] == expected_summary, "Unexpeted summary"
