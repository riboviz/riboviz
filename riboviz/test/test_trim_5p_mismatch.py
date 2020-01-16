"""
riboviz.trim_5p_mismatch test suite.
"""
import os
import tempfile
import pytest
from riboviz import trim_5p_mismatch
from riboviz.test import data


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


@pytest.mark.parametrize("test_case", TEST_5P_CASES + TEST_5POS5NEG_CASES)
def test_trim_5p_mismatch(test_case, trimmed_sam_file):
    """
    Run upgrade_config_file on previous versions of the simulated UMI
    data configuration file and compare to the current version.

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
