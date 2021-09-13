"""
:py:mod:`riboviz.sample_sheets` tests.
"""
import csv
import os
import tempfile
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from riboviz import sample_sheets
from riboviz.test import data

TEST_SAMPLES = [("Tag0", "ACG", 27),
                ("Tag1", "GAC", 27),
                ("Tag2", "CGA", 27),
                ("Tag3", "CCC", 0)]
"""
Test sample sheet data.
"""
TEST_SAMPLES_DEPLEX = TEST_SAMPLES + \
                      [(sample_sheets.UNASSIGNED_TAG,
                        sample_sheets.UNASSIGNED_READ, 9),
                       (sample_sheets.TOTAL_READS, "", 90)]
"""
Test sample sheet data with information on demultiplexed samples
including unassigned samples and total number of reads.
"""


def list_to_df(data, ignore_num_reads=False):
    """
    Convert a list of (sample_id, tag_read, num_reads) to a data
    frame.

    :param data: List
    :type data: list(tuple(str or unicode, str or unicode, int))
    :param ignore_num_reads: Ignore num_reads when creating data \
    frame?
    :type ignore_num_reads: bol
    :return: Data frame
    :rtype: pandas.core.frame.DataFrame
    """
    df_data = []
    for (sample_id, tag_read, num_reads) in data:
        entry = {sample_sheets.SAMPLE_ID: sample_id,
                 sample_sheets.TAG_READ: tag_read}
        if not ignore_num_reads:
            entry[sample_sheets.NUM_READS] = num_reads
        df_data.append(entry)
    df = pd.DataFrame(df_data)
    return df


@pytest.fixture(scope="function")
def tmp_file():
    """
    Create a temporary file.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp")
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


@pytest.mark.parametrize("file_name",
                         ["no_such_file.txt",
                          os.path.dirname(data.__file__)])
def test_load_sample_sheet_no_such_file(file_name):
    """
    Test :py:func:`riboviz.sample_sheets.load_sample_sheet` raises
    ``FileNotFoundError`` if the sample sheet does not exist or is a
    directory.

    :param file_name: File name
    :type file_name: str or unicode
    """
    with pytest.raises(FileNotFoundError):
        sample_sheets.load_sample_sheet(file_name)


@pytest.mark.parametrize("header",
                         [["A", "B"],
                          [sample_sheets.SAMPLE_ID, "A"],
                          [sample_sheets.TAG_READ, "A"]])
def test_load_sample_sheet_non_sample_sheet(tmp_file, header):
    """
    Test :py:func:`riboviz.sample_sheets.load_sample_sheet` raises
    an ``AssertionError`` if the sample sheet does not contain the
    columns :py:const:`sample_sheets.SAMPLE_ID` and
    :py:const:`sample_sheets.TAG_READ`.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    :param header: Header column names
    :type header: str or unicode
    """
    with open(tmp_file, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file,
                                delimiter="\t",
                                fieldnames=header)
        writer.writeheader()
    with pytest.raises(AssertionError):
        sample_sheets.load_sample_sheet(tmp_file)


def test_load_sample_sheet(tmp_file):
    """
    Test :py:func:`riboviz.sample_sheets.load_sample_sheet` with
    a valid sample sheet file.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    list_to_df(TEST_SAMPLES, True).to_csv(
        tmp_file, mode='w', sep="\t", index=False)
    df = sample_sheets.load_sample_sheet(tmp_file,
                                         comment="#")
    for column in [sample_sheets.SAMPLE_ID, sample_sheets.TAG_READ]:
        assert column in df.columns,\
            "Missing column {}".format(column)
    # Expect data frame content to be equivalent to that which was
    # saved.
    assert_frame_equal(df, list_to_df(TEST_SAMPLES, True))


def test_load_deplexed_sample_sheet(tmp_file):
    """
    Test :py:func:`riboviz.sample_sheets.load_deplexed_sample_sheet`
    with a valid sample sheet file with information on demultiplexed
    samples.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    list_to_df(TEST_SAMPLES_DEPLEX, False).to_csv(
        tmp_file, mode='w', sep="\t", index=False)
    df = sample_sheets.load_deplexed_sample_sheet(
        tmp_file, comment="#").fillna("")
    for column in [sample_sheets.SAMPLE_ID,
                   sample_sheets.TAG_READ,
                   sample_sheets.NUM_READS]:
        assert column in df.columns,\
            "Missing column {}".format(column)
    # Expect data frame content to be equivalent to that which was
    # saved.
    assert_frame_equal(df, list_to_df(TEST_SAMPLES_DEPLEX, False))


def test_save_deplexed_sample_sheet(tmp_file):
    """
    Test :py:func:`riboviz.sample_sheets.save_deplexed_sample_sheet`
    saves a sample sheet with the expected content.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    df = list_to_df(TEST_SAMPLES, False)
    # Save sample sheet with additional information about
    # demultiplexed samples then reload.
    sample_sheets.save_deplexed_sample_sheet(df, 9, tmp_file)
    df = sample_sheets.load_deplexed_sample_sheet(tmp_file,
                                                  comment="#").fillna("")
    # Validate against expected data frame.
    check_df = list_to_df(TEST_SAMPLES_DEPLEX, False)
    assert_frame_equal(df, check_df)


def test_get_non_zero_deplexed_samples():
    """
    Test :py:func:`riboviz.sample_sheets.get_non_zero_deplexed_samples`.
    """
    data = [("A", "CTAGA", 36067201),
            ("B", "GATCA", 48718085),
            ("C", "GCATA", 0),
            ("D", "TAGAC", 30848592),
            ("E", "TCTAG", 41407639),
            ("F", "ACTGA", 0),
            ("G", "CAGTA", 30643912),
            ("H", "GTACT", 28091942),
            ("I", "TGCAT", 0),
            ("Unassigned", "NNNNNNNNN", 8984320),
            ("Total", "", 264283118)]
    df = list_to_df(data, False)
    samples = sample_sheets.get_non_zero_deplexed_samples(df)
    expected_samples = [sample_id for sample_id, _, num_reads
                        in data[0:-2] if num_reads > 0]
    assert expected_samples == samples


def test_get_non_zero_deplexed_samples_all_zero():
    """
    Test :py:func:`riboviz.sample_sheets.get_non_zero_deplexed_samples`
    with sample sheet data for which the number of reads for each
    sample is zero.
    """
    data = [("A", "CTAGA", 0),
            ("B", "GATCA", 0),
            ("C", "GCATA", 0),
            ("Unassigned", "NNNNNNNNN", 12345),
            ("Total", "", 12345)]
    df = list_to_df(data, False)
    samples = sample_sheets.get_non_zero_deplexed_samples(df)
    assert len(samples) == 0
