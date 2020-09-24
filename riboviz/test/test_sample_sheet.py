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


@pytest.fixture(scope="function")
def sample_sheet_file():
    """
    Create a temporary sample sheet file using data from
    :py:const:`TEST_SAMPLES`.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp")
    with open(tmp_file, 'w') as csv_file:
        writer = csv.DictWriter(csv_file,
                                delimiter="\t",
                                fieldnames=[sample_sheets.SAMPLE_ID,
                                            sample_sheets.TAG_READ])
        writer.writeheader()
        for (sample_id, tag_read, _) in TEST_SAMPLES:
                writer.writerow(
                    { sample_sheets.SAMPLE_ID: sample_id,
                      sample_sheets.TAG_READ: tag_read })
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


@pytest.fixture(scope="function")
def deplex_sample_sheet_file():
    """
    Create a temporary sample sheet file with information on
    demultiplexed samples using data from
    :py:const:`TEST_SAMPLES_DEPLEX`.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp")
    with open(tmp_file, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file,
                                delimiter="\t",
                                fieldnames=[sample_sheets.SAMPLE_ID,
                                            sample_sheets.TAG_READ,
                                            sample_sheets.NUM_READS])
        writer.writeheader()
        for (sample_id, tag_read, num_reads) in TEST_SAMPLES_DEPLEX:
                writer.writerow(
                    { sample_sheets.SAMPLE_ID: sample_id,
                      sample_sheets.TAG_READ: tag_read,
                      sample_sheets.NUM_READS: num_reads })
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


def test_load_sample_sheet(sample_sheet_file):
    """
    Test :py:func:`riboviz.sample_sheets.load_sample_sheet` with
    a valid sample sheet file.

    :param sample_sheet_file: Sample sheet file
    :type sample_sheet_file: str or unicode
    """
    df = sample_sheets.load_sample_sheet(sample_sheet_file,
                                         comment="#")
    for column in [sample_sheets.SAMPLE_ID, sample_sheets.TAG_READ]:
        assert column in df.columns,\
            "Missing column {}".format(column)
    check_samples = TEST_SAMPLES
    num_rows, _ = df.shape
    assert num_rows == len(check_samples), \
        "Unexpected number of rows"
    # Validate rows
    for sample_id, tag_read, _ in check_samples:
        sample_df = df[(df[sample_sheets.SAMPLE_ID] == sample_id) &
                       (df[sample_sheets.TAG_READ] == tag_read)]
        assert not sample_df.empty, "Missing row for {}".format(sample_id)


def test_load_deplexed_sample_sheet(deplex_sample_sheet_file):
    """
    Test :py:func:`riboviz.sample_sheets.load_deplexed_sample_sheet`
    with a valid sample sheet file with information on demultiplexed
    samples.

    :param deplex_sample_sheet_file: Sample sheet file
    :type deplex_sample_sheet_file: str or unicode
    """
    df = sample_sheets.load_deplexed_sample_sheet(
        deplex_sample_sheet_file,
        comment="#").fillna("")
    for column in [sample_sheets.SAMPLE_ID,
                   sample_sheets.TAG_READ,
                   sample_sheets.NUM_READS]:
        assert column in df.columns,\
            "Missing column {}".format(column)
    num_rows, _ = df.shape
    assert num_rows == len(samples), \
        "Unexpected number of rows"
    for sample_id, tag_read, num_reads in samples:
        sample_df = df[(df[sample_sheets.SAMPLE_ID] == sample_id) &
                       (df[sample_sheets.TAG_READ] == tag_read) &
                       (df[sample_sheets.NUM_READS] == num_reads)]
        assert not sample_df.empty, "Missing row for {}".format(sample_id)



def test_save_deplexed_sample_sheet(tmp_file,
                                    sample_sheet_file,
                                    deplex_sample_sheet_file):
    """
    Test :py:func:`riboviz.sample_sheets.save_deplexed_sample_sheet`
    saves a sample sheet with the expected content validating it
    against a valid sample sheet file with information on
    demultiplexed samples.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    :param sample_sheet_file: Sample sheet file
    :type sample_sheet_file: str or unicode
    :param deplex_sample_sheet_file: Sample sheet file
    :type deplex_sample_sheet_file: str or unicode
    """
    # Load sample sheet.
    df = sample_sheets.load_sample_sheet(sample_sheet_file,
                                         comment="#")
    # Add number of reads for each sample.
    num_reads_list = [num_reads for _, _, num_reads in TEST_SAMPLES]
    df[sample_sheets.NUM_READS] = num_reads_list
    # Save sample sheet with additional information about
    # demultiplexed samples then reload.
    sample_sheets.save_deplexed_sample_sheet(df, 9, tmp_file)
    df = sample_sheets.load_deplexed_sample_sheet(tmp_file,
                                                  comment="#").fillna("")
    # Load sample sheet with expected content and compare.
    check_df = sample_sheets.load_deplexed_sample_sheet(
        deplex_sample_sheet_file,
        comment="#").fillna("")
    assert_frame_equal(df, check_df)


def test_get_non_zero_deplexed_samples():
    """
    Test :py:func:`riboviz.sample_sheets.get_non_zero_deplexed_samples`
    saves a sample sheet with the expected content.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    # def get_non_zero_deplexed_samples(sample_sheet):
    # The sample sheet is assumed to have columns, ``SampleID``,
    # ``TagRead`` and ``NumReads``. Rows whose ``SampleID`` values are
    # ``Unassigned`` or ``Total`` are ignored.
    #non_zero_samples = sample_sheet[
    #    ~sample_sheet[SAMPLE_ID].isin([UNASSIGNED_TAG, TOTAL_READS])
    #    & sample_sheet[NUM_READS] != 0]
    #return list(non_zero_samples[SAMPLE_ID])
    # TODO
    pass
