"""
riboviz.workflow_record test suite.
"""
import os
import shutil
import tempfile
import pytest
import pandas as pd
from riboviz import workflow_record


NUM_DIRECTORIES = 5
""" Number of test sub-directories """
NUM_FILES = 3
""" Number of test files """
DIRECTORY_FORMAT = "dir{}"
""" Directory name format """
FILE_FORMAT = "file{}_{}.txt"
""" File name format """
PROGRAM_FORMAT = "program{}"
""" Program name format """
DESCRIPTION_FORMAT = "description{}"
""" Description format """
SAMPLE_FORMAT = "sample{}"
""" Sample name format """


@pytest.fixture(scope="function")
def tsv_file():
    """
    Create a temporary TSV file to write the worfklow record to.

    :return: path to TSV file
    :rtype: str or unicode
    """
    _, tsv_file = tempfile.mkstemp(prefix="tmp", suffix=".tsv")
    yield tsv_file
    if os.path.exists(tsv_file):
        os.remove(tsv_file)


def record_test_files(directory, workflow_record_file):
    """
    Populate workflow record file with entries corresponding to
    directory and file structure created by tmp_test_dir.

    :param directory: Temporary directory with subdirectories and
    files
    :type directory: str or unicode
    :param workflow_record_file: path to TSV file
    :type workflow_record_file: str or unicode
    """
    workflow_record.create_record_file(workflow_record_file, delimiter="\t")
    for d in range(0, NUM_DIRECTORIES):
        sub_dir = os.path.join(directory,
                               DIRECTORY_FORMAT.format(d))
        files = [os.path.join(sub_dir, FILE_FORMAT.format(d, f))
                 for f in range(0, NUM_FILES)]
        workflow_record.record_step(workflow_record_file,
                                    PROGRAM_FORMAT.format(d),
                                    DESCRIPTION_FORMAT.format(d),
                                    files,
                                    files,
                                    SAMPLE_FORMAT.format(d))


@pytest.fixture(scope="function")
def tmp_test_dir(tsv_file):
    """
    Create a temporary directory with subdirectories and files.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :return: path to directory
    :rtype: str or unicode
    """
    tmp_test_dir = tempfile.mkdtemp(
        prefix="tmp_riboviz_test_workflow_record")
    for d in range(0, NUM_DIRECTORIES):
        sub_dir = os.path.join(tmp_test_dir,
                               DIRECTORY_FORMAT.format(d))
        os.mkdir(sub_dir)
        for f in range(0, NUM_FILES):
            file_name = os.path.join(sub_dir,
                                     FILE_FORMAT.format(d, f))
            open(file_name, 'w').close()  # Create empty file
    record_test_files(tmp_test_dir, tsv_file)
    yield tmp_test_dir
    shutil.rmtree(tmp_test_dir)


def test_create_record(tsv_file):
    """
    Test create_record_file and check file contents have a valid
    header only.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    """
    workflow_record.create_record_file(tsv_file, delimiter="\t")
    records_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    records = records_df.to_dict('records')
    assert len(records) == 0, "Expected 0 records"
    assert list(records_df.columns).sort() == workflow_record.HEADER.sort(), \
        "Headers do not match"


@pytest.mark.parametrize("read_or_write",
                         [workflow_record.READ, workflow_record.WRITE])
def test_get_record_row(read_or_write):
    """
    Test get_record_row returns a dict with the expected contents.

    :param read_or_write: READ or WRITE
    :type read_or_write: str or unicode
    """
    row = workflow_record.get_record_row(
        "some sample",
        "some description",
        "some program",
        "some file",
        read_or_write)
    assert row[workflow_record.SAMPLE_NAME] == "some sample", \
        "Unexpected SAMPLE_NAME"
    assert row[workflow_record.DESCRIPTION] == "some description", \
        "Unexpected DESCRIPTION"
    assert row[workflow_record.PROGRAM] == "some program", \
        "Unexpected PROGRAM"
    assert row[workflow_record.FILE] == "some file", \
        "Unexpected FILE"
    assert row[workflow_record.READ_WRITE] == read_or_write, \
        "Unexpected READ_WRITE"


def test_get_record_row_invalid_read_or_write():
    """
    Test get_record_row raises an error if given an invalid value for
    read_or_Write.
    """
    with pytest.raises(AssertionError):
        workflow_record.get_record_row(
            "some sample",
            "some description",
            "some program",
            "some file",
            "invalid")


@pytest.mark.parametrize("num_reads", [0, NUM_FILES])
@pytest.mark.parametrize("num_writes", [0, NUM_FILES])
@pytest.mark.parametrize("include_sample", [True, False])
def test_record_step(tsv_file, num_reads, num_writes, include_sample):
    """
    Test create_record then record_step and check file contents are
    as expected.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param num_reads: Number of read files to record
    :type num_reads: int
    :param num_writes: Number of written files to record
    :type num_writes: int
    :param include_sample: Record a sample ID?
    :type include_sample: bool
    """
    workflow_record.create_record_file(tsv_file, delimiter="\t")
    tag = "test"
    read_files = [FILE_FORMAT.format(tag, f)
                  for f in range(0, num_reads)]
    write_files = [FILE_FORMAT.format(tag, f)
                   for f in range(0, num_writes)]
    if include_sample:
        workflow_record.record_step(tsv_file,
                                    PROGRAM_FORMAT.format(tag),
                                    DESCRIPTION_FORMAT.format(tag),
                                    read_files,
                                    write_files,
                                    SAMPLE_FORMAT.format(tag))
    else:
        workflow_record.record_step(tsv_file,
                                    PROGRAM_FORMAT.format(tag),
                                    DESCRIPTION_FORMAT.format(tag),
                                    read_files,
                                    write_files)
    records_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    # Pandas treats missing cells as nan, so convert to ""
    records_df.fillna(value={workflow_record.SAMPLE_NAME: ""}, inplace=True)
    records = records_df.to_dict('records')

    assert len(records) == num_reads + num_writes, \
        "Expected {} records".format(len(records))
    for record in records:
        assert record[workflow_record.PROGRAM] == \
            PROGRAM_FORMAT.format(tag)
        assert record[workflow_record.DESCRIPTION] == \
            DESCRIPTION_FORMAT.format(tag)
        if include_sample:
            assert record[workflow_record.SAMPLE_NAME] == \
                SAMPLE_FORMAT.format(tag)
        else:
            assert record[workflow_record.SAMPLE_NAME] == ""
    read_records = records[0:num_reads]
    file_records = list(zip(range(num_reads), read_records))
    for n, record in file_records:
        assert record[workflow_record.FILE] == \
            FILE_FORMAT.format(tag, n)
        assert record[workflow_record.READ_WRITE] == workflow_record.READ
    write_records = records[num_reads:num_writes]
    file_records = list(zip(range(num_writes), write_records))
    for n, record in file_records:
        assert record[workflow_record.FILE] == \
            FILE_FORMAT.format(tag, n)
        assert record[workflow_record.READ_WRITE] == workflow_record.WRITE


def test_validate_records(tsv_file, tmp_test_dir):
    """
    Test validate_record.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param tmp_test_dir: Temporary directory with subdirectories and
    files
    :type tmp_test_dir: str or unicode
    """
    workflow_record.validate_records(tsv_file, [tmp_test_dir])


def test_validate_records_missing_record(tsv_file, tmp_test_dir):
    """
    Test validate_record where the directory contains a non-recorded
    file and check that an AssertionError is raised.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param tmp_test_dir: Temporary directory with subdirectories and
    files
    :type tmp_test_dir: str or unicode
    """
    sub_dir = os.path.join(tmp_test_dir,
                           DIRECTORY_FORMAT.format(1))
    file_name = os.path.join(sub_dir,
                             FILE_FORMAT.format(1, "unrecorded"))
    open(file_name, 'w').close()
    with pytest.raises(AssertionError):
        workflow_record.validate_records(tsv_file, [tmp_test_dir])


def test_validate_records_missing_file(tsv_file, tmp_test_dir):
    """
    Test validate_record where the workflow record contains a
    record for a file that does not exist check that an AssertionError
    is raised.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param tmp_test_dir: Temporary directory with subdirectories and
    files
    :type tmp_test_dir: str or unicode
    """
    workflow_record.record_step(tsv_file,
                                PROGRAM_FORMAT.format(123),
                                DESCRIPTION_FORMAT.format(123),
                                ["missing_file.txt"],
                                [])
    with pytest.raises(AssertionError):
        workflow_record.validate_records(tsv_file, [tmp_test_dir])
