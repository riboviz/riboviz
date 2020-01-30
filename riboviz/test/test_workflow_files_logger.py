"""
riboviz.workflow_files_logger test suite.
"""
import os
import shutil
import tempfile
import pytest
import pandas as pd
from riboviz import workflow_files_logger


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
    Create a temporary TSV file to write the worfklow logs to.

    :return: path to TSV file
    :rtype: str or unicode
    """
    _, tsv_file = tempfile.mkstemp(prefix="tmp", suffix=".tsv")
    yield tsv_file
    if os.path.exists(tsv_file):
        os.remove(tsv_file)


def log_test_files(directory, log_file):
    """
    Populate workflow files log file with entries corresponding to
    directory and file structure created by tmp_test_dir.

    :param directory: Temporary directory with subdirectories and
    files
    :type directory: str or unicode
    :param log_file: path to TSV file
    :type log_file: str or unicode
    """
    workflow_files_logger.create_log_file(log_file, delimiter="\t")
    for d in range(0, NUM_DIRECTORIES):
        sub_dir = os.path.join(directory,
                               DIRECTORY_FORMAT.format(d))
        files = [os.path.join(sub_dir, FILE_FORMAT.format(d, f))
                 for f in range(0, NUM_FILES)]
        workflow_files_logger.log_files(log_file,
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
        prefix="tmp_riboviz_test_workflow_files_logger")
    for d in range(0, NUM_DIRECTORIES):
        sub_dir = os.path.join(tmp_test_dir,
                               DIRECTORY_FORMAT.format(d))
        os.mkdir(sub_dir)
        for f in range(0, NUM_FILES):
            file_name = os.path.join(sub_dir,
                                     FILE_FORMAT.format(d, f))
            open(file_name, 'w').close()  # Create empty file
    log_test_files(tmp_test_dir, tsv_file)
    yield tmp_test_dir
    shutil.rmtree(tmp_test_dir)


def test_create_log_file(tsv_file):
    """
    Test create_log_file and check file contents have a valid
    header only.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    """
    workflow_files_logger.create_log_file(tsv_file, delimiter="\t")
    logs_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    logs = logs_df.to_dict('records')
    assert len(logs) == 0, "Expected 0 logs"
    assert list(logs_df.columns).sort() == \
        workflow_files_logger.HEADER.sort(), \
        "Headers do not match"


@pytest.mark.parametrize("read_or_write",
                         [workflow_files_logger.READ,
                          workflow_files_logger.WRITE])
def test_get_log_entry(read_or_write):
    """
    Test get_log_entry returns a dict with the expected contents.

    :param read_or_write: READ or WRITE
    :type read_or_write: str or unicode
    """
    row = workflow_files_logger.get_log_entry(
        "some sample",
        "some description",
        "some program",
        "some file",
        read_or_write)
    assert row[workflow_files_logger.SAMPLE_NAME] == "some sample", \
        "Unexpected SAMPLE_NAME"
    assert row[workflow_files_logger.DESCRIPTION] == "some description", \
        "Unexpected DESCRIPTION"
    assert row[workflow_files_logger.PROGRAM] == "some program", \
        "Unexpected PROGRAM"
    assert row[workflow_files_logger.FILE] == "some file", \
        "Unexpected FILE"
    assert row[workflow_files_logger.READ_WRITE] == read_or_write, \
        "Unexpected READ_WRITE"


def test_get_log_entry_invalid_read_or_write():
    """
    Test get_log_entry raises an error if given an invalid value for
    read_or_Write.
    """
    with pytest.raises(AssertionError):
        workflow_files_logger.get_log_entry("some sample",
                                            "some description",
                                            "some program",
                                            "some file",
                                            "invalid")


@pytest.mark.parametrize("num_reads", [0, NUM_FILES])
@pytest.mark.parametrize("num_writes", [0, NUM_FILES])
@pytest.mark.parametrize("include_sample", [True, False])
def test_log_files(tsv_file, num_reads, num_writes, include_sample):
    """
    Test create_log_file then log_files and check file contents are
    as expected.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param num_reads: Number of read files to log
    :type num_reads: int
    :param num_writes: Number of written files to log
    :type num_writes: int
    :param include_sample: Log a sample ID?
    :type include_sample: bool
    """
    workflow_files_logger.create_log_file(tsv_file, delimiter="\t")
    tag = "test"
    read_files = [FILE_FORMAT.format(tag, f)
                  for f in range(0, num_reads)]
    write_files = [FILE_FORMAT.format(tag, f)
                   for f in range(0, num_writes)]
    if include_sample:
        workflow_files_logger.log_files(
            tsv_file,
            PROGRAM_FORMAT.format(tag),
            DESCRIPTION_FORMAT.format(tag),
            read_files,
            write_files,
            SAMPLE_FORMAT.format(tag))
    else:
        workflow_files_logger.log_files(
            tsv_file,
            PROGRAM_FORMAT.format(tag),
            DESCRIPTION_FORMAT.format(tag),
            read_files,
            write_files)
    logs_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    # Pandas treats missing cells as nan, so convert to ""
    logs_df.fillna(value={workflow_files_logger.SAMPLE_NAME: ""},
                   inplace=True)
    logs = logs_df.to_dict('records')

    assert len(logs) == num_reads + num_writes, \
        "Expected {} logs".format(len(logs))
    for log in logs:
        assert log[workflow_files_logger.PROGRAM] == \
            PROGRAM_FORMAT.format(tag)
        assert log[workflow_files_logger.DESCRIPTION] == \
            DESCRIPTION_FORMAT.format(tag)
        if include_sample:
            assert log[workflow_files_logger.SAMPLE_NAME] == \
                SAMPLE_FORMAT.format(tag)
        else:
            assert log[workflow_files_logger.SAMPLE_NAME] == ""
    read_logs = logs[0:num_reads]
    file_logs = list(zip(range(num_reads), read_logs))
    for n, log in file_logs:
        assert log[workflow_files_logger.FILE] == \
            FILE_FORMAT.format(tag, n)
        assert log[workflow_files_logger.READ_WRITE] == \
            workflow_files_logger.READ
    write_logs = logs[num_reads:num_writes]
    file_logs = list(zip(range(num_writes), write_logs))
    for n, log in file_logs:
        assert log[workflow_files_logger.FILE] == \
            FILE_FORMAT.format(tag, n)
        assert log[workflow_files_logger.READ_WRITE] == \
            workflow_files_logger.WRITE


@pytest.mark.parametrize("num_inputs", [0, NUM_FILES])
@pytest.mark.parametrize("include_sample", [True, False])
def test_log_input_files(tsv_file, num_inputs, include_sample):
    """
    Test create_log_file then log_input_files and check file contents
    are as expected.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param num_inputs: Number of read files to log
    :type num_inputs: int
    :param include_sample: Log a sample ID?
    :type include_sample: bool
    """
    workflow_files_logger.create_log_file(tsv_file, delimiter="\t")
    tag = "test"
    input_files = [FILE_FORMAT.format(tag, f)
                   for f in range(0, num_inputs)]
    if include_sample:
        workflow_files_logger.log_input_files(
            tsv_file,
            input_files,
            SAMPLE_FORMAT.format(tag))
    else:
        workflow_files_logger.log_input_files(
            tsv_file,
            input_files)
    logs_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    # Pandas treats missing cells as nan, so convert to ""
    logs_df.fillna(value={workflow_files_logger.SAMPLE_NAME: ""},
                   inplace=True)
    logs_df.fillna(value={workflow_files_logger.DESCRIPTION: ""},
                   inplace=True)
    logs = logs_df.to_dict('records')

    assert len(logs) == num_inputs, \
        "Expected {} logs".format(len(logs))
    for log in logs:
        assert log[workflow_files_logger.PROGRAM] == \
            workflow_files_logger.INPUT
        assert log[workflow_files_logger.DESCRIPTION] == ""
        if include_sample:
            assert log[workflow_files_logger.SAMPLE_NAME] == \
                SAMPLE_FORMAT.format(tag)
        else:
            assert log[workflow_files_logger.SAMPLE_NAME] == ""
    input_logs = logs[0:num_inputs]
    file_logs = list(zip(range(num_inputs), input_logs))
    for n, log in file_logs:
        assert log[workflow_files_logger.FILE] == \
            FILE_FORMAT.format(tag, n)
        assert log[workflow_files_logger.READ_WRITE] == \
            workflow_files_logger.READ


def test_validate_log_file(tsv_file, tmp_test_dir):
    """
    Test validate_log_file.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param tmp_test_dir: Temporary directory with subdirectories and
    files
    :type tmp_test_dir: str or unicode
    """
    workflow_files_logger.validate_log_file(tsv_file, [tmp_test_dir])


def test_validate_log_file_unlogged_file(tsv_file, tmp_test_dir):
    """
    Test validate_log where the directory contains an unogged
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
                             FILE_FORMAT.format(1, "unlogged"))
    open(file_name, 'w').close()
    with pytest.raises(AssertionError):
        workflow_files_logger.validate_log_file(tsv_file, [tmp_test_dir])


def test_validate_log_file_missing_file(tsv_file, tmp_test_dir):
    """
    Test validate_log_file where the workflow record contains a
    record for a file that does not exist check that an AssertionError
    is raised.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param tmp_test_dir: Temporary directory with subdirectories and
    files
    :type tmp_test_dir: str or unicode
    """
    workflow_files_logger.log_files(
        tsv_file,
        PROGRAM_FORMAT.format(123),
        DESCRIPTION_FORMAT.format(123),
        ["missing_file.txt"],
        [])
    with pytest.raises(AssertionError):
        workflow_files_logger.validate_log_file(tsv_file, [tmp_test_dir])
