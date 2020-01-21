"""
riboviz.workflow_record test suite.
"""
import os
import tempfile
import pytest
import pandas as pd
from riboviz import workflow_record


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


def test_record_step(tsv_file):
    """
    Test create_record then record_step and check file contents are
    as expected.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    """
    workflow_record.create_record_file(tsv_file, delimiter="\t")
    workflow_record.record_step(
        tsv_file,
        "program0",
        "description0",
        ["rfile0a", "rfile0b"],
        ["wfile0a", "wfile0b"],
        "sample0")
    workflow_record.record_step(
        tsv_file,
        "program1",
        "description1",
        [],
        ["wfile1a", "wfile1b"],
        "sample1")
    workflow_record.record_step(
        tsv_file,
        "program2",
        "description2",
        ["rfile2a", "rfile2b"],
        ["wfile2a", "wfile2b"])
    workflow_record.record_step(
        tsv_file,
        "program3",
        "description3",
        ["rfile3a", "rfile3b"],
        [],
        "sample3")
    records_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    # Pandas treats missing cells as nan, so convert to ""
    records_df.fillna(value={workflow_record.SAMPLE_NAME: ""}, inplace=True)
    records = records_df.to_dict('records')
    assert len(records) == 12, "Expected 12 records"

    for record in records[0:4]:
        assert record[workflow_record.PROGRAM] == "program0"
        assert record[workflow_record.DESCRIPTION] == "description0"
        assert record[workflow_record.SAMPLE_NAME] == "sample0"
    for record in records[0:2]:
        assert record[workflow_record.READ_WRITE] == workflow_record.READ
    for record in records[2:4]:
        assert record[workflow_record.READ_WRITE] == workflow_record.WRITE
    assert records[0][workflow_record.FILE] == "rfile0a"
    assert records[1][workflow_record.FILE] == "rfile0b"
    assert records[2][workflow_record.FILE] == "wfile0a"
    assert records[3][workflow_record.FILE] == "wfile0b"

    for record in records[4:6]:
        assert record[workflow_record.PROGRAM] == "program1"
        assert record[workflow_record.DESCRIPTION] == "description1"
        assert record[workflow_record.SAMPLE_NAME] == "sample1"
    for record in records[4:6]:
        assert record[workflow_record.READ_WRITE] == workflow_record.WRITE
    assert records[4][workflow_record.FILE] == "wfile1a"
    assert records[5][workflow_record.FILE] == "wfile1b"

    for record in records[6:10]:
        assert record[workflow_record.PROGRAM] == "program2"
        assert record[workflow_record.DESCRIPTION] == "description2"
        assert record[workflow_record.SAMPLE_NAME] == ""
    for record in records[6:8]:
        assert record[workflow_record.READ_WRITE] == workflow_record.READ
    for record in records[8:10]:
        assert record[workflow_record.READ_WRITE] == workflow_record.WRITE
    assert records[6][workflow_record.FILE] == "rfile2a"
    assert records[7][workflow_record.FILE] == "rfile2b"
    assert records[8][workflow_record.FILE] == "wfile2a"
    assert records[9][workflow_record.FILE] == "wfile2b"

    for record in records[10:]:
        assert record[workflow_record.PROGRAM] == "program3"
        assert record[workflow_record.DESCRIPTION] == "description3"
        assert record[workflow_record.SAMPLE_NAME] == "sample3"
    for record in records[10:]:
        assert record[workflow_record.READ_WRITE] == workflow_record.READ
    assert records[10][workflow_record.FILE] == "rfile3a"
    assert records[11][workflow_record.FILE] == "rfile3b"


def test_record_step_no_inputs_outputs(tsv_file):
    """
    Test create_record then record_step with empty input and output
    lists and check file contents are as expected.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    """
    workflow_record.create_record_file(tsv_file, delimiter="\t")
    workflow_record.record_step(
        tsv_file,
        "some program",
        "some description",
        [],
        [],
        "some sample")
    records_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    records = records_df.to_dict('records')
    assert len(records) == 0, "Expected 0 records"
