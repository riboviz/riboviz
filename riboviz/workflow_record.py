#!/usr/bin/env python
"""
RiboViz utilities for recording workflow steps.
"""
import os
import numpy as np
import pandas as pd
from riboviz import provenance


SAMPLE_NAME = "SampleName"
""" Sample name column name """
DESCRIPTION = "Description"
""" Description column name """
PROGRAM = "Program"
""" Program column name """
FILE = "File"
""" File column name """
READ_WRITE = "Read/Write"
""" Read/write file column name """
READ = "read"
""" Read file value """
WRITE = "write"
""" Write file value """
HEADER = [SAMPLE_NAME, DESCRIPTION, PROGRAM, FILE, READ_WRITE]
""" Workflow record file header """
INPUT = "input"
"""
Special value for PROGRAM field to denote input files that do not
originate from any step in a workflow
"""


def create_record_file(file_name, delimiter="\t"):
    """
    Create a workflow record file with provenance comments and a
    header row, HEADER.

    :param file_name: Workflow record file
    :type file_name: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    provenance.write_provenance_header(__file__, file_name)
    data = pd.DataFrame(columns=HEADER)
    data.to_csv(file_name, mode='a', sep=delimiter, index=False)


def record_step(file_name,
                program,
                description,
                files_read,
                files_written,
                sample_name=None,
                delimiter="\t"):
    """
    Append a workflow record to the given file. If both files_read and
    files_written are [] then this is a no-op.

    :param file_name: Workflow record file
    :type file_name: str or unicode
    :param program: Program invoked during step
    :type program: str or unicode
    :param description: Description of step
    :type description: str or unicode
    :param files_read: Files read by program
    :type files_read: list(str or unicode)
    :param files_written: Files written by program
    :type files_written: list(str or unicode)
    :param sample_name: Sample name (optional)
    :type sample_name: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    data = pd.DataFrame(columns=HEADER)
    rows = []
    if len(files_read) == 0 and len(files_written) == 0:
        return
    for read in files_read:
        rows.append(get_record_row(sample_name, description, program,
                                   read, READ))
    for write in files_written:
        rows.append(get_record_row(sample_name, description, program,
                                   write, WRITE))
    data = data.append(rows)
    data.to_csv(file_name, mode='a', sep=delimiter, index=False, header=False)


def get_record_row(sample_name,
                   description,
                   program,
                   file_name,
                   read_or_write):
    """
    Get a workflow record row constructed from the given parameters.

    :param sample_name: Sample name
    :type sample_name: str or unicode
    :param description: Description of step
    :type description: str or unicode
    :param program: Program invoked during step
    :type program: str or unicode
    :param file_name: File name
    :type file_name: str or unicode
    :param read_or_write: READ or WRITE
    :type read_or_write: str or unicode
    :return: Row with SAMPLE_NAME, DESCRIPTION, PROGRAM, FILE,
    READ_WRITE keys
    :type row: dict
    :raises AssertionError: if read_or_write is not READ or WRITE
    """
    assert read_or_write in [READ, WRITE]
    row = {SAMPLE_NAME: sample_name,
           DESCRIPTION: description,
           PROGRAM: program,
           FILE: file_name,
           READ_WRITE: read_or_write}
    return row


def validate_records(workflow_record_file, dirs):
    """
    Check each workflow record file exists and check that every file
    in the given directories is recorded in the workflow record file.

    :param workflow_record_file: Workflow record file
    :type workflow_record_file: str or unicode
    :param dirs: Directories
    :type dirs: list(str or unicode)
    :raises AssertionError: if any file does not exist or any file
    in the given directories is not recorded in the workflow record
    file
    """
    workflow_records = pd.read_csv(workflow_record_file,
                                   comment="#",
                                   delimiter="\t")
    recorded_files = list(np.unique(workflow_records[FILE].to_numpy()))
    actual_files = [
        os.path.join(dir_path, file_name)
        for dir in dirs
        for (dir_path, dir_name, file_names) in
        os.walk(os.path.expanduser(dir))
        for file_name in file_names
    ]
    for file_name in recorded_files:
        assert os.path.exists(file_name), \
            "Workflow record file not found: {}".format(file_name)
    additional_files = list(set(actual_files) - set(recorded_files))
    assert len(additional_files) == 0, \
        "{} are not in the workflow record".format(additional_files)
