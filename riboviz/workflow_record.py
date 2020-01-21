#!/usr/bin/env python
"""
RiboViz utilities for recording workflow steps.
"""
import csv
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
    with open(file_name, 'a', newline='') as csv_file:
        writer = csv.DictWriter(csv_file,
                                fieldnames=HEADER,
                                delimiter=delimiter)
        writer.writeheader()


def record_step(file_name,
                program,
                description,
                files_read,
                files_written,
                sample_name=None,
                delimiter="\t"):
    """
    Append a workflow record to the given file.

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
    with open(file_name, 'a', newline='') as csv_file:
        writer = csv.DictWriter(csv_file,
                                fieldnames=HEADER,
                                delimiter=delimiter)
        read_rows = [get_record_row(sample_name, description, program,
                                    read, READ)
                     for read in files_read]
        write_rows = [get_record_row(sample_name, description, program,
                                     write, WRITE)
                      for write in files_written]
        writer.writerows(read_rows)
        writer.writerows(write_rows)


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
    :throws AssertionError: if read_or_write is not READ or WRITE
    """
    assert read_or_write in [READ, WRITE]
    row = {SAMPLE_NAME: sample_name,
           DESCRIPTION: description,
           PROGRAM: program,
           FILE: file_name,
           READ_WRITE: read_or_write}
    return row
