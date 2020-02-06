#!/usr/bin/env python
"""
RiboViz utilities for recording information about which sample names
correspond to which sample files.
"""
import os.path
import pandas as pd
from riboviz import provenance


SAMPLE_NAME = "SampleName"
""" Sample name column name """
FILE = "File"
""" File column name """
HEADER = [SAMPLE_NAME, FILE]
""" Sample files header """


def create_file(file_name, delimiter="\t"):
    """
    Create a samples names file with provenance comments and a
    header row, HEADER.

    :param file_name: Samples file
    :type file_name: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    provenance.write_provenance_header(__file__, file_name)
    data = pd.DataFrame(columns=HEADER)
    data.to_csv(file_name, mode='w', sep=delimiter, index=False)


def record_sample(file_name, sample_name, sample_file, delimiter="\t"):
    """
    Append an entry for a sample to a samples file. If file_name does
    not exist it is created.

    :param file_name: Samples file
    :type file_name: str or unicode
    :param sample_name: Sample name
    :type sample_name: str or unicode
    :param sample_file: Sample file name
    :type sample_file: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    if not os.path.exists(file_name) or \
       (os.path.isfile(file_name) and os.path.getsize(file_name) == 0):
        create_file(file_name)
    data = pd.DataFrame(columns=HEADER)
    data = data.append([{SAMPLE_NAME: sample_name, FILE: sample_file}])
    data.to_csv(file_name, mode='a', sep=delimiter, index=False,
                header=False)
