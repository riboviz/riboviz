"""
riboviz.sample_names_files test suite.
"""
import os
import tempfile
import pytest
import pandas as pd
from riboviz import sample_names_files


@pytest.fixture(scope="function")
def tsv_file():
    """
    Create a temporary TSV file to write the data to.

    :return: path to TSV file
    :rtype: str or unicode
    """
    _, tsv_file = tempfile.mkstemp(prefix="tmp", suffix=".tsv")
    yield tsv_file
    if os.path.exists(tsv_file):
        os.remove(tsv_file)


def test_create_file(tsv_file):
    """
    Test create_file and check file contents have a valid
    header only.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    """
    sample_names_files.create_file(tsv_file, delimiter="\t")
    records_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    records = records_df.to_dict('records')
    assert len(records) == 0, "Expected 0 logs"
    assert list(records_df.columns).sort() == \
        sample_names_files.HEADER.sort(), \
        "Headers do not match"


@pytest.mark.parametrize("create_file", [True, False])
def test_record_sample(tsv_file, create_file):
    """
    Test create_file then record_sample and check file contents are as
    expected.

    :param tsv_file: path to TSV file
    :type tsv_file: str or unicode
    :param create_file: Create file before recording?
    :type create_file: bool
    """
    if create_file:
        sample_names_files.create_file(tsv_file, delimiter="\t")
    prefixes = ["A", "B", "C"]
    file_format = "{}.dat"
    for prefix in prefixes:
        sample_names_files.record_sample(tsv_file,
                                         prefix,
                                         file_format.format(prefix))
    records_df = pd.read_csv(tsv_file, sep="\t", comment="#")
    records = records_df.to_dict('records')
    assert len(records) == len(prefixes), \
        "Expected {} records".format(len(records))
    for record in records:
        assert record[sample_names_files.SAMPLE_NAME] in prefixes
        prefixes.remove(record[sample_names_files.SAMPLE_NAME])
        assert record[sample_names_files.FILE] == \
            file_format.format(record[sample_names_files.SAMPLE_NAME])
