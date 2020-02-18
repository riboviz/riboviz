"""
Sample sheet-related constants and functions.
"""
import errno
import os
import pandas as pd
from riboviz import provenance


SAMPLE_ID = "SampleID"
""" Column name. """
TAG_READ = "TagRead"
""" Column name. """
NUM_READS = "NumReads"
"""
Column name, used for sample sheets with information on
demultiplexed samples.
"""
TOTAL_READS = "Total"
""" ``SampleID`` value for total number of reads row. """
UNASSIGNED_TAG = "Unassigned"
""" ``SampleID`` value for number of unassigned reads row. """
UNASSIGNED_READ = "NNNNNNNNN"
""" ``TagRead`` value for number of unassigned reads row. """


def load_sample_sheet(file_name, delimiter="\t", comment="#"):
    """
    Load a sample sheet from a file. The sample sheet is assumed to
    have a header with column names ``SampleID`` and ``TagRead``.

    :param file_name: File name
    :type file_name: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    :param comment: Comment prefix
    :type comment: str or unicode
    :return: Sample sheet
    :rtype: pandas.core.frame.DataFrame
    :raise FileNotFoundError: If the file cannot be found or is \
    not a file
    :raise AssertionError: If there is no header with ``SampleID`` \
    and ``TagRead`` columns.
    """
    if not os.path.exists(file_name) or not os.path.isfile(file_name):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                file_name)
    sample_sheet = pd.read_csv(file_name,
                               comment=comment,
                               delimiter=delimiter)
    for column in [SAMPLE_ID, TAG_READ]:
        assert column in sample_sheet.columns,\
            "Missing column {} in {}".format(column, file_name)
    return sample_sheet


def load_deplexed_sample_sheet(file_name, delimiter="\t", comment="#"):
    """
    Load a sample sheet, with information about demultiplexed samples,
    from a file. The sample sheet is assumed to have a header with
    column names ``SampleID``, ``TagRead`` and ``NumReads``. See also
    :py:func:`load_sample_sheet`.

    :param file_name: File name
    :type file_name: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    :param comment: Comment prefix
    :type comment: str or unicode
    :return: Sample sheet
    :rtype: pandas.core.frame.DataFrame
    :raise FileNotFoundError: If the file cannot be found or is \
    not a file
    :raise AssertionError: If there is no header with ``SampleID``, \
    ``TagRead`` and ``NumReads`` columns.
    """
    sample_sheet = load_sample_sheet(file_name, delimiter, comment)
    assert NUM_READS in sample_sheet.columns,\
        "Missing column {} in {}".format(NUM_READS, file_name)
    return sample_sheet


def save_deplexed_sample_sheet(sample_sheet,
                               num_unassigned_reads,
                               file_name,
                               delimiter="\t"):
    """
    Save a sample sheet, with information about demultiplexed samples,
    to a file. The sample sheet is assumed to have columns
    ``SampleID``, ``TagRead`` and ``NumReads``. Two rows are appended
    to the sample sheet before saving::

        Unassigned NNNNNNNNN <num_unassigned_reads>
        Total                <total_reads>

    :param sample_sheet: Sample sheet
    :type sample_sheet: pandas.core.frame.DataFrame
    :param num_unassigned_reads: Number of unassigned reads
    :type num_unassigned_reads: int
    :param file_name: File name
    :type file_name: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    deplexed_sample_sheet = sample_sheet[[
        SAMPLE_ID,
        TAG_READ,
        NUM_READS
    ]]
    unassigned_row = pd.DataFrame([[UNASSIGNED_TAG,
                                    UNASSIGNED_READ,
                                    num_unassigned_reads]],
                                  columns=deplexed_sample_sheet.columns)
    deplexed_sample_sheet = deplexed_sample_sheet.append(unassigned_row,
                                                         ignore_index=True)
    total_reads = deplexed_sample_sheet[NUM_READS].sum()
    total_row = pd.DataFrame([[TOTAL_READS, "", total_reads]],
                             columns=deplexed_sample_sheet.columns)
    deplexed_sample_sheet = deplexed_sample_sheet.append(total_row,
                                                         ignore_index=True)
    provenance.write_provenance_header(__file__, file_name)
    deplexed_sample_sheet[list(deplexed_sample_sheet.columns)].to_csv(
        file_name, mode='a', sep=delimiter, index=False)


def get_non_zero_deplexed_samples(sample_sheet):
    """
    Given a sample sheet return the names of the samples for which one
    or more reads were found.

    The sample sheet is assumed to have columns, ``SampleID``,
    ``TagRead`` and ``NumReads``. Rows whose ``SampleID`` values are
    ``Unassigned`` or ``Total`` are ignored.

    :param sample_sheet: Sample sheet
    :type sample_sheet: pandas.core.frame.DataFrame
    :return samples: Samples
    :rtype samples: list(str or unicode)
    """
    non_zero_samples = sample_sheet[
        ~sample_sheet[SAMPLE_ID].isin([UNASSIGNED_TAG, TOTAL_READS])
        & sample_sheet[NUM_READS] != 0]
    return list(non_zero_samples[SAMPLE_ID])
