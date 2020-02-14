"""
Sample sheet-related constants and functions.
"""
import errno
import os
import pandas as pd
from riboviz import provenance


SAMPLE_ID = "SampleID"
""" Sample sheet SampleID column name """
TAG_READ = "TagRead"
""" Sample sheet TagRead column name """
NUM_READS = "NumReads"
""" Sample sheet NumReads column name """
TOTAL_READS = "Total"
""" Sample sheet total reads tag """
UNASSIGNED_TAG = "Unassigned"
""" Sample sheet unassigned reads tag """
UNASSIGNED_READ = "NNNNNNNNN"
""" Sample sheet unassigned read marker """


def load_sample_sheet(filename, delimiter="\t", comment="#"):
    """
    Load sample sheet from a tab-separated values file. The sample
    sheet is assumed to have columns, "SampleID" and "TagRead".

    :param filename: File name
    :type filename: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    :param comment: Comment prefix
    :type comment: str or unicode
    :return: Sample sheet
    :rtype: pandas.core.frame.DataFrame
    :raise FileNotFoundError: if filename cannot be found or is not \
    a file
    :raise AssertionError: if the sample-sheet does not contain \
    "SampleID" or "TagRead" columns
    """
    if not os.path.exists(filename) or not os.path.isfile(filename):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                filename)
    sample_sheet = pd.read_csv(filename,
                               comment=comment,
                               delimiter=delimiter)
    for column in [SAMPLE_ID, TAG_READ]:
        assert column in sample_sheet.columns,\
            "Missing column {} in {}".format(column, filename)
    return sample_sheet


def save_deplexed_sample_sheet(sample_sheet,
                               num_unassigned_reads,
                               filename,
                               delimiter="\t"):
    """
    Save sample sheet with data on demultiplexed reads as a
    tab-separated values file. The sample sheet is assumed to have
    columns, "SampleID", "TagRead" and "NumReads". Two rows are
    appended to the sample sheet before saving:

        Unassigned NNNNNNNNN <num_unassigned_reads>
        Total                <total_reads>

    :param sample_sheet: Sample sheet
    :type sample_sheet: pandas.core.frame.DataFrame
    :param num_unassigned_reads: Number of unassigned reads
    :type num_unassigned_reads: int
    :param filename: File name
    :type filename: str or unicode
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
    provenance.write_provenance_header(__file__, filename)
    deplexed_sample_sheet[list(deplexed_sample_sheet.columns)].to_csv(
        filename, mode='a', sep=delimiter, index=False)


def load_deplexed_sample_sheet(filename, delimiter="\t", comment="#"):
    """
    Load demultiplexed sample sheet from a tab-separated values
    file. The sample sheet is assumed to have columns, "SampleID" and
    "TagRead" and "NumReads".

    :param filename: File name
    :type filename: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    :param comment: Comment prefix
    :type comment: str or unicode
    :return: Sample sheet
    :rtype: pandas.core.frame.DataFrame
    :raise FileNotFoundError: if filename cannot be found or is not \
    a file
    :raise AssertionError: if the sample-sheet does not contain \
    "SampleID" or "TagRead" columns
    """
    if not os.path.exists(filename) or not os.path.isfile(filename):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                filename)
    sample_sheet = pd.read_csv(filename,
                               comment=comment,
                               delimiter=delimiter)
    for column in [SAMPLE_ID, TAG_READ, NUM_READS]:
        assert column in sample_sheet.columns,\
            "Missing column {} in {}".format(column, filename)
    return sample_sheet


def get_non_zero_deplexed_samples(sample_sheet):
    """
    Given a sample sheet with data on demultiplexed reads as a
    tab-separated values file return the names of the samples for
    which one or more reads were found.

    The sample sheet is assumed to have columns, "SampleID", "TagRead"
    and "NumReads". Rows whose "SampleID" values are "Unassigned" and
    "Total" are ignored.

    :param sample_sheet: Sample sheet
    :type sample_sheet: pandas.core.frame.DataFrame
    :return samples: Samples
    :rtype samples: list(str or unicode)
    """
    non_zero_samples = sample_sheet[
        ~sample_sheet[SAMPLE_ID].isin([UNASSIGNED_TAG, TOTAL_READS])
        & sample_sheet[NUM_READS] != 0]
    return list(non_zero_samples[SAMPLE_ID])
