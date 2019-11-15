"""
Utilities.
"""
import csv
import errno
import itertools
import os
import pandas as pd


NUCLEOTIDES = "ACGT"
""" Nucleotide letters """
FASTQ_NAME = "{}.fastq"
""" .fastq file name format string """
FASTQ_GZ_NAME = FASTQ_NAME + ".gz"
""" .fastq.gz file name format string """
FASTQ_FORMAT = "fastq"
""" Format string for use with Bio.SeqIO.write. """
BARCODE_DELIMITER = "_"
""" Default barcode delmiter in fastq header """
UMI_DELIMITER = "_"
""" Default UMI delmiter in fastq header """
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


def value_in_dict(key, dictionary, allow_false_empty=False):
    """
    Check that a value is in a dictionary and the value is not None.

    If dictionary is

        {"A":1,"B":None,"C":{},"D":[],"E":[1],"F":True,"G":False}

    then

        value_in_dict("A", dictionary): True
        value_in_dict("B", dictionary): False
        value_in_dict("C", dictionary): False
        value_in_dict("D", dictionary): False
        value_in_dict("E", dictionary): True
        value_in_dict("F", dictionary): True
        value_in_dict("G", dictionary): False

        value_in_dict("A", dictionary, True): True
        value_in_dict("B", dictionary, True): False
        value_in_dict("C", dictionary, True): True
        value_in_dict("D", dictionary, True): True
        value_in_dict("E", dictionary, True): True
        value_in_dict("F", dictionary, True): True
        value_in_dict("G", dictionary, True): True

    :param key: Key
    :type key: -
    :param dictionary: Dictionary
    :type dictionary: dict
    :param allow_false_empty: Consider False, empty string, list, dict
    to be existant
    :type allow_false_empty: bool
    :return: True or False
    :rtype: bool
    """
    is_in = key in dictionary and dictionary[key] is not None
    if not allow_false_empty:
        is_in = is_in and bool(dictionary[key])
    return is_in


def list_to_str(lst):
    """
    Convert list to space-delimited string.

    :param lst: list
    :type lst: list
    :return: list as string
    :rtype: str or unicode
    """
    return ' '.join(map(str, lst))


def hamming_distance(str1, str2):
    """
    Returns the hamming distance between two strings.

    :param str1: String
    :type str1: str or unicode
    :param str2: String
    :type str2: str or unicode
    :return: hamming distance
    :rtype: int
    """
    return sum(1 for (a, b) in zip(str1, str2) if a != b)


def barcode_matches(record,
                    barcode,
                    mismatches=0,
                    delimiter=BARCODE_DELIMITER):
    """
    Returns True if fastq record header includes barcode, up to a given
    number of mismatches. The header is assumed to be of form:

        @...<DELIMITER><BARCODE><DELIMITER>...

    If <BARCODE> differs from barcode by greater than mismatches, o
    there is no <BARCODE> or <BARCODE> has a different length than
    barcode then False is returned.

    :param record: Record
    :type record: str or unicode
    :param barcode: Barcode
    :type barcode: str or unicode
    :param mismatches: Number of mismatches
    :type mismatches: int
    :param delimiter: Barcode delimiter
    :type delimiter: str or unicode
    :returns: True or False
    :rtype: bool
    """
    chunks = record.split(delimiter)
    if len(chunks) == 1:
        return False
    candidate = chunks[1]
    if len(candidate) != len(barcode):
        return False
    return hamming_distance(candidate, barcode) <= mismatches


def generate_barcode_pairs(filename, length=1):
    """
    Generate barcode pairs and write each pair plus the Hamming
    distance between then to a file of tab-separated values.

    :param filename: Filename
    :type filename: str or unicode
    :param length: Barcode length
    :type length: int
    """
    barcodes = [''.join(i) for i in itertools.product(NUCLEOTIDES,
                                                      repeat=length)]
    with open(filename, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for (a, b) in itertools.product(barcodes, repeat=2):
            distance = hamming_distance(a, b)
            writer.writerow([a, b, distance])


def get_fastq_filename(tag, is_gz=False):
    """
    Given a tag return a fastq[.gz] file name e.g. given "tag01"
    return "tag01.fastq".

    :param tag: Tag
    :type tag: str or unicode
    :param is_gz: If True, add .fastq.gz extension, else add .fastq
    extension
    :type is_gz: bool
    :return: filename
    :rtype: str or unicode
    """
    if is_gz:
        return FASTQ_GZ_NAME.format(tag.lower())
    return FASTQ_NAME.format(tag.lower())


def get_fastq_filenames(tags, is_gz=False):
    """
    Given a list of tags return fastq[.gz] file names.

    :param tags: Tags
    :type tags: list(str or unicode)
    :param is_gz: If True, add .fastq.gz extension, else add .fastq
    extension
    :type is_gz: bool
    :return: filenames
    :rtype: list(str or unicode)
    """
    return [get_fastq_filename(tag, is_gz) for tag in tags]


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
    :raise FileNotFoundError: if filename cannot be found or is not a
    file
    :raise AssertionError: if the sample-sheet does not contain
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
    deplexed_sample_sheet[list(deplexed_sample_sheet.columns)].to_csv(
        filename, sep=delimiter, index=False)


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
    :raise FileNotFoundError: if filename cannot be found or is not a
    file
    :raise AssertionError: if the sample-sheet does not contain
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
