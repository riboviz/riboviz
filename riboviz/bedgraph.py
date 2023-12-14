"""
Bedgraph-related constants and functions.
"""
import pandas as pd

BEDGRAPH_EXT = "bedgraph"
""" File extension. """
BEDGRAPH_FORMAT = "{}." + BEDGRAPH_EXT
""" File name format. """
TRACK_PREFIX = "track type=bedGraph"
""" Track line prefix. """
COLUMNS = ["Chromosome", "Start", "End", "Data"]
""" Column names. """


def load_bedgraph(bed_file):
    """
    Load a bedGraph file. A bedGraph file has format:

    * Track definition line: ``track type=bedGraph ...``
    * Track data rows, each in 4 column BED format:
        - Chromosome (string)
        - Chromosome start coordinate (integer)
        - Chromosome end coordinate (integer)
        - Data value (integer, as riboviz uses bedGraphs for counts)

    The ``DataFrame`` returned has four column names: ``Chromosome``,
    ``Start``, ``End``, ``Data``.

    :param bed_file: File name
    :type bed_file: str or unicode
    :return: bedGraph track definition line and data
    :rtype: tuple(str or unicode, pandas.core.frame.DataFrame)
    :raise AssertionError: If the first line of the file does \
    not start with ``track type=bedGraph`` or the rest of the file \
    does not have 4 columns
    :raise Exception: if any problems arise
    """
    with open(bed_file) as f:
        track = f.readline()
    assert track.startswith(TRACK_PREFIX),\
        "Invalid bedgraph file: %s. Invalid track line: %s"\
        % (bed_file, track)
    data = pd.read_csv(bed_file, sep="\t", header=None, skiprows=1)
    assert data.shape[1] == len(COLUMNS),\
        "Invalid bedgraph file: %s. Expected 4 columns, found %d"\
        % (bed_file, data.shape[1])
    data.columns = COLUMNS
    return (track, data)


def equal_bedgraph(file1, file2):
    """
    Compare two bedGraph files for equality.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: If the files are invalid bedGraph files \
    or their contents differ
    :raise Exception: If problems arise when loading the files
    """
    (track1, data1) = load_bedgraph(file1)
    (track2, data2) = load_bedgraph(file2)
    assert track1 == track2,\
        "Unequal bedGraph tracks: %s (%s), %s (%s)"\
        % (file1, track1, file2, track2)
    assert data1.shape[0] == data2.shape[0],\
        "Unequal bedGraph rows: %s (%d), %s (%d)"\
        % (file1, data1.shape[0], file2, data2.shape[0])
    assert data1.equals(data2),\
        "Unequal bedGraph data: %s, %s" % (file1, file2)
