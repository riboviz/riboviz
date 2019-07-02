"""
Helper methods for comparing files of different types.
"""
import os
import os.path
import subprocess
import sys
import pandas


def equal_sizes(file1, file2):
    """
    Compare sizes of two files.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if file sizes differ
    """
    print("CHECK: equal_sizes")
    stat1 = os.stat(file1)
    stat2 = os.stat(file2)
    assert stat1.st_size == stat2.st_size,\
        "File sizes differ: %s, %s" % (file1, file2)


def equal_names(file1, file2):
    """
    Compare local names of two files.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if file names differ
    """
    print("CHECK: equal_names")
    local_file1 = os.path.split(file1)[1]
    local_file2 = os.path.split(file2)[1]
    assert local_file1 == local_file2,\
        "File names differ: %s, %s" % (local_file1, local_file2)


def equal_h5(file1, file2):
    """
    Compare two HDF5 files for equality, using `h5diff`.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files differ
    """
    print("CHECK: equal_h5")
    # TODO Implement in-Python comparison.
    cmd = ["h5diff", file1, file2]
    return_code = subprocess.call(cmd)
    assert return_code == 0,\
        "Non-zero return code (%d) from %s" % (
            return_code, ' '.join(map(str, cmd)))


BEDGRAPH_COLUMNS = ["Chromosome", "Start", "End", "Data"]


def load_bedgraph(file):
    """
    Load bedGraph file. A bedGraph file has format:

    * Track definition line.
      - `track type=bedGraph ...` (optional)
    * Track data rows, each in 4 column BED format:
      - Chromosone (string)
      - Chromosone start coordinate (integer)
      - Chromosone end coordinate (integer)
      - Data value (integer, as RiboViz uses bedGraphs for counts)

    The DataFrame returned has four column names: "Chromosome",
    "Start", "End", "Data".

    :param file: File name
    :type file: str or unicode
    :return: bedGraph data
    :rtype: pandas.core.frame.DataFrame
    :raise AssertionError: if file does not have 4 columns
    """
    # TODO handle track-line, if provided.
    data = pandas.read_csv(file, sep="\t", header=None)
    assert data.shape[1] == len(BEDGRAPH_COLUMNS),\
        "Non-bedgraph file: %s. Expected 4 columns, found %d" % (
            str(file), data.shape[1])
    data.columns = BEDGRAPH_COLUMNS
    return data


def equal_bedgraph(file1, file2):
    """
    Compare two bedGraph files for equality.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files are not bedGraphs or they
    differ in their data
    """
    print("CHECK: equal_bedgraph")
    data1 = load_bedgraph(file1)
    data2 = load_bedgraph(file2)
    assert data1.equals(data2),\
        "BedGraph file data differs: %s, %s" % (file1, file2)


def compare(file1, file2):
    """
    Compare two bedgraph files for equality.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files differ
    """
    assert os.path.exists(file1)
    assert os.path.exists(file2)
    assert os.path.isfile(file1)
    assert os.path.isfile(file2)
#    equal_names(file1, file2)
    ext = os.path.splitext(file1)[1]
    if ext in [".pdf", ".ht2"]:
        equal_sizes(file1, file2)
    if ext in [".h5"]:
        equal_h5(file1, file2)
    if ext in [".bedgraph"]:
        equal_bedgraph(file1, file2)


if __name__ == "__main__":
    compare(sys.argv[1], sys.argv[2])
