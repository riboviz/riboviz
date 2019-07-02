"""
Helper methods for comparing files of different types.
"""
import os
import os.path
import subprocess
import sys
import pandas
from pysam import AlignmentFile


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
        "Unequal file sizes: %s, %s" % (file1, file2)


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
        "Unequal file names: %s, %s" % (local_file1, local_file2)


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
        "Invalid bedgraph file: %s. Expected 4 columns, found %d" % (
            file, data.shape[1])
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
        "Unequal bedGraph data: %s, %s" % (file1, file2)


def equal_bam(file1, file2, has_index=True):
    """
    Compare two BAM files for equality.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :param has_index: are files expected to have complementary BAI
    files?
    :type has_index: bool
    :raise AssertionError: if files differ in their data
    :raise Exception: if problems arise when loading the files or, if
    applicable, their BAI files
    """
    print(("CHECK: equal_bam(has_index=" + str(has_index) + ")"))
    with AlignmentFile(file1, mode="rb") as bam_file1,\
            AlignmentFile(file2, mode="rb") as bam_file2:
        assert bam_file1.is_bam, "Non-BAM file: %s" % file1
        assert bam_file2.is_bam, "Non-BAM file: %s" % file2
        if has_index:
            # If has_index() is True then an index file was found
            # and opened.
            assert bam_file1.has_index(), "No BAM index: %s" % file1
            assert bam_file2.has_index(), "No BAM index: %s" % file2
        stats1 = bam_file1.get_index_statistics()
        stats2 = bam_file1.get_index_statistics()
        assert stats1 == stats2, "Unequal index statistics: %s, %s"\
            % (file1, file2)
        assert bam_file1.nreferences == bam_file2.nreferences,\
            "Unequal number of reference sequences: %s (%d), %s (%d)"\
            % (file1, bam_file1.nreferences, file2, bam_file2.nreferences)
        assert bam_file1.references == bam_file2.references,\
            "Unequal reference sequence names: %s, %s"\
            % (file1, file2)
        assert bam_file1.lengths == bam_file2.lengths,\
            "Unequal reference sequence lengths: %s, %s"\
            % (file1, file2)
        assert bam_file1.nocoordinate == bam_file2.nocoordinate,\
            "Unequal number of reads without coordinates: %s (%d), %s (%d)"\
            % (file1, bam_file1.nocoordinate, file2, bam_file2.nocoordinate)
        assert bam_file1.mapped == bam_file2.mapped,\
            "Unequal number of mapped alignments: %s (%d), %s (%d)"\
            % (file1, bam_file1.mapped, file2, bam_file2.mapped)
        assert bam_file1.unmapped == bam_file2.unmapped,\
            "Unequal number of unmapped alignments: %s (%d), %s (%d)"\
            % (file1, bam_file1.unmapped, file2, bam_file2.unmapped)

        keys1 = list(bam_file1.header.keys())
        keys2 = list(bam_file2.header.keys())
        assert keys1 == keys2,\
            "Unequal header keys: %s, %s" % (file1, file2)
        for key in keys1:
            if key == "PG":
                # Skip PG as this holds command-line invocation that
                # created the file and paths within will differ.
                continue
            assert bam_file1.header[key] == bam_file2.header[key],\
                "Unequal values for key %s: %s (%s), %s (%s)"\
                % (key, file1, str(bam_file1.header[key]),
                   file1, str(bam_file2.header[key]))
        # In RiboViz order matters.
        reads1 = [read for read in bam_file1.fetch()]
        reads2 = [read for read in bam_file2.fetch()]
        i = 0
        for seg1, seg2 in zip(reads1, reads2):
            i = i + 1
            # compare returns -1 if seg1 < seg2, 0 if =, 1 if >
            comparison = seg1.compare(seg2)
            if comparison != 0:
                print("Unequal segments:")
                print(("Pair: " + str(i) + " Compare:" + str(comparison)))
                print((str(seg1) + str(seg2)))
# TODO uncomment when resolved how best to compare these.
#            assert comparison == 0,\
#                "Unequal segments: %s (%s), %s (%s)"\
#                % (file1, str(seg1), file2, str(seg2))


def compare(file1, file2):
    """
    Compare two files for equality.

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
    if ext in [".bam"]:
        equal_bam(file1, file2)


if __name__ == "__main__":
    compare(sys.argv[1], sys.argv[2])
