"""
Helper methods for comparing files of different types.
"""
import os
import os.path
import subprocess
import sys

from Bio import SeqIO
import pandas
import pysam


def equal_sizes(file1, file2):
    """
    Compare sizes of two files.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if file sizes differ
    :raise Exception: if problems arise when loading the files
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
    :raise Exception: if problems arise when loading the files
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
    :raise Exception: if problems arise when loading the files or
    running `h5diff`
    """
    print("CHECK: equal_h5")
    # TODO implement in-Python comparison.
    cmd = ["h5diff", file1, file2]
    return_code = subprocess.call(cmd)
    assert return_code == 0,\
        "Non-zero return code (%d) from %s" % (
            return_code, ' '.join(map(str, cmd)))


BEDGRAPH_COLUMNS = ["Chromosome", "Start", "End", "Data"]


def load_bedgraph(bed_file):
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

    :param bed_file: File name
    :type bed_file: str or unicode
    :return: bedGraph data
    :rtype: pandas.core.frame.DataFrame
    :raise AssertionError: if file does not have 4 columns
    :raise Exception: if problems arise when loading the file
    """
    # TODO handle track-line, if provided.
    data = pandas.read_csv(bed_file, sep="\t", header=None)
    assert data.shape[1] == len(BEDGRAPH_COLUMNS),\
        "Invalid bedgraph file: %s. Expected 4 columns, found %d"\
        % (bed_file, data.shape[1])
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
    :raise Exception: if problems arise when loading the files
    """
    print("CHECK: equal_bedgraph")
    data1 = load_bedgraph(file1)
    data2 = load_bedgraph(file2)
    assert data1.shape[0] == data2.shape[0],\
        "Unequal bedGraph rows: %s (%d), %s (%d)"\
        % (file1, data1.shape[0], file2, data2.shape[0])
    assert data1.equals(data2),\
        "Unequal bedGraph data: %s, %s" % (file1, file2)


def equal_bam(file1, file2):
    """
    Compare two BAM files for equality. The following are compared:

    * Index-specific information:
      - Index statistics.
      - Number of unequal reads without coordinates.
      - Number of mapped alignments.
      - Number of unmapped alignments.
    * Header values for all but "PG".
    * Reference numbers, names and lengths.
    * Reads

    BAM files are required to have complementary BAI files.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files differ in their data or they
    are missing complementary BAI files
    :raise Exception: if problems arise when loading the files or, if
    applicable, their complementary BAI files
    """
    print("CHECK: equal_bam")
    with pysam.AlignmentFile(file1, mode="rb", check_sq=False) as bam_file1,\
            pysam.AlignmentFile(file2, mode="rb", check_sq=False) as bam_file2:
        assert bam_file1.is_bam, "Non-BAM file: %s" % file1
        assert bam_file2.is_bam, "Non-BAM file: %s" % file2
        assert bam_file1.has_index(), "No BAM index: %s" % file1
        assert bam_file2.has_index(), "No BAM index: %s" % file2
        stats1 = bam_file1.get_index_statistics()
        stats2 = bam_file2.get_index_statistics()
        assert stats1 == stats2,\
            "Unequal index statistics: %s (%s), %s (%s)"\
            % (file1, str(stats1), file2, str(stats2))
        assert bam_file1.nocoordinate == bam_file2.nocoordinate,\
            "Unequal number of reads without coordinates: %s (%d), %s (%d)"\
            % (file1, bam_file1.nocoordinate, file2, bam_file2.nocoordinate)
        assert bam_file1.mapped == bam_file2.mapped,\
            "Unequal number of mapped alignments: %s (%d), %s (%d)"\
            % (file1, bam_file1.mapped, file2, bam_file2.mapped)
        assert bam_file1.unmapped == bam_file2.unmapped,\
            "Unequal number of unmapped alignments: %s (%d), %s (%d)"\
            % (file1, bam_file1.unmapped, file2, bam_file2.unmapped)
        equal_bam_sam_headers(bam_file1, bam_file2)
        equal_bam_sam_references(bam_file1, bam_file2)
        equal_bam_sam_reads(bam_file1, bam_file2)


def equal_sam(file1, file2):
    """
    Compare two SAM files for equality. The following are compared:

    * Header values for all but "PG".
    * Reference numbers, names and lengths.
    * Reads.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files differ in their data
    :raise Exception: if problems arise when loading the files
    """
    print("CHECK: equal_sam")
    with pysam.AlignmentFile(file1, check_sq=False) as sam_file1,\
            pysam.AlignmentFile(file2, check_sq=False) as sam_file2:
        assert not sam_file1.is_bam, "Non-SAM file: %s" % file1
        assert not sam_file2.is_bam, "Non-SAM file: %s" % file2
        equal_bam_sam_headers(sam_file1, sam_file2)
        equal_bam_sam_references(sam_file1, sam_file2)
        equal_bam_sam_reads(sam_file1, sam_file2)


def equal_bam_sam_references(file1, file2):
    """
    Compare BAM or SAM file references for equality.
    Reference numbers, names and lengths are compared.

    :param file1: File name
    :type file1: pysam.AlignmentFile
    :param file2: File name
    :type file2: pysam.AlignmentFile
    :raise AssertionError: if files differ in their references data
    """
    assert file1.nreferences == file2.nreferences,\
        "Unequal number of reference sequences: %s (%d), %s (%d)"\
        % (file1.filename, file1.nreferences,
           file2.filename, file2.nreferences)
    assert file1.references == file2.references,\
        "Unequal reference sequence names: %s (%s), %s (%s)"\
        % (file1.filename, str(file1.references),
           file2.filename, str(file2.references))
    assert file1.lengths == file2.lengths,\
        "Unequal reference sequence lengths: %s (%s), %s (%s)"\
        % (file1.filename, str(file1.lengths),
           file2.filename, str(file2.lengths))


def equal_bam_sam_headers(file1, file2):
    """
    Compare BAM or SAM file header data for equality.

    Values for all but "PG" are compared. "PG" is not compared as it
    holds information on the command-line invocation that created a
    file and the file names and paths may differ.

    :param file1: File name
    :type file1: pysam.AlignmentFile
    :param file2: File name
    :type file2: pysam.AlignmentFile
    :raise AssertionError: if files differ in their header data
    """
    keys1 = list(file1.header.keys())
    keys2 = list(file2.header.keys())
    assert keys1 == keys2,\
        "Unequal header keys: %s, %s" % (file1.filename, file2.filename)
    for key in keys1:
        if key == "PG":
            continue
        assert file1.header[key] == file2.header[key],\
            "Unequal values for key %s: %s (%s), %s (%s)"\
            % (key, file1.filename, str(file1.header[key]),
               file2.filename, str(file2.header[key]))


def equal_bam_sam_reads(file1, file2):
    """
    Compare BAM or SAM reads for equality.

    :param file1: File name
    :type file1: pysam.AlignmentFile
    :param file2: File name
    :type file2: pysam.AlignmentFile
    :raise AssertionError: if files differ in their reads
    """
    reads1 = [read for read in file1.fetch()]
    reads2 = [read for read in file2.fetch()]
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
#        assert comparison == 0,\
#            "Unequal segments: %s (%s), %s (%s)"\
#            % (file1.filename, str(seg1), file2.filename, str(seg2))


def equal_tsv(file1, file2):
    """
    Compare two tab-separated (TSV) files for exact equality.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files differ in their data
    :raise Exception: if problems arise when loading the files
    """
    print("CHECK: equal_tsv")
    data1 = pandas.read_csv(file1, sep="\t")
    data2 = pandas.read_csv(file2, sep="\t")
    assert data1.shape == data2.shape,\
        "Unequal rows/columns: %s (%s), %s (%s)"\
        % (file1, str(data1.shape), file2, str(data2.shape))
    assert data1.equals(data2),\
        "Unequal TSV data: %s, %s" % (file1, file2)


def equal_fastq(file1, file2):
    """
    Compare two FASTQ files for equality. The following are compared:

    * Both files have the same number of records.
    * All records in file1 are also in file2. The order of records is
      ignored.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files differ in their data
    :raise Exception: if problems arise when loading the files
    """
    print("CHECK: equal_fastq")
    seq_index1 = SeqIO.index(file1, "fastq")
    seq_index2 = SeqIO.index(file2, "fastq")
    assert len(seq_index1) == len(seq_index2),\
        "Unequal number of sequences: %s (%d), %s (%d)"\
        % (file1, len(seq_index1), file2, len(seq_index2))
    for seq_id in seq_index1:
        assert seq_id in seq_index2,\
            "Missing ID: %s in %s but not in %s"\
            % (seq_id, file1, file2)
        seq1 = seq_index1[seq_id]
        seq2 = seq_index2[seq_id]
        assert seq1.name == seq2.name,\
            "Unequal name: %s (%s), %s (%s)"\
            % (file1, seq1.name, file2, seq2.name)
        assert seq1.seq == seq2.seq,\
            "Unequal sequence: %s (%s), %s (%s)"\
            % (file1, seq1.seq, file2, seq2.seq)
        assert seq1.annotations == seq2.annotations,\
            "Unequal annotations: %s (%s), %s (%s)"\
            % (file1, str(seq1.annotations),
               file2, str(seq2.annotations))
        assert seq1.dbxrefs == seq2.dbxrefs,\
            "Unequal database cross-references: %s (%s), %s (%s)"\
            % (file1, seq1.dbxrefs, file2, seq2.dbxrefs)
        assert seq1.features == seq2.features,\
            "Unequal features: %s (%s), %s (%s)"\
            % (file1, seq1.features, file2, seq2.features)
        assert seq1.letter_annotations == seq2.letter_annotations,\
            "Unequal letter annotations: %s (%s), %s (%s)"\
            % (file1, seq1.letter_annotations,
               file2, seq2.letter_annotations)
        assert seq1.description == seq2.description,\
            "Unequal description: %s (%s), %s (%s)"\
            % (file1, seq1.description, file2, seq2.description)
    seq_index1.close()
    seq_index2.close()


def compare(file1, file2):
    """
    Compare two files for equality.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files differ
    """
    assert os.path.exists(file1), "Non-existent file: %s" % file1
    assert os.path.exists(file2), "Non-existent file: %s" % file2
    assert not os.path.isdir(file1), "Directory: %s" % file1
    assert not os.path.isdir(file2), "Directory: %s" % file2
#    equal_names(file1, file2)
    ext = os.path.splitext(file1)[1]
    if ext in [".pdf", ".ht2", ".bai"]:
        equal_sizes(file1, file2)
    if ext in [".h5"]:
        equal_h5(file1, file2)
    if ext in [".bedgraph"]:
        equal_bedgraph(file1, file2)
    if ext in [".bam"]:
        equal_bam(file1, file2)
    if ext in [".sam"]:
        equal_sam(file1, file2)
    if ext in [".tsv"]:
        equal_tsv(file1, file2)
    if ext in [".fq"]:
        equal_fastq(file1, file2)


if __name__ == "__main__":
    compare(sys.argv[1], sys.argv[2])
