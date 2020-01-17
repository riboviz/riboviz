"""
Helper methods for comparing and validating files of different types.
"""
import os
import os.path
import subprocess
from Bio import SeqIO
import numpy as np
import pandas as pd
import pysam
from riboviz import bedgraph
from riboviz import sam_bam


def equal_names(file1, file2):
    """
    Compare local names of two files each of which must exist
    and be a file.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if file names differ
    :raise Exception: if problems arise when loading the files
    """
    local_file1 = os.path.split(file1)[1]
    local_file2 = os.path.split(file2)[1]
    assert os.path.exists(file1) and os.path.isfile(file1),\
        "File %s does not exist or is not a file"
    assert os.path.exists(file2) and os.path.isfile(file2),\
        "File %s does not exist or is not a file"
    assert local_file1 == local_file2,\
        "Unequal file names: %s, %s" % (local_file1, local_file2)


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
    stat1 = os.stat(file1)
    stat2 = os.stat(file2)
    assert stat1.st_size == stat2.st_size,\
        "Unequal file sizes: %s, %s" % (file1, file2)


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
    # TODO implement in-Python comparison.
    cmd = ["h5diff", "-q", file1, file2]
    return_code = subprocess.call(cmd)
    assert return_code == 0,\
        "Non-zero return code (%d) from %s" % (
            return_code, ' '.join(map(str, cmd)))


def load_bedgraph(bed_file):
    """
    Load bedGraph file. A bedGraph file has format:

    * Track definition line.
      - `track type=bedGraph ...`
    * Track data rows, each in 4 column BED format:
      - Chromosone (string)
      - Chromosone start coordinate (integer)
      - Chromosone end coordinate (integer)
      - Data value (integer, as RiboViz uses bedGraphs for counts)

    The DataFrame returned has four column names: "Chromosome",
    "Start", "End", "Data".

    :param bed_file: File name
    :type bed_file: str or unicode
    :return: bedGraph track and data
    :rtype: tuple(str or unicode, pandas.core.frame.DataFrame)
    :raise AssertionError: if the first line of the file does
    not start with "track type=bedGraph" or the rest of the file
    does not have 4 columns
    :raise Exception: if problems arise when loading the file
    """
    with open(bed_file) as f:
        track = f.readline()
    assert track.startswith(bedgraph.TRACK_PREFIX),\
        "Invalid bedgraph file: %s. Invalid track line: %s"\
        % (bed_file, track)
    data = pd.read_csv(bed_file, sep="\t", header=None, skiprows=1)
    assert data.shape[1] == len(bedgraph.COLUMNS),\
        "Invalid bedgraph file: %s. Expected 4 columns, found %d"\
        % (bed_file, data.shape[1])
    data.columns = bedgraph.COLUMNS
    return (track, data)


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

    BAM files are assumed to have been sorted by their leftmost
    coordinate position.

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
    with pysam.AlignmentFile(file1, mode="rb") as bam_file1,\
            pysam.AlignmentFile(file2, mode="rb") as bam_file2:
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
        equal_bam_sam_metadata(bam_file1, bam_file2)
        equal_bam_sam_headers(bam_file1, bam_file2)
        equal_bam_sam_references(bam_file1, bam_file2)
        equal_bam_sam_reads(bam_file1, bam_file2)


def equal_sam(file1, file2):
    """
    Compare two SAM files for equality. The following are compared:

    * Header values for all but "PG".
    * Reference numbers, names and lengths.
    * Reads.

    SAM files are assumed to have been sorted by their leftmost
    coordinate position.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: if files differ in their data
    :raise Exception: if problems arise when loading the files
    """
    with pysam.AlignmentFile(file1) as sam_file1,\
            pysam.AlignmentFile(file2) as sam_file2:
        assert sam_file1.is_sam, "Non-SAM file: %s" % file1
        assert sam_file2.is_sam, "Non-SAM file: %s" % file2
        equal_bam_sam_metadata(sam_file1, sam_file2)
        equal_bam_sam_headers(sam_file1, sam_file2)
        equal_bam_sam_references(sam_file1, sam_file2)
        equal_bam_sam_reads(sam_file1, sam_file2)


def equal_bam_sam_metadata(file1, file2):
    """
    Compare BAM or SAM file metadata for equality.
    Category, version, compression, description are compared.

    :param file1: File name
    :type file1: pysam.AlignmentFile
    :param file2: File name
    :type file2: pysam.AlignmentFile
    :raise AssertionError: if files differ in their metadata
    """
    assert file1.category == file2.category,\
        "Unequal category: %s (%s), %s (%s)"\
        % (file1.filename, file1.category,
           file2.filename, file2.category)
    assert file1.version == file2.version,\
        "Unequal version: %s (%s), %s (%s)"\
        % (file1.filename, str(file1.version),
           file2.filename, str(file2.version))
    assert file1.compression == file2.compression,\
        "Unequal compression: %s (%s), %s (%s)"\
        % (file1.filename, file1.compression,
           file2.filename, file2.compression)
    assert file1.description == file2.description,\
        "Unequal description: %s (%s), %s (%s)"\
        % (file1.filename, file1.description,
           file2.filename, file2.description)


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
        if key == sam_bam.PG_TAG:
            continue
        assert file1.header[key] == file2.header[key],\
            "Unequal values for key %s: %s (%s), %s (%s)"\
            % (key, file1.filename, str(file1.header[key]),
               file2.filename, str(file2.header[key]))


def equal_bam_sam_reads(file1, file2):
    """
    Compare BAM or SAM reads for equality.

    BAM/SAM files are assumed to have been sorted by their leftmost
    coordinate position.

    :param file1: File name
    :type file1: pysam.AlignmentFile
    :param file2: File name
    :type file2: pysam.AlignmentFile
    :raise AssertionError: if files differ in their reads
    """
    # Get total number of reads in each file.
    with pysam.AlignmentFile(file1.filename) as f1:
        num_reads1 = f1.count()
    with pysam.AlignmentFile(file2.filename) as f2:
        num_reads2 = f2.count()
    assert num_reads1 == num_reads2,\
        "Unequal read counts: %s (%d), %s (%d)"\
        % (file1.filename, num_reads1, file2.filename, num_reads2)
    if num_reads1 == 0:
        return
    # Iterate through files in batches to reduce memory overheads.
    # Exit when all reads have been processed.
    iter1 = file1.fetch(until_eof=True)
    iter2 = file2.fetch(until_eof=True)
    read1 = next(iter1)  # Guaranteed to be at least one.
    read_pos = read1.pos
    num_reads_done1 = 1
    num_reads_done2 = 0
    all_done = False
    while not all_done:
        reads1 = [read1]
        current_pos = read_pos
        next_pos_read = None
        # 1. Iterate through file1 to get all references with the same
        #    leftmost coordinate position.
        while read_pos == current_pos:
            try:
                read1 = next(iter1)
                reads1.append(read1)
                read_pos = read1.pos
                num_reads_done1 = num_reads_done1 + 1
            except StopIteration:
                all_done = True
                break
        positions1 = {r1.pos for r1 in reads1}
        # Internal check to ensure that logic above:
        # * Either returns a batch with two positions i.e. a batch
        #   where all reads share a common position, plus one extra
        #   read representing the start of the next
        #   batch. Unfortunately we can't look ahead in the iterator
        #   so need to handle the extra read from the next batch.
        # * Or returns a batch with one position i.e. a batch that
        #   is terminated by the end of the file.
        assert len(positions1) <= 2,\
            "More than two positions in batch: %s (%s)"\
            % (file1.filename, str(positions1))
        # 3. If we've read the first read of the next batch then
        #    remove its read from reads1 and keep for the next
        #    iteration.
        if len(positions1) == 2:
            next_pos_read = reads1.pop()
        # 2. Read the same number of reads from file2. This
        #    should not fail as it is already known that both files
        #    have the same number of reads.
        reads2 = []
        for _ in range(len(reads1)):
            read2 = next(iter2)
            reads2.append(read2)
            num_reads_done2 = num_reads_done2 + 1
        positions2 = {r2.pos for r2 in reads2}
        # 3. Check that the reads from file2 all share a common
        #    position and that this position is the same as for
        #    file1. If not, then the files differ.
        assert len(positions2) == 1,\
            "Unequal positions: %s (%s), %s (%s)"\
            % (file1.filename, current_pos,
               file2.filename, str(positions2))
        position2 = positions2.pop()  # Get only member.
        assert current_pos == position2,\
            "Unequal positions: %s (%s), %s (%s)"\
            % (file1.filename, str(current_pos),
               file2.filename, str(position2))
        # 4. Check that the reads from each file are the same.
        reads1.sort(key=get_segment_qname)
        reads2.sort(key=get_segment_qname)
        assert reads1 == reads2,\
            "Unequal reads at position %s: %s, %s"\
            % (str(current_pos), file1.filename, file2.filename)
        # 5. Initialise for next iteration.
        if next_pos_read is not None:
            read1 = next_pos_read
        else:
            all_done = True
    # Internal checks to ensure that logic above does check all reads.
    assert num_reads_done1 == num_reads1,\
        "Reads processed (%d) not equal to reads (%d) in file (%s)"\
        % (num_reads_done1, num_reads1, file1.filename)
    assert num_reads_done2 == num_reads2,\
        "Reads processed (%d) not equal to reads (%d) in file (%s)"\
        % (num_reads_done2, num_reads2, file2.filename)


def get_segment_qname(segment):
    """
    Return the qualified name of a read segment

    :param segment: read segment
    :type segment: pysam.libcalignedsegment.AlignedSegment
    :return: qualified name
    :rtype: str or unicode
    """
    return segment.qname


def equal_dataframes(data1, data2, tolerance=0.0001):
    """
    Compare two Pandas dataframes for equality.

    The dataframes are expected to be two dimensional i.e. rows and
    columns.

    Dataframes are compared column-by-column:

    * float64 columns are converted to numpy arrays then tested for
      equality to within the given tolerance using
      numpy.allclose. This is used instead of
      pandas.testing.assert_frame_equal as there is an issue with how
      that function handles precision (see
      pandas.testing.assert_frame_equal doesn't do precision according
      to the doc #25068,
      https://github.com/pandas-dev/pandas/issues/25068). In addition,
      "NAN" values are considered to be equal.
    * All other columns (object, int64, bool, datetime64, timedelta)
      are compared for exact equality using
      pandas.core.series.Series.equals.

    :param data1: dataframe
    :type data1: pandas.core.frame.DataFrame
    :param data2: dataframe
    :type data2: pandas.core.frame.DataFrame
    :param tolerance: tolerance for floating point comparisons
    :type tolerance: float
    :raise AssertionError: if DataFrames differ in their data
    """
    assert data1.shape == data2.shape,\
        "Unequal shape: %s, %s"\
        % (str(data1.shape), str(data2.shape))
    assert data1.columns.equals(data2.columns),\
        "Unequal column names: %s, %s"\
        % (str(data1.columns), str(data2.columns))
    for column in data1.columns:
        column1 = data1[column]
        column2 = data2[column]
        if column1.dtype in (int, float) and column2.dtype in (int, float):
            column_data1 = column1.to_numpy()
            column_data2 = column2.to_numpy()
            assert np.allclose(column_data1,
                               column_data2,
                               rtol=0,
                               atol=tolerance,
                               equal_nan=True),\
                "Unequal column values: %s" % column
        else:
            assert column1.equals(column2),\
                "Unequal column values: %s" % column


def equal_tsv(file1, file2, tolerance=0.0001, comment="#"):
    """
    Compare two tab-separated (TSV) files for equality.

    See equal_dataframes.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :param tolerance: tolerance for floating point comparisons
    :type tolerance: float
    :param comment: Comment prefix
    :type comment: str or unicode
    :raise AssertionError: if files differ in their data
    :raise Exception: if problems arise when loading the files
    """
    data1 = pd.read_csv(file1, sep="\t", comment=comment)
    data2 = pd.read_csv(file2, sep="\t", comment=comment)
    try:
        equal_dataframes(data1, data2, tolerance)
    except AssertionError as e:
        # Add file names to error message.
        message = e.args[0]
        message += " in file: " + str(file1) + ":" + str(file2)
        e.args = (message,)
        raise


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
    seqs1 = {}
    for seq1 in SeqIO.parse(file1, "fastq"):
        seqs1[seq1.name] = seq1
    for seq2 in SeqIO.parse(file2, "fastq"):
        assert seq2.name in seqs1,\
            "Missing ID: %s in %s but not in %s"\
            % (seq2.name, file2, file1)
        seq1 = seqs1[seq2.name]
        del seqs1[seq2.name]
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
    assert not seqs1,\
        "Missing IDs: %s in %s but not in %s"\
        % (str(seqs1), file1, file2)


def compare(file1, file2, compare_names=True):
    """
    Compare two files for equality. The following functions are used
    to compare each type of file:

    * .pdf: equal_names(file1, file2)
    * .ht2, .bai: equal_sizes(file1, file2)
    * .h5: equal_h5(file1, file2)
    * .bedgraph: equal_bedgraph(file1, file2)
    * .bam: equal_bam(file1, file2)
    * .sam: equal_sam(file1, file2)
    * .tsv: equal_tsv(file1, file2)
    * .fq: equal_fastq(file1, file2)

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :param compare_names: compare file names?
    :type: bool
    :raise AssertionError: if files differ
    """
    assert os.path.exists(file1), "Non-existent file: %s" % file1
    assert os.path.exists(file2), "Non-existent file: %s" % file2
    assert not os.path.isdir(file1), "Directory: %s" % file1
    assert not os.path.isdir(file2), "Directory: %s" % file2
    if compare_names:
        equal_names(file1, file2)
    ext = os.path.splitext(file1)[1]
    if ext in [".pdf"]:
        equal_names(file1, file2)
    if ext in [".ht2", ".bai"]:
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
