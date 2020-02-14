"""
SAM and BAM-related constants and functions.
"""
import pysam
from riboviz import utils

PG_TAG = "PG"
""" SAM file PG (program) tag """
SAM_EXT = "sam"
""" SAM file extension """
BAM_EXT = "bam"
""" BAM file extension """
BAI_EXT = "bai"
""" BAI file extension """
SAM_FORMAT = "{}." + SAM_EXT
""" Format string for SAM files """
BAM_FORMAT = "{}." + BAM_EXT
""" Format string for BAM files """
BAI_FORMAT = "{}." + BAI_EXT
"""
Format string for BAM index file based on default behaviour of
"samtools index"
"""


def is_bam(file_name):
    """
    Does the given file end with .bam or .BAM?

    :param file_name: File name
    :type file_name: str or unicode
    :return: True if file_name ends with .bam or .BAM, False otherwise
    :rtype: bool
    """
    ext = utils.get_file_ext(file_name)
    return ext.lower() == BAM_EXT


def is_sam(file_name):
    """
    Does the given file end with .sam or .SAM?

    :param file_name: File name
    :type file_name: str or unicode
    :return: True if file_name ends with .sam or .SAM, False otherwise
    :rtype: bool
    """
    ext = utils.get_file_ext(file_name)
    return ext.lower() == SAM_EXT


def count_sequences(file_name):
    """
    Count number of sequences and mapped (primary aligned) sequences
    in a SAM or BAM file.

    :param file_name: SAM/BAM file name
    :type file_name: str or unicode
    :return: (number of sequences, number of mapped sequences)
    :rtype: tuple(int, int)
    """
    if is_bam(file_name):
        mode = "rb"
    else:
        mode = "r"
    num_sequences = 0
    num_mapped_sequences = 0
    with pysam.AlignmentFile(file_name, mode=mode) as f:
        for sequence in f:
            num_sequences = num_sequences + 1
            if sequence.is_unmapped or sequence.is_secondary:
                continue
            num_mapped_sequences = num_mapped_sequences + 1
    return (num_sequences, num_mapped_sequences)


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
    :raise AssertionError: if files differ in their data or they \
    are missing complementary BAI files
    :raise Exception: if problems arise when loading the files or, \
    if applicable, their complementary BAI files
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
        if key == PG_TAG:
            continue
        assert file1.header[key] == file2.header[key],\
            "Unequal values for key %s: %s (%s), %s (%s)"\
            % (key, file1.filename, str(file1.header[key]),
               file2.filename, str(file2.header[key]))


def get_segment_qname(segment):
    """
    Return the qualified name of a read segment

    :param segment: read segment
    :type segment: pysam.libcalignedsegment.AlignedSegment
    :return: qualified name
    :rtype: str or unicode
    """
    return segment.qname


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
