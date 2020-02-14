"""
FASTQ file-related utilities.
"""
import gzip
import os.path
from Bio import SeqIO
from riboviz import utils

FASTQ_EXT = "fastq"
""" fastq file extension """
FQ_EXT = "fq"
""" fq file extension """
FASTQ_FORMAT = "{}." + FASTQ_EXT
""" fastq file name format string """
FQ_FORMAT = "{}." + FQ_EXT
""" fq file name format string """
FASTQ_GZ_EXT = FASTQ_EXT + ".gz"
""" fastq.gz file extension """
FQ_GZ_EXT = FQ_EXT + ".gz"
""" fq.gz file extension """
FASTQ_GZ_FORMAT = "{}." + FASTQ_GZ_EXT
""" fastq.gz file name format string """
FQ_GZ_FORMAT = "{}." + FQ_GZ_EXT
""" fq.gz file name format string """

FASTQ_EXTS = [FASTQ_EXT, FQ_EXT]
""" FASTQ file extensions """
FASTQ_GZ_EXTS = [FASTQ_GZ_EXT, FQ_GZ_EXT]
""" FASTQ GZ file extensions """
FASTQ_ALL_EXTS = FASTQ_EXTS + FASTQ_GZ_EXTS
""" All FASTQ file extensions """

FASTQ_FORMATS = {FASTQ_EXT: FASTQ_FORMAT,
                 FQ_EXT: FQ_FORMAT,
                 FASTQ_GZ_EXT: FASTQ_GZ_FORMAT,
                 FQ_GZ_EXT: FQ_GZ_FORMAT}
""" Map from file extensions to file name formats """


def is_fastq_gz(file_name):
    """
    Does the given file end with .gz or .GZ or .gzip or .GZIP?

    :param file_name: File name
    :type file_name: str or unicode
    :return: True if file_name ends with .gz or .GZ
    :rtype: bool
    """
    ext = utils.get_file_ext(file_name)
    return ext.lower() in FASTQ_GZ_EXTS


def strip_fastq_gz(file_name):
    """
    If the given file ends with .gz or .GZ then remove the .gz or .GZ
    extension. If the file doesn't end with this extension then the
    original name is returned.

    :param file_name: File name
    :type file_name: str or unicode
    :return: File name without .gz or .GZ extension
    :rtype: str or unicode
    """
    if is_fastq_gz(file_name):
        return os.path.splitext(file_name)[0]
    return file_name


def count_sequences(file_name):
    """
    Count number of sequences in a FASTQ file. Both FASTQ and FASTQ.GZ
    files are handled.

    :param file_name: File name
    :type file_name: str or unicode
    :return: number of sequences
    :rtype: int
    """
    num_sequences = 0
    if is_fastq_gz(file_name):
        open_file = gzip.open
    else:
        open_file = open
    with open_file(file_name, "rt") as f:
        for _ in SeqIO.parse(f, "fastq"):
            num_sequences = num_sequences + 1
    return num_sequences


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
