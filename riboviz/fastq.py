"""
FASTQ file-related utilities.
"""
import gzip
import os.path
from Bio import SeqIO


EXTENSIONS = [".fq", ".fastq"]
""" FASTQ file extensions """
FASTQ_FORMAT = "{}.fastq"
""" .fastq file name format string """
FASTQ_GZ_FORMAT = FASTQ_FORMAT + ".gz"
""" .fastq.gz file name format string """


def is_fastq_gz(file_name):
    """
    Does the given file end with .gz or .GZ or .gzip or .GZIP?

    :param file_name: File name
    :type file_name: str or unicode
    :return: True if file_name ends with .gz or .GZ, False
    otherwise
    :rtype: bool
    """
    _, ext = os.path.splitext(os.path.basename(file_name))
    return ext.lower() in [".gz", ".gzip"]


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
        return FASTQ_GZ_FORMAT.format(tag)
    return FASTQ_FORMAT.format(tag)


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


def count_records(file_name):
    """
    Count number of records in FASTQ file. Both FASTQ and FASTQ.GZ
    files are handled.

    :param file_name: File name
    :type file_name: str or unicode
    :return: number of records
    :rtype: int
    """

    num_records = 0
    if is_fastq_gz(file_name):
        open_file = gzip.open
    else:
        open_file = open
    with open_file(file_name, "rt") as f:
        for _ in SeqIO.parse(f, "fastq"):
            num_records = num_records + 1
    return num_records
