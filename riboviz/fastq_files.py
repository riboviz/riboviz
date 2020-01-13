"""
FASTQ file-related utilities.
"""


FASTQ_NAME = "{}.fastq"
""" .fastq file name format string """
FASTQ_GZ_NAME = FASTQ_NAME + ".gz"
""" .fastq.gz file name format string """
FASTQ_FORMAT = "fastq"
""" Format string for use with Bio.SeqIO.write. """


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
        return FASTQ_GZ_NAME.format(tag)
    return FASTQ_NAME.format(tag)


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
