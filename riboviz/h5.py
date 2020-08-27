"""
H5-related constants and functions.
"""
import subprocess

H5_EXT = "h5"
""" File extension. """
H5_FORMAT = "{}." + H5_EXT
""" File name format. """


def equal_h5(file1, file2):
    """
    Compare two H5 files for equality, using the ``h5diff``
    command-line tool.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: If the file contents differ
    :raise Exception: If problems arise when loading the files or \
    running ``h5diff``
    """
    # TODO implement in-Python comparison.
    cmd = ["h5diff", "-q", file1, file2]
    return_code = subprocess.call(cmd)
    assert return_code == 0,\
        "Non-zero return code (%d) from %s" % (
            return_code, ' '.join(map(str, cmd)))
