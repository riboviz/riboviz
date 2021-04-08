"""
html-related constants and functions.
"""
import subprocess

HTML_EXT = "html"
""" File extension. """
HTML_FORMAT = "{}_output_report." + HTML_EXT
""" File name format. """


def equal_html(file1, file2):
    """
    Compare two html files for equality, using ``utils.equal_file_names``.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: If the file contents differ
    :raise Exception: If problems arise when loading the files or \
    running ``utils.equal_file_names``
    """
    # TODO implement in-Python comparison.
    cmd = [utils.equal_file_names(file1, file2)]
    return_code = subprocess.call(cmd)
    assert return_code == 0,\
        "Non-zero return code (%d) from %s" % (
            return_code, ' '.join(map(str, cmd)))
