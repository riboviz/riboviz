"""
Provenance-related utilities.
"""
from datetime import datetime
import os.path
from io import StringIO
import git
import git.exc


def get_version(file_path=__file__):
    """
    Get RiboViz version information.

    If the file is within the scope of a Git repository then a
    message including the Git commit hash and date of HEAD is
    returned.

    If the file is not within the scope of a Git repository then an
    "unknown" message is returned.

    :param file_path: Python file path
    :type file_path: str or unicode
    :return: messsage
    :rtype: str or unicode
    """
    location = os.path.dirname(os.path.abspath(file_path))
    try:
        repository = git.Repo(location, search_parent_directories=True)
        sha = repository.head.object.hexsha
        time = repository.head.commit.authored_datetime
        version = "commit {} date {}".format(sha, str(time))
    except git.exc.InvalidGitRepositoryError:  # pylint: disable=E1101
        version = "unknown"
    return version


def write_provenance(file_handle, file_path, prefix="# ", eol="\n"):
    """
    Write a provenance header to a file handle with RiboViz date and
    version information.

    :param file_path: Path to Python file creating the data file.
    :type file_path: str or unicode
    :param data_file: Path to file into which header is to be written.
    :type data_file: str or unicode
    :param prefix: Prefix for each line (e.g. a comment specifier)
    :type prefix: str or unicode
    :param eol: End of line character
    :type eol: str or unicode
    """
    file_handle.write("{}Created by: RiboViz{}".format(prefix, eol))
    file_handle.write("{}Date: {}{}".format(prefix, datetime.today(), eol))
    import __main__
    if hasattr(__main__, "__file__"):
        file_handle.write("{}Command-line tool: {}{}".format(
            prefix, __main__.__file__, eol))
    file_handle.write("{}File: {}{}".format(prefix, file_path, eol))
    file_handle.write("{}Version: {}{}".format(
        prefix, get_version(file_path), eol))


def write_provenance_header(file_path, data_file, prefix="# "):
    """
    Write a provenance header to a file with RiboViz date and version
    information.

    :param file_path: Path to Python file invoking this function
    :type file_path: str or unicode
    :param data_file: Path to file into which header is to be written.
    :type data_file: str or unicode
    :param prefix: Prefix for each line (e.g. a comment specifier)
    :type prefix: str or unicode
    """
    with open(data_file, 'w') as f:
        write_provenance(f, file_path, prefix)


def get_provenance_str(file_path, eol="\n"):
    """
    Get provenance information in a string with RiboViz date and
    version information.

    :param file_path: Path to Python file invoking this function
    :type file_path: str or unicode
    :param eol: End of line character
    :type eol: str or unicode
    """
    with StringIO() as s:
        write_provenance(s, file_path, "", eol)
        provenance = s.getvalue()
    return provenance
