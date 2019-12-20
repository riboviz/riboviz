"""
Provenance-related utilities.
"""
from datetime import datetime
import os.path
import git
import git.exc


def get_version(file_path=__file__):
    """
    Get RiboViz version information.

    If the file is within the scope of a Git repository then a
    message including the given file name, Git commit hash and date of
    HEAD is returned.

    If the file is not within the scope of a Git repository then an
    "unknown" message is returned.

    :param file_path: Python file path
    :type file_path: str or unicode
    :return: message
    :rtype: str or unicode
    """
    file_name = os.path.basename(file_path)
    location = os.path.dirname(os.path.abspath(file_path))
    try:
        repository = git.Repo(location, search_parent_directories=True)
        sha = repository.head.object.hexsha
        time = repository.head.commit.authored_datetime
        version = "commit {} date {}".format(sha, str(time))
    except git.exc.InvalidGitRepositoryError:  # pylint: disable=E1101
        version = "unknown"
    return "{} version {}".format(file_name, version)


def write_metadata_header(file_path, data_file):
    """
    Write a metadata header to a file with RiboViz date and version
    information.

    :param file_path: Path to Python file creating the data file.
    :type file_path: str or unicode
    :param data_file: Path to file into which header is to be written.
    :type data_file: str or unicode
    """
    with open(data_file, 'w') as f:
        f.write("# Created by: RiboViz\n")
        f.write("# Date: {}\n".format(datetime.today()))
        f.write("# Component/version: {}\n".format(get_version(file_path)))
