"""
Provenance-related utilities.
"""
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
