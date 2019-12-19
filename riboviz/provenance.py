"""
Provenance-related utilities.
"""
import os.path
import git
import git.exc


def get_version(tool=os.path.basename(__file__)):
    """
    Get RiboViz version information.

    If this file is within the scope of a Git repository then a
    message including the given tool name, Git commit hash and date of
    HEAD is returned.

    If this file is not within the scope of a Git repository then an
    "unknown" message is returned.

    :param tool: Tool name
    :type tool: str or unicode
    :return: message
    :rtype: str or unicode
    """
    location = os.path.dirname(os.path.abspath(__file__))
    try:
        repository = git.Repo(location, search_parent_directories=True)
        sha = repository.head.object.hexsha
        time = repository.head.commit.authored_datetime
        version = "commit {} date {}".format(sha, str(time))
    except git.exc.InvalidGitRepositoryError:  # pylint: disable=E1101
        version = "unknown"
    return "{} version {}".format(tool, version)
