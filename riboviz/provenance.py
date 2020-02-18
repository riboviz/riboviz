"""
Provenance-related functions.
"""
from datetime import datetime
import os.path
from io import StringIO
import git
import git.exc


def get_version(file_path=__file__):
    """
    Get version information about ``file_path`` using the ``git``
    package.

    If ``file_path`` is within the scope of a Git repository then a
    string including the Git commit hash and date of ``HEAD`` is
    returned. The message has format ``commit <HASH> date <DATE>``.

    If ``file_path`` is not within the scope of a Git repository then
    the string ``unknown`` is returned.

    :param file_path: File path
    :type file_path: str or unicode
    :return: Version information
    :rtype: str or unicode
    """
    location = os.path.dirname(os.path.abspath(file_path))
    try:
        repository = git.Repo(location,
                              search_parent_directories=True)
        sha = repository.head.object.hexsha
        time = repository.head.commit.authored_datetime
        version = "commit {} date {}".format(sha, str(time))
    except git.exc.InvalidGitRepositoryError:  # pylint: disable=E1101
        version = "unknown"
    return version


def write_provenance(file_handle, file_path, prefix="# ", eol="\n"):
    """
    Write a provenance header to a file including version information
    about ``file_path`` obtained using the ``git`` package. See also
    :py:func:`get_version`.

    The header has form::

        <prefix> Created by: RiboViz
        <prefix> Date: <YYYY>-<MM>-<DD> <HH>:<MM>:<SS>.<SSSSS>
        <prefix> Command-line tool: <MAIN_MODULE>
        <prefix> File: <file_path>
        <prefix> Version: commit <HASH> date <DATE>

    where ``<MAIN_MODULE>`` is the Python module defined as
    ``__main__`` (i.e. the identify of the currently running
    program).

    :param file_handle: File handle to which content is to be written
    :type file_handle: _io.TextIOWrapper
    :param file_path: File path
    :type file_path: str or unicode
    :param prefix: Prefix for each line e.g. a comment symbol
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


def write_provenance_header(file_path, provenance_file, prefix="# "):
    """
    Write a provenance header to a file including version information
    about ``file_path`` obtained using the ``git`` package. See
    :py:func:`write_provenance`.

    :param file_path: File path
    :type file_path: str or unicode
    :param provenance_file: File to which content is to be written
    :type provenance: _io.TextIOWrapper
    :param prefix: Prefix for each line e.g. a comment symbol
    :type prefix: str or unicode
    """
    with open(provenance_file, 'w') as f:
        write_provenance(f, file_path, prefix)


def get_provenance_str(file_path, eol="\n"):
    """
    Get a provenance header including version information about
    ``file_path`` obtained using the ``git`` package. See
    :py:func:`write_provenance`.

    :param file_path: File path
    :type file_path: str or unicode
    :param eol: End of line character
    :type eol: str or unicode
    :return: Provenance header as a string
    :rtype: str or unicode
    """
    with StringIO() as s:
        write_provenance(s, file_path, "", eol)
        provenance = s.getvalue()
    return provenance
