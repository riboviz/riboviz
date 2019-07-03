"""
pytest plugin file.

See https://docs.pytest.org/en/latest/writing_plugins.html
"""
import os.path
import pytest

EXPECTED = "--expected"
""" Command-line flag for directory with expected data files """
ACTUAL = "--actual"
"""
Command-line flag for directory to be validated against directory with
expected data files
"""


def pytest_addoption(parser):
    """
    pytest configutation hook.

    See
    https://docs.pytest.org/en/latest/reference.html#_pytest.hookspec.pytest_addoption

    :param parser: command-line parser
    :type parser: _pytest.config.argparsing.Parser
    """
    parser.addoption(
        EXPECTED,
        action="store",
        required=True,
        help="Directory with expected data files"
    )
    parser.addoption(
        ACTUAL,
        action="store",
        default="vignette",
        help="Directory to be validated against directory with expected data files"
    )


@pytest.fixture
def command_option(request):
    """
    Gets values for "--expected" and "--actual" command-line
    options.

    :param request:
    :type request: _pytest.fixtures.SubRequest
    :return: directory with expected data files, directory to be
    validated against directory with expected data files
    :rtype: tuple(str or unicode, str or unicode)
    :raise AssertionError: if either option has a value that is
    not a directory
    """
    expected = request.config.getoption(EXPECTED)
    actual = request.config.getoption(ACTUAL)
    assert os.path.exists(expected) and os.path.isdir(expected),\
        "No such directory: %s" % expected
    assert os.path.exists(expected) and os.path.isdir(actual),\
        "No such directory: %s" % actual
    return (expected, actual)
