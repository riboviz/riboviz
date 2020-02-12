"""
pytest plugin file.

See https://docs.pytest.org/en/latest/writing_plugins.html
"""
import os.path
import pytest

EXPECTED = "--expected"
""" Command-line flag for directory with expected data files """
SKIP_WORKFLOW = "--skip-workflow"
""" Command-line flag to specify that workflow should not be run """
CHECK_INDEX_TMP = "--check-index-tmp"
"""
Command-line flag to specify that index and temprary files should be
checked too (default is that only output files are checked.
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
        SKIP_WORKFLOW,
        action="store_true",
        required=False,
        help="Workflow should not be run"
    )
    parser.addoption(
        CHECK_INDEX_TMP,
        action="store_true",
        required=False,
        help="Index and temporary files should be checked too"
    )


@pytest.fixture
def expected(request):
    """
    Gets value for "--expected" command-line option.

    :param request: request
    :type request: _pytest.fixtures.SubRequest
    :return: directory with expected data files
    :rtype: str or unicode
    :raise AssertionError: if the option has a value that is
    not a directory
    """
    expected_dir = request.config.getoption(EXPECTED)
    assert os.path.exists(expected_dir) and os.path.isdir(expected_dir),\
        "No such directory: %s" % expected_dir
    return expected_dir


@pytest.fixture(scope="module")
def skip_workflow(request):
    """
    Gets value for "--skip-workflow" command-line option.

    :param request: request
    :type request: _pytest.fixtures.SubRequest
    :return: whether workflow should be skipped
    :rtype: bool
    """
    return request.config.getoption(SKIP_WORKFLOW)


@pytest.fixture(scope="module")
def skip_index_tmp(request):
    """
    Gets value for "--check-index-tmp" command-line option and, if
    False, or undefined, invokes "pytest.skip" to skip a test.

    :param request: request
    :type request: _pytest.fixtures.SubRequest
    :return: whether temporary files should be checked
    :rtype: bool
    """
    if not request.config.getoption(CHECK_INDEX_TMP):
        pytest.skip('Checking output files only')
