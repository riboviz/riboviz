"""
pytest plugin file for regression tests. See
https://docs.pytest.org/en/latest/writing_plugins.html

This allows pytest to take in additional command-line parameters to
pass onto the regression test modules:

* ``--expected=<DIRECTORY>``: Directory with expected data files,
  against which files in ``vignette/`` will be checked.
* ``--skip-workflow``: Workflow will not be run prior to checking data
  files.
* `--check-index-tmp`: Check index and temporary files (default is
  that only the output files are checked).
"""
import os.path
import pytest

EXPECTED = "--expected"
""" Directory with expected data files command-line flag."""
SKIP_WORKFLOW = "--skip-workflow"
""" Do not run workflow command-line flag. """
CHECK_INDEX_TMP = "--check-index-tmp"
""" Check index and temporary files command-line flag. """


def pytest_addoption(parser):
    """
    pytest configuration hook. See
    https://docs.pytest.org/en/latest/reference.html#_pytest.hookspec.pytest_addoption

    :param parser: command-line parser
    :type parser: _pytest.config.argparsing.Parser
    """
    parser.addoption(EXPECTED,
                     action="store",
                     required=True,
                     help="Directory with expected data files")
    parser.addoption(SKIP_WORKFLOW,
                     action="store_true",
                     required=False,
                     help="Do not run workflow")
    parser.addoption(CHECK_INDEX_TMP,
                     action="store_true",
                     required=False,
                     help="Check index and temporary files")


@pytest.fixture
def expected_fixture(request):
    """
    Gets value for ``--expected`` command-line option.

    :param request: request
    :type request: _pytest.fixtures.SubRequest
    :return: directory
    :rtype: str or unicode
    :raise AssertionError: if the option has a value that is \
    not a directory
    """
    expected_dir = request.config.getoption(EXPECTED)
    assert os.path.exists(expected_dir) and os.path.isdir(expected_dir),\
        "No such directory: %s" % expected_dir
    return expected_dir


@pytest.fixture(scope="module")
def skip_workflow_fixture(request):
    """
    Gets value for ``--skip-workflow`` command-line option.

    :param request: request
    :type request: _pytest.fixtures.SubRequest
    :return: flag
    :rtype: bool
    """
    return request.config.getoption(SKIP_WORKFLOW)


@pytest.fixture(scope="module")
def skip_index_tmp_fixture(request):
    """
    Gets value for `--check-index-tmp` command-line option. If
    ``False``, or undefined, invokes ``pytest.skip`` to skip
    test.

    :param request: request
    :type request: _pytest.fixtures.SubRequest
    :return: flag
    :rtype: bool
    """
    if not request.config.getoption(CHECK_INDEX_TMP):
        pytest.skip('Skipped index and temporary files tests')
