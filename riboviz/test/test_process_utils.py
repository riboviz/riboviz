"""
process_utils tests.
"""
import os.path
import subprocess
import sys
import tempfile
import pytest
from riboviz import process_utils


@pytest.fixture(scope="function")
def log_out():
    """
    Create a temporary file for stdout logs and delete when done.

    :return: file
    :rtype: str or unicdo(dict, str or unicode)
    """
    _, out_file = tempfile.mkstemp(prefix="log_out_", suffix=".txt")
    yield out_file
    if os.path.exists(out_file):
        os.remove(out_file)


@pytest.fixture(scope="function")        
def log_err():
    """
    Create a temporary file for stderr logs and delete when done.

    :return: file
    :rtype: str or unicdo(dict, str or unicode)
    """
    _, err_file = tempfile.mkstemp(prefix="log_err_", suffix=".txt")
    yield err_file
    if os.path.exists(err_file):
        os.remove(err_file)


@pytest.fixture(scope="function")
def redirect():
    """
    Create a temporary file to capture redirected output and delete
    when done.

    :return: file
    :rtype: str or unicdo(dict, str or unicode)
    """
    _, redirect_file = tempfile.mkstemp(prefix="redirect_", suffix=".txt")
    yield redirect_file
    if os.path.exists(redirect_file):
        os.remove(redirect_file)


def test_run_command_stdout_stderr():
    """
    Test writing output and errors to stdout and stderr. This
    test ensures no unexpected exceptions are thrown.
    """
    cmd = ["ls", "-l", "riboviz", "xxx", "rscripts"]
    try:
        process_utils.run_command(cmd)
    except AssertionError:
        pass
def test_run_command_log_out_err(log_out, log_err):
    """
    Test writing output and errors to log files.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    :param log_err: Error log file (via test fixture)
    :type log_err: str or unicode
    """
    cmd = ["ls", "-l", "riboviz", "xxx", "rscripts"]
    with open(log_out, 'w') as out, open(log_err, 'w') as err:
        try:
            process_utils.run_command(cmd, out, err)
        except AssertionError:
            pass
    # TODO validate log_out
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert 1 == len(lines)
    assert "ls: cannot access 'xxx': No such file or directory" == lines[0]
def test_run_command_log_out_error_one_file(log_out):
    """
    Test writing output and errors to a single log file.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    cmd = ["ls", "-l", "riboviz", "xxx", "rscripts"]
    with open(log_out, "w") as out_err:
        try:
            process_utils.run_command(cmd, out_err, out_err)
        except AssertionError:
            pass
    # TODO validate log_out
    # lines = [line.rstrip('\n') for line in open(log_out)]
    # assert 1 == len(lines)
    # assert "ls: cannot access 'xxx': No such file or directory" == lines[0]

def test_run_command_stdout_stderr2():
    """
    Test writing output and errors to stdout and stderr. This
    test ensures no unexpected exceptions are thrown.
    """
    cmd = ["du", "-sb", "README.md", "xxx", "LICENSE"]
    try:
        process_utils.run_command(cmd)
    except AssertionError:
        pass
def test_run_command_log_out_err2(log_out, log_err):
    """
    Test writing output and errors to log files.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    :param log_err: Error log file (via test fixture)
    :type log_err: str or unicode
    """
    cmd = ["du", "-sb", "README.md", "xxx", "LICENSE"]
    with open(log_out, 'w') as out, open(log_err, 'w') as err:
        try:
            process_utils.run_command(cmd, out, err)
        except AssertionError:
            pass
    # TODO validate log_out
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert 1 == len(lines)
    assert "du: cannot access 'xxx': No such file or directory" == lines[0]
def test_run_command_log_out_error_one_file2(log_out):
    """
    Test writing output and errors to a single log file.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    cmd = ["du", "-sb", "README.md", "xxx", "LICENSE"]
    with open(log_out, "w") as out_err:
        try:
            process_utils.run_command(cmd, out_err, out_err)
        except AssertionError:
            pass
    # TODO validate log_out
    # lines = [line.rstrip('\n') for line in open(log_out)]
    # assert 1 == len(lines)
    # assert "du: cannot access 'xxx': No such file or directory" == lines[0]

def test_run_redirect_command_stdout(redirect):
    """
    Test writing errors to stderr. This test ensures no unexpected
    exceptions are thrown.

    :param redirect: File for redirected output (via test fixture)
    :type redirect: str or unicode
    """
    cmd = ["cat", "README.md", "xxx", "README.md"]
    try:
        process_utils.run_redirect_command(cmd, redirect)
    except AssertionError:
        pass
    # TODO validate redirect
def test_run_redirect_command_log_err(log_err, redirect):
    """
    Test writing errors to a log file.

    :param log_err: Error log file (via test fixture)
    :type log_err: str or unicode
    :param redirect: File for redirected output (via test fixture)
    :type redirect: str or unicode
    """
    cmd = ["cat", "README.md", "xxx", "README.md"]
    with open(log_err, "w") as err:
        try:
            process_utils.run_redirect_command(cmd, redirect, err)
        except AssertionError:
            pass
    # TODO validate redirect
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert 1 == len(lines)
    assert "cat: xxx: No such file or directory" == lines[0]
        
def test_run_pipe_command_stdout_stderr():
    """
    Test writing output and errors to stdout and stderr where
    the first command in the pipeline logs some error.

    This test ensures no unexpected exceptions are thrown.
    """
    cmd1 = ["cat", "README.md", "xxx", "LICENSE"]
    cmd2 = ["wc", "-l"]
    try:
        process_utils.run_pipe_command(cmd1, cmd2)
    except AssertionError:
        pass
def test_run_pipe_command_log_out_err(log_out, log_err):
    """
    Test writing output and errors to stdout and stderr where
    the first command in the pipeline logs some error.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    :param log_err: Error log file (via test fixture)
    :type log_err: str or unicode
    """
    cmd1 = ["cat", "README.md", "xxx", "LICENSE"]
    cmd2 = ["wc", "-l"]
    with open(log_out, 'w') as out, open(log_err, 'w') as err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out, err)
        except AssertionError:
            pass
    # TODO validate output
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert 1 == len(lines)
    assert "cat: xxx: No such file or directory" == lines[0]
def test_run_pipe_command_log_out_err_one_file(log_out):
    """
    Test writing output and errors to a single log file where
    the first command in the pipeline logs some error.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    cmd1 = ["cat", "README.md", "xxx", "LICENSE"]
    cmd2 = ["wc", "-l"]
    with open(log_out, 'w') as out_err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out_err, out_err)
        except AssertionError:
            pass
    # TODO validate output
    # lines = [line.rstrip('\n') for line in open(log_out)]
    # assert 1 == len(lines)
    # assert "cat: xxx: No such file or directory" == lines[0]

def test_run_pipe_command_stdout_stderr2():
    """
    Test writing output and errors to stdout and stderr where
    the second command in the pipeline logs some error.

    This test ensures no unexpected exceptions are thrown.
    """
    cmd1 = ["cat", "README.md", "xxx", "LICENSE"]
    cmd2 = ["wc", "-l", "-x"]
    try:
        process_utils.run_pipe_command(cmd1, cmd2)
    except AssertionError:
        pass
def test_run_pipe_command_log_out_err2(log_out, log_err):
    """
    Test writing output and errors to log files where
    the second command in the pipeline logs some error.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    :param log_err: Error log file (via test fixture)
    :type log_err: str or unicode
    """
    cmd1 = ["cat", "README.md", "xxx", "LICENSE"]
    cmd2 = ["wc", "-l", "-x"]
    with open(log_out, 'w') as out, open(log_err, 'w') as err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_err)]
    # TODO validate output
    assert 3 == len(lines)
    assert "cat: xxx: No such file or directory" == lines[0]
    assert "wc: invalid option -- 'x'" == lines[1]
    assert "Try 'wc --help' for more information." == lines[2]

def test_run_pipe_command_log_out_err_one_file2(log_out):
    """
    Test writing output and errors to a single log file where
    the second command in the pipeline logs some error.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    cmd1 = ["cat", "README.md", "xxx", "LICENSE"]
    cmd2 = ["wc", "-l", "-x"]
    with open(log_out, 'w') as out_err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out_err, out_err)
        except AssertionError:
            pass
    # TODO validate output
    # lines = [line.rstrip('\n') for line in open(log_out)]
    # assert 3 == len(lines)
    # assert "cat: xxx: No such file or directory" == lines[0]
    # assert "wc: invalid option -- 'x'" == lines[1]
    # assert "Try 'wc --help' for more information." == lines[2]
