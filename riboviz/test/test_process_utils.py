"""
process_utils tests.
"""
import os.path
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
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
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
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    with open(log_out, 'w') as out, open(log_err, 'w') as err:
        try:
            process_utils.run_command(cmd, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 2
    assert lines[0] == path  # Output from ls
    assert lines[1] == path  # Output from ls
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert len(lines) == 1
    assert lines[0] == \
        "ls: cannot access 'no-such-file.txt': No such file or directory"


def test_run_command_log_out_error_one_file(log_out):
    """
    Test writing output and errors to a single log file.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    with open(log_out, "w") as out_err:
        try:
            process_utils.run_command(cmd, out_err, out_err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 3
    assert lines[0] == \
        "ls: cannot access 'no-such-file.txt': No such file or directory"
    assert lines[1] == path  # Output from ls
    assert lines[2] == path  # Output from ls


def test_run_command_log_out_err_alt(log_out, log_err):
    """
    Test another example of writing output and errors to log files.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    :param log_err: Error log file (via test fixture)
    :type log_err: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd = ["wc", "-l", path, "no-such-file.txt", path]
    with open(log_out, 'w') as out, open(log_err, 'w') as err:
        try:
            process_utils.run_command(cmd, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 3
    assert lines[0] == "%5d %s" % (num_lines, path)  # Output from wc
    assert lines[1] == "%5d %s" % (num_lines, path)  # Output from wc
    assert lines[2] == "%5d total" % (num_lines * 2)  # Output from wc
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert len(lines) == 1
    assert lines[0] == \
        "wc: no-such-file.txt: No such file or directory"


def test_run_command_log_out_error_one_file_alt(log_out):
    """
    Test another example of writing output and errors to a single log
    file.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd = ["wc", "-l", path, "no-such-file.txt", path]
    with open(log_out, "w") as out_err:
        try:
            process_utils.run_command(cmd, out_err, out_err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 4
    assert lines[0] == "%5d %s" % (num_lines, path)  # Output from wc
    assert lines[1] == "wc: no-such-file.txt: No such file or directory"
    assert lines[2] == "%5d %s" % (num_lines, path)  # Output from wc
    assert lines[3] == "%5d total" % (num_lines * 2)  # Output from wc


def test_run_redirect_command_stdout(redirect):
    """
    Test writing errors to stderr. This test ensures no unexpected
    exceptions are thrown.

    :param redirect: File for redirected output (via test fixture)
    :type redirect: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["cat", path, "no-such-file.txt"]
    try:
        process_utils.run_redirect_command(cmd, redirect)
    except AssertionError:
        pass
    # Compare path to captured redirect.
    with open(path) as expected, open(redirect) as actual:
        for line1, line2 in zip(expected, actual):
            assert line1 == line2


def test_run_redirect_command_log_err(log_err, redirect):
    """
    Test writing errors to a log file.

    :param log_err: Error log file (via test fixture)
    :type log_err: str or unicode
    :param redirect: File for redirected output (via test fixture)
    :type redirect: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["cat", path, "no-such-file.txt", path]
    with open(log_err, "w") as err:
        try:
            process_utils.run_redirect_command(cmd, redirect, err)
        except AssertionError:
            pass
    # Compare path to captured redirect.
    with open(path) as expected, open(redirect) as actual:
        for line1, line2 in zip(expected, actual):
            assert line1 == line2
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert len(lines) == 1
    assert lines[0] == \
        "cat: no-such-file.txt: No such file or directory"


def test_run_pipe_command_stdout_stderr():
    """
    Test writing output and errors to stdout and stderr where
    the first command in the pipeline logs some error.

    This test ensures no unexpected exceptions are thrown.
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
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
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    with open(log_out, 'w') as out, open(log_err, 'w') as err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 1
    assert lines[0] == str(num_lines * 2)  # Output from wc
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert len(lines) == 1
    assert lines[0] == "cat: no-such-file: No such file or directory"


def test_run_pipe_command_log_out_err_one_file(log_out):
    """
    Test writing output and errors to a single log file where
    the first command in the pipeline logs some error.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    with open(log_out, 'w') as out_err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out_err, out_err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 2
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert str(num_lines * 2) == lines[1]  # Output from wc


def test_run_pipe_command_stdout_stderr2():
    """
    Test writing output and errors to stdout and stderr where
    the second command in the pipeline logs some error.

    This test ensures no unexpected exceptions are thrown.
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
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
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l", "-x"]
    with open(log_out, 'w') as out, open(log_err, 'w') as err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 0  # Expect output to be empty
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert len(lines) == 3
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert lines[1] == "wc: invalid option -- 'x'"
    assert lines[2] == "Try 'wc --help' for more information."


def test_run_pipe_command_log_out_err_one_file2(log_out):
    """
    Test writing output and errors to a single log file where
    the second command in the pipeline logs some error.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l", "-x"]
    with open(log_out, 'w') as out_err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out_err, out_err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 3
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert lines[1] == "wc: invalid option -- 'x'"
    assert lines[2] == "Try 'wc --help' for more information."


def test_run_logged_command(log_out):
    """
    Test writing output and errors to a single log file.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    try:
        process_utils.run_logged_command(cmd, log_out)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 3
    assert lines[0] == \
        "ls: cannot access 'no-such-file.txt': No such file or directory"
    assert lines[1] == path  # Output from ls
    assert lines[2] == path  # Output from ls


def test_run_logged_redirect_command(log_err, redirect):
    """
    Test writing errors to a log file.

    :param log_err: Error log file (via test fixture)
    :type log_err: str or unicode
    :param redirect: File for redirected output (via test fixture)
    :type redirect: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["cat", path, "no-such-file.txt", path]
    try:
        process_utils.run_logged_redirect_command(cmd, redirect, log_err)
    except AssertionError:
        pass
    # Compare path to captured redirect.
    with open(path) as expected, open(redirect) as actual:
        for line1, line2 in zip(expected, actual):
            assert line1 == line2
    lines = [line.rstrip('\n') for line in open(log_err)]
    assert len(lines) == 1
    assert lines[0] == \
        "cat: no-such-file.txt: No such file or directory"


def test_run_logged_pipe_command_log(log_out):
    """
    Test writing output and errors to a single log file where
    the first command in the pipeline logs some error.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    try:
        process_utils.run_logged_pipe_command(cmd1, cmd2, log_out)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 2
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert str(num_lines * 2) == lines[1]  # Output from wc


def test_run_logged_pipe_command2(log_out):
    """
    Test writing output and errors to a single log file where
    the second command in the pipeline logs some error.

    :param log_out: Output log file (via test fixture)
    :type log_out: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l", "-x"]
    try:
        process_utils.run_logged_pipe_command(cmd1, cmd2, log_out)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(log_out)]
    assert len(lines) == 3
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert lines[1] == "wc: invalid option -- 'x'"
    assert lines[2] == "Try 'wc --help' for more information."
