"""
:py:mod:`riboviz.process_utils` tests.
"""
import os.path
import tempfile
import pytest
from riboviz import process_utils
from riboviz import utils


@pytest.fixture(scope="function")
def tmp_stdout_file():
    """
    Create a temporary file with a ``log`` suffix.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_stdout_file = tempfile.mkstemp(prefix="tmp_stdout", suffix=".log")
    yield tmp_stdout_file
    if os.path.exists(tmp_stdout_file):
        os.remove(tmp_stdout_file)


@pytest.fixture(scope="function")
def tmp_stderr_file():
    """
    Create a temporary file with a ``log`` suffix.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_stderr_file = tempfile.mkstemp(prefix="tmp_stderr", suffix=".log")
    yield tmp_stderr_file
    if os.path.exists(tmp_stderr_file):
        os.remove(tmp_stderr_file)


@pytest.fixture(scope="function")
def tmp_redirect_file():
    """
    Create a temporary file with a ``txt`` suffix.

    :return: file
    :rtype: str or unicdo(dict, str or unicode)
    """
    _, tmp_redirect_file = tempfile.mkstemp(prefix="tmp_redirect_",
                                            suffix=".txt")
    yield tmp_redirect_file
    if os.path.exists(tmp_redirect_file):
        os.remove(tmp_redirect_file)


@pytest.fixture(scope="function")
def tmp_cmd_file():
    """
    Create a temporary file with a ``sh`` suffix.

    :return: file
    :rtype: str or unicdo(dict, str or unicode)
    """
    _, tmp_cmd_file = tempfile.mkstemp(prefix="tmp_cmd_", suffix=".sh")
    yield tmp_cmd_file
    if os.path.exists(tmp_cmd_file):
        os.remove(tmp_cmd_file)


def test_run_command_stdout_stderr():
    """
    Test :py:func:`riboviz.process_utils.run_command` using standard
    output and standard error.
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    try:
        process_utils.run_command(cmd)
    except AssertionError:
        pass


def test_run_command_log_out_err(tmp_stdout_file, tmp_stderr_file):
    """
    Test :py:func:`riboviz.process_utils.run_command` using files to
    capture standard output and standard error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_stderr_file: Error log file
    :type tmp_stderr_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    with open(tmp_stdout_file, 'w') as out, \
         open(tmp_stderr_file, 'w') as err:
        try:
            process_utils.run_command(cmd, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 2
    assert lines[0] == path  # Output from ls
    assert lines[1] == path  # Output from ls
    lines = [line.rstrip('\n') for line in open(tmp_stderr_file)]
    assert len(lines) == 1
    assert lines[0] == \
        "ls: cannot access 'no-such-file.txt': No such file or directory" \
        or lines[0] == \
        "ls: cannot access no-such-file.txt: No such file or directory"


def test_run_command_log_out_error_one_file(tmp_stdout_file):
    """
    Test :py:func:`riboviz.process_utils.run_command` using a single
    file to capture both standard output and standard error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    with open(tmp_stdout_file, "w") as out_err:
        try:
            process_utils.run_command(cmd, out_err, out_err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 3
    assert lines[0] == \
        "ls: cannot access 'no-such-file.txt': No such file or directory" \
        or lines[0] == \
        "ls: cannot access no-such-file.txt: No such file or directory"
    assert lines[1] == path  # Output from ls
    assert lines[2] == path  # Output from ls


def test_run_command_log_out_err_alt(tmp_stdout_file, tmp_stderr_file):
    """
    Test :py:func:`riboviz.process_utils.run_command` using files to
    capture standard output and standard error. Different commands are
    submitted to the operating system.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_stderr_file: Error log file
    :type tmp_stderr_file: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd = ["wc", "-l", path, "no-such-file.txt", path]
    with open(tmp_stdout_file, 'w') as out, \
        open(tmp_stderr_file, 'w') as err:
        try:
            process_utils.run_command(cmd, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 3
    assert lines[0] == "%5d %s" % (num_lines, path)  # Output from wc
    assert lines[1] == "%5d %s" % (num_lines, path)  # Output from wc
    assert lines[2] == "%5d total" % (num_lines * 2)  # Output from wc
    lines = [line.rstrip('\n') for line in open(tmp_stderr_file)]
    assert len(lines) == 1
    assert lines[0] == \
        "wc: no-such-file.txt: No such file or directory"


def test_run_command_log_out_error_one_file_alt(tmp_stdout_file):
    """
    Test :py:func:`riboviz.process_utils.run_command` using a single
    file to capture both standard output and standard error. Different
    commands are submitted to the operating system.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd = ["wc", "-l", path, "no-such-file.txt", path]
    with open(tmp_stdout_file, "w") as out_err:
        try:
            process_utils.run_command(cmd, out_err, out_err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 4
    assert lines[0] == "%5d %s" % (num_lines, path)  # Output from wc
    assert lines[1] == "wc: no-such-file.txt: No such file or directory"
    assert lines[2] == "%5d %s" % (num_lines, path)  # Output from wc
    assert lines[3] == "%5d total" % (num_lines * 2)  # Output from wc


def test_run_redirect_command_stdout(tmp_redirect_file):
    """
    Test :py:func:`riboviz.process_utils.run_redirect_command` using
    standard output.

    :param tmp_redirect_file: File for redirected output
    :type tmp_redirect_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["cat", path, "no-such-file.txt"]
    try:
        process_utils.run_redirect_command(cmd, tmp_redirect_file)
    except AssertionError:
        pass
    # Compare path to captured redirect.
    with open(path) as expected, open(tmp_redirect_file) as actual:
        for line1, line2 in zip(expected, actual):
            assert line1 == line2


def test_run_redirect_command_tmp_stderr_file(tmp_redirect_file,
                                              tmp_stderr_file):
    """
    Test :py:func:`riboviz.process_utils.run_redirect_command` using
    a file to capture standard error.

    :param tmp_redirect_file: File for redirected output
    :type tmp_redirect_file: str or unicode
    :param tmp_stderr_file: Error log file
    :type tmp_stderr_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["cat", path, "no-such-file.txt", path]
    with open(tmp_stderr_file, "w") as err:
        try:
            process_utils.run_redirect_command(cmd, tmp_redirect_file, err)
        except AssertionError:
            pass
    # Compare path to captured redirect.
    with open(path) as expected, open(tmp_redirect_file) as actual:
        for line1, line2 in zip(expected, actual):
            assert line1 == line2
    lines = [line.rstrip('\n') for line in open(tmp_stderr_file)]
    assert len(lines) == 1
    assert lines[0] == \
        "cat: no-such-file.txt: No such file or directory"


def test_run_pipe_command_stdout_stderr():
    """
    Test :py:func:`riboviz.process_utils.run_pipe_command` using
    standard output and standard error.
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    try:
        process_utils.run_pipe_command(cmd1, cmd2)
    except AssertionError:
        pass


def test_run_pipe_command_log_out_err(tmp_stdout_file, tmp_stderr_file):
    """
    Test :py:func:`riboviz.process_utils.run_pipe_command` using files
    to capture standard output and standard error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_stderr_file: Error log file
    :type tmp_stderr_file: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    with open(tmp_stdout_file, 'w') as out, open(tmp_stderr_file, 'w') as err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 1
    assert lines[0] == str(num_lines * 2)  # Output from wc
    lines = [line.rstrip('\n') for line in open(tmp_stderr_file)]
    assert len(lines) == 1
    assert lines[0] == "cat: no-such-file: No such file or directory"


def test_run_pipe_command_log_out_err_one_file(tmp_stdout_file):
    """
    Test :py:func:`riboviz.process_utils.run_pipe_command` using a
    single file to capture both standard output and standard error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    with open(tmp_stdout_file, 'w') as out_err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out_err, out_err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 2
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert str(num_lines * 2) == lines[1]  # Output from wc


def test_run_pipe_command_stdout_stderr_error():
    """
    Test :py:func:`riboviz.process_utils.run_pipe_command` using
    standard output and standard error, where the second command in
    the pipeline includes an error.
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l", "-x"]
    try:
        process_utils.run_pipe_command(cmd1, cmd2)
    except AssertionError:
        pass


def test_run_pipe_command_log_out_err_error(tmp_stdout_file,
                                            tmp_stderr_file):
    """
    Test :py:func:`riboviz.process_utils.run_pipe_command` using files
    to capture standard output and standard error, where the second
    command in the pipeline includes an error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_stderr_file: Error log file
    :type tmp_stderr_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l", "-x"]
    with open(tmp_stdout_file, 'w') as out, open(tmp_stderr_file, 'w') as err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out, err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 0  # Expect output to be empty
    lines = [line.rstrip('\n') for line in open(tmp_stderr_file)]
    assert len(lines) == 3
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert lines[1] == "wc: invalid option -- 'x'"
    assert lines[2] == "Try 'wc --help' for more information."


def test_run_pipe_command_log_out_err_one_file_error(tmp_stdout_file):
    """
    Test :py:func:`riboviz.process_utils.run_pipe_command` using a
    single file to capture both standard output and standard error,
    where the second command in the pipeline includes an error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l", "-x"]
    with open(tmp_stdout_file, 'w') as out_err:
        try:
            process_utils.run_pipe_command(cmd1, cmd2, out_err, out_err)
        except AssertionError:
            pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 3
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert lines[1] == "wc: invalid option -- 'x'"
    assert lines[2] == "Try 'wc --help' for more information."


def test_run_logged_command(tmp_stdout_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_command` using a
    single file to capture both standard output and standard error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    try:
        process_utils.run_logged_command(cmd, tmp_stdout_file)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 3
    assert lines[0] == \
        "ls: cannot access 'no-such-file.txt': No such file or directory" \
        or lines[0] == \
        "ls: cannot access no-such-file.txt: No such file or directory"
    assert lines[1] == path  # Output from ls
    assert lines[2] == path  # Output from ls


def test_run_logged_command_cmd_file(tmp_stdout_file, tmp_cmd_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_command` using a
    single file to capture both standard output and standard error and
    a file to capture commands sent to the operating system.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_cmd_file: Command file
    :type tmp_cmd_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    try:
        process_utils.run_logged_command(cmd, tmp_stdout_file,
                                         tmp_cmd_file)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 3
    assert lines[0] == \
        "ls: cannot access 'no-such-file.txt': No such file or directory" \
        or lines[0] == \
        "ls: cannot access no-such-file.txt: No such file or directory"
    assert lines[1] == path  # Output from ls
    assert lines[2] == path  # Output from ls
    with open(tmp_cmd_file) as f:
        actual_cmds = f.readlines()
    assert len(actual_cmds) == 1
    assert actual_cmds[0].rstrip('\n') == utils.list_to_str(cmd)


def test_run_logged_command_cmd_file_cmd_to_log(tmp_stdout_file,
                                                tmp_cmd_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_command` using a
    single file to capture both standard output and standard error and
    a file to capture commands sent to the operating system, where the
    command to be logged differs from that submitted.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_cmd_file: Command file
    :type tmp_cmd_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    cmd_to_log = ["ls", path, "'no-such-file.txt'", path]
    try:
        process_utils.run_logged_command(cmd,
                                         tmp_stdout_file,
                                         tmp_cmd_file,
                                         cmd_to_log=cmd_to_log)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 3
    assert lines[0] == \
        "ls: cannot access 'no-such-file.txt': No such file or directory" \
        or lines[0] == \
        "ls: cannot access no-such-file.txt: No such file or directory"
    assert lines[1] == path  # Output from ls
    assert lines[2] == path  # Output from ls
    with open(tmp_cmd_file) as f:
        actual_cmds = f.readlines()
    assert len(actual_cmds) == 1
    assert actual_cmds[0].rstrip('\n') == utils.list_to_str(cmd_to_log)


def test_run_logged_command_cmd_file_dry_run(tmp_stdout_file, tmp_cmd_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_command` using a
    single file to capture both standard output and standard error and
    a file to capture commands sent to the operating system, with
    the ``dry_run`` parameter set to ``True``.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_cmd_file: Command file
    :type tmp_cmd_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["ls", path, "no-such-file.txt", path]
    process_utils.run_logged_command(cmd, tmp_stdout_file,
                                     tmp_cmd_file, True)
    with open(tmp_stdout_file) as f:
        lines = f.readlines()
    assert len(lines) == 0
    with open(tmp_cmd_file) as f:
        actual_cmds = f.readlines()
    assert len(actual_cmds) == 1
    assert actual_cmds[0].rstrip('\n') == utils.list_to_str(cmd)


def test_run_logged_redirect_command(tmp_stderr_file, tmp_redirect_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_redirect_command`
    using a file to capture standard error.

    :param tmp_stderr_file: Error log file
    :type tmp_stderr_file: str or unicode
    :param tmp_redirect_file: File for redirected output
    :type tmp_redirect_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["cat", path, "no-such-file.txt", path]
    try:
        process_utils.run_logged_redirect_command(cmd, tmp_redirect_file, tmp_stderr_file)
    except AssertionError:
        pass
    # Compare path to captured redirect.
    with open(path) as expected, open(tmp_redirect_file) as actual:
        for line1, line2 in zip(expected, actual):
            assert line1 == line2
    lines = [line.rstrip('\n') for line in open(tmp_stderr_file)]
    assert len(lines) == 1
    assert lines[0] == \
        "cat: no-such-file.txt: No such file or directory"


def test_run_logged_redirect_command_cmd_file(
        tmp_stderr_file, tmp_redirect_file, tmp_cmd_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_redirect_command`
    using a file to capture standard error and a file to capture
    commands sent to the operating system.

    :param tmp_stderr_file: Error log file
    :type tmp_stderr_file: str or unicode
    :param tmp_redirect_file: File for redirected output
    :type tmp_redirect_file: str or unicode
    :param tmp_cmd_file: Command file
    :type tmp_cmd_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["cat", path, "no-such-file.txt", path]
    try:
        process_utils.run_logged_redirect_command(cmd,
                                                  tmp_redirect_file,
                                                  tmp_stderr_file,
                                                  tmp_cmd_file)
    except AssertionError:
        pass
    # Compare path to captured redirect.
    with open(path) as expected, open(tmp_redirect_file) as actual:
        for line1, line2 in zip(expected, actual):
            assert line1 == line2
    lines = [line.rstrip('\n') for line in open(tmp_stderr_file)]
    assert len(lines) == 1
    assert lines[0] == \
        "cat: no-such-file.txt: No such file or directory"
    with open(tmp_cmd_file) as f:
        actual_cmds = f.readlines()
    assert len(actual_cmds) == 1
    expected_cmd = "%s > %s" % (utils.list_to_str(cmd), tmp_redirect_file)
    assert actual_cmds[0].rstrip('\n') == expected_cmd


def test_run_logged_redirect_command_cmd_file_dry_run(
        tmp_stderr_file, tmp_redirect_file, tmp_cmd_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_redirect_command`
    using a file to capture standard error and a file to capture
    commands sent to the operating system, with the ``dry_run``
    parameter set to ``True``.

    :param tmp_stderr_file: Error log file
    :type tmp_stderr_file: str or unicode
    :param tmp_redirect_file: File for redirected output
    :type tmp_redirect_file: str or unicode
    :param tmp_cmd_file: Command file
    :type tmp_cmd_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd = ["cat", path, "no-such-file.txt", path]
    process_utils.run_logged_redirect_command(cmd,
                                              tmp_redirect_file,
                                              tmp_stderr_file,
                                              tmp_cmd_file,
                                              True)
    with open(tmp_redirect_file) as f:
        lines = f.readlines()
    assert len(lines) == 0
    with open(tmp_stderr_file) as f:
        lines = f.readlines()
    assert len(lines) == 0
    with open(tmp_cmd_file) as f:
        actual_cmds = f.readlines()
    assert len(actual_cmds) == 1
    expected_cmd = "%s > %s" % (utils.list_to_str(cmd), tmp_redirect_file)
    assert actual_cmds[0].rstrip('\n') == expected_cmd


def test_run_logged_pipe_command_log(tmp_stdout_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_pipe_command`
    using a single file to capture both standard output and standard
    error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    try:
        process_utils.run_logged_pipe_command(cmd1, cmd2, tmp_stdout_file)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 2
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert str(num_lines * 2) == lines[1]  # Output from wc


def test_run_logged_pipe_command_log_cmd_file(tmp_stdout_file,
                                              tmp_cmd_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_pipe_command`
    using a single file to capture both standard output and standard
    error and a file to capture commands sent to the operating
    system.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_cmd_file: Command file
    :type tmp_cmd_file: str or unicode
    """
    path = os.path.realpath(__file__)
    num_lines = len([line for line in open(path)])
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    try:
        process_utils.run_logged_pipe_command(cmd1,
                                              cmd2,
                                              tmp_stdout_file,
                                              tmp_cmd_file)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 2
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert str(num_lines * 2) == lines[1]  # Output from wc
    with open(tmp_cmd_file) as f:
        actual_cmds = f.readlines()
    assert len(actual_cmds) == 1
    expected_cmd = "%s | %s" % (utils.list_to_str(cmd1),
                                utils.list_to_str(cmd2))
    assert actual_cmds[0].rstrip('\n') == expected_cmd


def test_run_logged_pipe_command_log_cmd_file_dry_run(tmp_stdout_file,
                                                      tmp_cmd_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_pipe_command`
    using a single file to capture both standard output and standard
    error and a file to capture commands sent to the operating
    system, with the ``dry_run`` parameter set to ``True``.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    :param tmp_cmd_file: Command file
    :type tmp_cmd_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l"]
    process_utils.run_logged_pipe_command(cmd1,
                                          cmd2,
                                          tmp_stdout_file,
                                          tmp_cmd_file,
                                          True)
    with open(tmp_stdout_file) as f:
        lines = f.readlines()
    assert len(lines) == 0
    with open(tmp_cmd_file) as f:
        actual_cmds = f.readlines()
    assert len(actual_cmds) == 1
    expected_cmd = "%s | %s" % (utils.list_to_str(cmd1),
                                utils.list_to_str(cmd2))
    assert actual_cmds[0].rstrip('\n') == expected_cmd


def test_run_logged_pipe_command_error(tmp_stdout_file):
    """
    Test :py:func:`riboviz.process_utils.run_logged_pipe_command`
    using a single file to capture both standard output and standard
    error, where the first command in the pipeline includes an error.

    :param tmp_stdout_file: Output log file
    :type tmp_stdout_file: str or unicode
    """
    path = os.path.realpath(__file__)
    cmd1 = ["cat", path, "no-such-file", path]
    cmd2 = ["wc", "-l", "-x"]
    try:
        process_utils.run_logged_pipe_command(cmd1, cmd2, tmp_stdout_file)
    except AssertionError:
        pass
    lines = [line.rstrip('\n') for line in open(tmp_stdout_file)]
    assert len(lines) == 3
    assert lines[0] == "cat: no-such-file: No such file or directory"
    assert lines[1] == "wc: invalid option -- 'x'"
    assert lines[2] == "Try 'wc --help' for more information."
