"""
Python ``subprocess``-related functions.
"""
import subprocess
import sys
from riboviz import utils


def run_command(cmd, out=sys.stdout, err=sys.stderr):
    """
    Run operating system command via Python ``subprocess``.

    :param cmd: Commnand and arguments
    :type cmd: list(str or unicode)
    :param out: Standard output desination (``sys.stdout`` or file)
    :type out: _io.TextIOWrapper
    :param err: Standard error desination (``sys.stderr`` or file)
    :type err: _io.TextIOWrapper
    :raise AssertionError: If the command returns a non-zero exit code
    """
    exit_code = subprocess.call(cmd, stdout=out, stderr=err)
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)


def run_redirect_command(cmd, out, err=sys.stderr):
    """
    Run operating system command via Python ``subprocess`` and
    redirect output to a file. Uses a pattern suggested by:
    https://www.saltycrane.com/blog/2008/09/how-get-stdout-and-stderr-using-python-subprocess-module/

    :param cmd: Commnand and arguments
    :type cmd: list(str or unicode)
    :param out: Output file name
    :type out: str or unicode
    :param err: Standard error desination (``sys.stderr`` or file)
    :type err: _io.TextIOWrapper
    :raise FileNotFoundError: if the command to run cannot be found
    :raise AssertionError: If the command returns a non-zero exit code
    """
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p_out, p_err = p.communicate()
    exit_code = p.returncode
    with open(out, "wb") as f:
        f.write(p_out)
    # communicate() returns bytes, so convert to string.
    err.write(p_err.decode('utf-8'))
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)


def run_pipe_command(cmd1, cmd2, out=sys.stdout, err=sys.stderr):
    """
    Run operating system command via Python ``subprocess`` and pipe
    output into another command. Uses pattern suggested by:
    https://docs.python.org/2/library/subprocess.html#replacing-shell-pipeline

    :param cmd1: Commnand and arguments
    :type cmd1: list(str or unicode)
    :param cmd2: Commnand and arguments
    :type cmd2: list(str or unicode)
    :param out: Standard output desination (``sys.stdout`` or file)
    :type out: _io.TextIOWrapper
    :param err: Standard error desination (``sys.stderr`` or file)
    :type err: _io.TextIOWrapper
    :raise FileNotFoundError: if the commands to run cannot be found
    :raise AssertionError: If the commands return a non-zero exit code
    """
    process1 = subprocess.Popen(cmd1,
                                stdout=subprocess.PIPE,
                                stderr=err)
    process2 = subprocess.Popen(cmd2,
                                stdin=process1.stdout,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process1.stdout.close()
    p_out, p_err = process2.communicate()
    out.write(p_out.decode('utf-8'))
    err.write(p_err.decode('utf-8'))
    exit_code = process2.returncode
    assert exit_code == 0, ("%s | %s failed with exit code %d"
                            % (cmd1, cmd2, exit_code))


def run_logged_command(cmd,
                       log_file,
                       cmd_file=None,
                       dry_run=False,
                       cmd_to_log=None):
    """
    Run operating system command via Python ``subprocess`` and capture
    standard output and standard error into a log file. Uses
    :py:func:`run_command`.

    If ``cmd_file`` is not ``None`` then the command submitted to the
    operating system is recorded into ``cmd_file``.

    If ``dry_run`` is ``True`` then the command will not be submitted
    to the operating system. Using this with ``cmd_file`` allows a
    record of the command that *would* be submitted to be made.

    For cases where the command to submit may differ depending on
    whether it is run via the command-line or via Python
    ``subprocess``, ``cmd_to_log`` can be used to provide the version
    of the command that needs to be inserted into the ``cmd_file``.

    :param cmd: Commnand and arguments
    :type cmd: list(str or unicode)
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_file: Bash commands file
    :type cmd_file: str or unicode
    :param dry_run: Do not submit command to operating system?
    :type dry_run: bool
    :param cmd_to_log: Command to log
    :type cmd_to_log: list(str or unicode)
    :raise FileNotFoundError: if the command to run cannot be found
    :raise AssertionError: If the command returns a non-zero exit code
    """
    if cmd_file is not None:
        if cmd_to_log is not None:
            cmd_to_log_str = utils.list_to_str(cmd_to_log)
        else:
            cmd_to_log_str = utils.list_to_str(cmd)
        with open(cmd_file, "a") as f:
            f.write(cmd_to_log_str + "\n")
    if dry_run:
        return
    with open(log_file, "a") as f:
        run_command(cmd, f, f)


def run_logged_redirect_command(cmd,
                                out,
                                log_file,
                                cmd_file=None,
                                dry_run=False):
    """
    Run operating system command via Python ``subprocess`` and
    redirect output to a file and capture standard error into a log
    file. Uses :py:func:`run_redirect_command`.

    If ``cmd_file`` is not ``None`` then the command submitted to the
    operating system is recorded into ``cmd_file``.

    If ``dry_run`` is ``True`` then the command will not be submitted
    to the operating system. Using this with ``cmd_file`` allows a
    record of the command that *would* be submitted to be made.

    :param cmd: Commnand and arguments
    :type cmd: list(str or unicode)
    :param out: Output file name
    :type out: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_file: Bash commands file
    :type cmd_file: str or unicode
    :param dry_run: Do not submit command to operating system?
    :type dry_run: bool
    :raise FileNotFoundError: if the command to run cannot be found
    :raise AssertionError: If the command returns a non-zero exit code
    """
    if cmd_file is not None:
        with open(cmd_file, "a") as f:
            f.write(("%s > %s\n" % (utils.list_to_str(cmd), out)))
    if dry_run:
        return
    with open(log_file, "a") as f:
        run_redirect_command(cmd, out, f)


def run_logged_pipe_command(cmd1,
                            cmd2,
                            log_file,
                            cmd_file=None,
                            dry_run=False):
    """
    Run operating system command via Python ``subprocess`` and pipe
    output into another command and capture standard output and
    standard error into a log file. Uses :py:func:`run_pipe_command`.

    If ``cmd_file`` is not ``None`` then the command submitted to the
    operating system is recorded into ``cmd_file``.

    If ``dry_run`` is ``True`` then the command will not be submitted
    to the operating system. Using this with ``cmd_file`` allows a
    record of the command that *would* be submitted to be made.

    :param cmd1: Commnand and arguments
    :type cmd1: list(str or unicode)
    :param cmd2: Commnand and arguments
    :type cmd2: list(str or unicode)
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_file: Bash commands file
    :type cmd_file: str or unicode
    :param dry_run: Do not submit command to operating system?
    :type dry_run: bool
    :raise FileNotFoundError: if the commands to run cannot be found
    :raise AssertionError: If the commands return a non-zero exit code
    """
    if cmd_file is not None:
        with open(cmd_file, "a") as f:
            f.write(("%s | %s\n" % (utils.list_to_str(cmd1),
                                    utils.list_to_str(cmd2))))
    if dry_run:
        return
    with open(log_file, "a") as f:
        run_pipe_command(cmd1, cmd2, f, f)
