"""
subprocess-related utilities.
"""
import subprocess
import sys
from riboviz import utils


def run_command(cmd, out=sys.stdout, err=sys.stderr):
    """
    Helper function to run shell command.

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param out: Standard output desination (stdout or file)
    :type out: _io.TextIOWrapper
    :param err: Standard error desination (stderr or file)
    :type err: _io.TextIOWrapper
    :raise AssertionError: if the command returns a non-zero exit code
    """
    exit_code = subprocess.call(cmd, stdout=out, stderr=err)
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)


def run_redirect_command(cmd, out, err=sys.stderr):
    """
    Helper function to run shell command and redirect output to a file.

    Uses pattern suggested by:
    https://www.saltycrane.com/blog/2008/09/how-get-stdout-and-stderr-using-python-subprocess-module/

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param out: Output file name
    :type out: str or unicode
    :param err: Standard error desination (stderr or file)
    :type err: _io.TextIOWrapper
    :raise OSError: if the command being run cannot be found
    (Python 2)
    :raise FileNotFoundError: if the command being run cannot be found
    (Python 3)
    :raise AssertionError: if the command returns a non-zero exit code
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
    Helper function to run shell command and pipe output into another.

    Uses pattern suggested by:
    https://docs.python.org/2/library/subprocess.html#replacing-shell-pipeline

    :param cmd1: Commnand to run and its arguments
    :type cmd1: list(str or unicode)
    :param cmd2: Commnand to run and its arguments
    :type cmd2: list(str or unicode)
    :param out: Standard output desination (stdout or file)
    :type out: _io.TextIOWrapper
    :param err: Standard error desination (stderr or file)
    :type err: _io.TextIOWrapper
    :raise OSError: if the command being run cannot be found
    (Python 2)
    :raise FileNotFoundError: if the command being run cannot be found
    (Python 3)
    :raise AssertionError: if the commands returns a non-zero exit code
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


def run_logged_command(cmd, log_file, cmd_file=None, dry_run=False):
    """
    Helper function to run shell command and capture stdout and stderr
    in a log file. The command submitted to the shell is also
    recorded, if cmd_file is not None.

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param log_file: File to log stdout and stderr.
    :type log_file: str or unicode
    :param cmd_file: File to log command to
    :type cmd_file: str or unicode
    :param dry_run: Do not submit command to shell - use with cmd_file
    to log commands that would be run
    :type dry_run: bool
    :raise OSError: if the command being run cannot be found
    (Python 2)
    :raise FileNotFoundError: if the command being run cannot be found
    (Python 3)
    :raise AssertionError: if the command returns a non-zero exit code
    """
    if cmd_file is not None:
        with open(cmd_file, "a") as f:
            f.write(utils.list_to_str(cmd) + "\n")
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
    Helper function to run shell command and redirect output to a file
    and capture stderr in a log file. The command submitted to the
    shell is also recorded, if cmd_file is not None.

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param out: Output file name
    :type out: str or unicode
    :param log_file: File to log stdout and stderr.
    :type log_file: str or unicode
    :param cmd_file: File to log command to
    :type cmd_file: str or unicode
    :param dry_run: Do not submit command to shell - use with cmd_file
    to log commands that would be run
    :type dry_run: bool
    :raise OSError: if the command being run cannot be found
    (Python 2)
    :raise FileNotFoundError: if the command being run cannot be found
    (Python 3)
    :raise AssertionError: if the command returns a non-zero exit code
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
    Helper function to run shell command and pipe output into another
    and capture stdout and stderr in a log file. The command submitted
    to the shell is also recorded, if cmd_file is not None.

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param log_file: File to log stdout and stderr.
    :type log_file: str or unicode
    :param cmd_file: File to log command to
    :type cmd_file: str or unicode
    :param dry_run: Do not submit command to shell - use with cmd_file
    to log commands that would be run
    :type dry_run: bool
    :raise OSError: if the command being run cannot be found
    (Python 2)
    :raise FileNotFoundError: if the command being run cannot be found
    (Python 3)
    :raise AssertionError: if the commands returns a non-zero exit code
    """
    if cmd_file is not None:
        with open(cmd_file, "a") as f:
            f.write(("%s | %s\n" % (utils.list_to_str(cmd1),
                                    utils.list_to_str(cmd2))))
    if dry_run:
        return
    with open(log_file, "a") as f:
        run_pipe_command(cmd1, cmd2, f, f)
