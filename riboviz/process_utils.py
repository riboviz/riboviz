"""
subprocess-related utilities.
"""

import subprocess


def run_command(cmd):
    """
    Helper function to run shell command.

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :raise FileNotFoundError: if the command being run cannot be found
    :raise AssertionError: if the command returns a non-zero exit code
    """
    exit_code = subprocess.call(cmd)
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)


def run_redirect_command(cmd, out_file):
    """
    Helper function to run shell command and redirect output to a file.

    Use pattern suggested by:
    https://www.saltycrane.com/blog/2008/09/how-get-stdout-and-stderr-using-python-subprocess-module/

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param out_file: Output file
    :type out_file: str or unicode
    :raise FileNotFoundError: if the command being run cannot be found
    :raise AssertionError: if the command returns a non-zero exit code
    """
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdout, _ = p.communicate()
    exit_code = p.returncode
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)
    with open(out_file, "wb") as f:
        f.write(stdout)


def run_pipe_command(cmd1, cmd2):
    """
    Helper function to run shell command and pipe output into another.

    Use pattern suggested by:
    https://docs.python.org/2/library/subprocess.html#replacing-shell-pipeline
    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :raise FileNotFoundError: if the commands being run cannot be found
    :raise AssertionError: if the commands returns a non-zero exit code
    """
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    process2 = subprocess.Popen(cmd2, stdin=process1.stdout)
    process1.stdout.close()
    _, _ = process2.communicate()
    exit_code = process2.returncode
    assert exit_code == 0, ("%s | %s failed with exit code %d"
                            % (cmd1, cmd2, exit_code))
