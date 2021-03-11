"""
Nextflow-related test fixtures.
"""
import os
import shutil
import tempfile
import subprocess
import pytest
import yaml
import riboviz.process_utils
import riboviz.tools
from riboviz import test
from riboviz import params


def configuration_fixture(request):
    """
    Create a temporary configuration file and temporary directories.

    Load a configuration file, create temporary index, temporary and
    output directories, update the configuration to reference these
    files, save the configuration into a temporary file then yield
    configuration and temporary configuration file.

    On completion delete temporary files and directories.

    Any test module using fixtures that call this function are
    must define a ``TEST_CONFIG_FILE`` parameter specifying the
    configuration file to be copied.

    :param request: Test module using this test fixture
    :type request: _pytest.fixtures.SubRequest
    :return: temporary configuration and configuration file
    :rtype: tuple(dict, str or unicode)
    """
    config_yaml = getattr(request.module, "TEST_CONFIG_FILE")
    with open(config_yaml, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    _, config_file = tempfile.mkstemp(prefix="tmp", suffix=".yaml")
    index_dir = tempfile.mkdtemp("test_index")
    tmp_dir = tempfile.mkdtemp("test_tmp")
    out_dir = tempfile.mkdtemp("test_output")
    config[params.INDEX_DIR] = index_dir
    config[params.TMP_DIR] = tmp_dir
    config[params.OUTPUT_DIR] = out_dir
    with open(config_file, 'w') as f:
        yaml.dump(config, f)

    yield (config, config_file)

    index_dir = config[params.INDEX_DIR]
    tmp_dir = config[params.TMP_DIR]
    out_dir = config[params.OUTPUT_DIR]
    shutil.rmtree(index_dir)
    shutil.rmtree(tmp_dir)
    shutil.rmtree(out_dir)
    if os.path.exists(config_file):
        os.remove(config_file)


configuration = pytest.fixture(scope='function')(configuration_fixture)
""" Function-level fixture for :py:func:`configuration_fixture` """

configuration_module = pytest.fixture(scope='module')(configuration_fixture)
""" Module-level fixture for :py:func:`configuration_fixture` """


def run_nextflow(config_file, envs={}, validate_only=False):
    """
    Run ``nextflow run`` using a given conflguration file.

    :param config_file: Configuration file
    :type config_file: str or unicode
    :param envs: Environment variables
    :type envs: dict
    :param validate_only: Run with \
    :py:const:`riboviz.params.VALIDATE_ONLY` (``--validate_only``)?
    :type validate_only: bool
    :return: exit code
    :rtype: int
    """
    cmd = ["nextflow", "run", test.NEXTFLOW_WORKFLOW,
           "-params-file", config_file, "-ansi-log", "false"]
    if validate_only:
        cmd.append("--{}".format(params.VALIDATE_ONLY))
    env = dict(os.environ)
    for variable, value in envs.items():
        env[variable] = value
    exit_code = subprocess.call(cmd, env=env)
    return exit_code


@pytest.fixture(scope="module")
def nextflow_fixture(configuration_module):
    """
    Run ``nextflow run`` using a given conflguration file.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :raise AssertionError: if the exit code is non-zero
    """
    _, config_file = configuration_module
    exit_code = run_nextflow(config_file)
    assert exit_code == 0, \
        "'nextflow run' returned non-zero exit code %d" % exit_code
