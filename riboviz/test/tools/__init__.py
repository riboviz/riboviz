"""
prep_riboviz.py-related test fixtures.
"""
import os
import shutil
import tempfile
import pytest
import yaml
import riboviz
import riboviz.test
import riboviz.process_utils
import riboviz.validation
import riboviz.tools
from riboviz import params
from riboviz.tools import prep_riboviz


def configuration_fixture(request):
    """
    Create a temporary configuration file and temporary directories.

    Load a configuration file, create temporary index, tmp, and output
    directories, update configuration to reference these files, save
    as a temporary file, yield configuration to callers and delete all
    temporary files on test completion.

    The test module using fixtures that call this function are
    expected to define a TEST_CONFIG_FILE parameter specifying the
    configuration file to be copied.
    config_yaml = getattr(request.module, "TEST_CONFIG_FILE")

    :param request: Test module using this test fixture
    :type request: _pytest.fixtures.SubRequest
    :return: configuration, path to configuration file
    :rtype: tuple(dict, str or unicode)
    """
    config_yaml = getattr(request.module, "TEST_CONFIG_FILE")
    with open(config_yaml, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    _, config_path = tempfile.mkstemp(prefix="tmp", suffix=".yaml")
    index_dir = tempfile.mkdtemp("tmp_riboviz_test_tools_index")
    tmp_dir = tempfile.mkdtemp("tmp_riboviz_test_tools_tmp")
    out_dir = tempfile.mkdtemp("tmp_riboviz_test_tools_out")
    logs_dir = tempfile.mkdtemp("tmp_riboviz_test_tools_logs")
    _, cmd_file = tempfile.mkstemp(prefix="tmp", suffix=".sh")
    config[params.INDEX_DIR] = index_dir
    config[params.TMP_DIR] = tmp_dir
    config[params.OUTPUT_DIR] = out_dir
    config[params.LOGS_DIR] = logs_dir
    config[params.CMD_FILE] = cmd_file
    with open(config_path, 'w') as f:
        yaml.dump(config, f)

    yield (config, config_path)

    index_dir = config[params.INDEX_DIR]
    tmp_dir = config[params.TMP_DIR]
    out_dir = config[params.OUTPUT_DIR]
    logs_dir = config[params.LOGS_DIR]
    cmd_file = config[params.CMD_FILE]
    shutil.rmtree(index_dir)
    shutil.rmtree(tmp_dir)
    shutil.rmtree(out_dir)
    shutil.rmtree(logs_dir)
    if os.path.exists(cmd_file):
        os.remove(cmd_file)
    if os.path.exists(config_path):
        os.remove(config_path)


configuration = pytest.fixture(scope='function')(configuration_fixture)
""" Function-level fixture for configuration_fixture """

configuration_module = pytest.fixture(scope='module')(configuration_fixture)
""" Module-level fixture for configuration_fixture """
