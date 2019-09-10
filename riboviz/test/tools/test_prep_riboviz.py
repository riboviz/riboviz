"""
prep_riboviz.py test suite to test error handling and exit codes.
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
from riboviz.tools import prep_riboviz


@pytest.fixture(scope="module")
def arguments():
    """
    Fixture to create default locations of prep_riboviz.py
    arguments.

    :return: Python scripts directory, R scripts directory,
    data directory, vignette configuration file
    :rtype: tuple(str or unicode, str or unicode,
    str or unicode, str or unicode)
    """
    arguments = (riboviz.test.PY_SCRIPTS,
                 riboviz.test.R_SCRIPTS,
                 riboviz.test.DATA_DIR,
                 riboviz.test.VIGNETTE_CONFIG)
    yield arguments


@pytest.fixture(scope="function")
def configuration(arguments):
    """
    Create a temporary copy of the vignette configuration file,
    create temporary index, tmp, and output directories, save
    updated configuration and delete all after tests complete.

    :param arguments: Python scripts directory, R scripts directory,
    data directory, vignette configuration file
    (pytest fixture defined in this module)
    :type arguments: tuple(str or unicode, str or unicode,
    str or unicode, str or unicode)
    :return: configuration, path to configuration file
    :rtype: tuple(dict, str or unicode)
    """
    _, _, _, config_yaml = arguments
    with open(config_yaml, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)

    _, path = tempfile.mkstemp(prefix="tmp", suffix=".yaml")
    index_dir = tempfile.mkdtemp("tmp_test_prep_riboviz_index")
    tmp_dir = tempfile.mkdtemp("tmp_test_prep_riboviz_tmp")
    out_dir = tempfile.mkdtemp("tmp_test_prep_riboviz_out")
    logs_dir = tempfile.mkdtemp("tmp_test_prep_riboviz_logs")

    config["dir_index"] = index_dir
    config["dir_tmp"] = tmp_dir
    config["dir_out"] = out_dir
    config["dir_logs"] = logs_dir

    with open(path, 'w') as f:
        yaml.dump(config, f)
    yield (config, path)
    if os.path.exists(path):
        os.remove(path)
    shutil.rmtree(index_dir)
    shutil.rmtree(tmp_dir)
    shutil.rmtree(out_dir)
    shutil.rmtree(logs_dir)


def test_config_error_missing_config_file(arguments):
    """
    Test that a non-existent configuration file causes
    EXIT_CONFIG_ERROR to be returned.

    :param arguments: Python scripts directory, R scripts directory,
    data directory, vignette configuration file
    (pytest fixture defined in this module)
    :type arguments: tuple(str or unicode, str or unicode,
    str or unicode, str or unicode)
    """
    py_scripts, r_scripts, data_dir, _ = arguments
    exit_code = prep_riboviz.prep_riboviz(py_scripts,
                                          r_scripts,
                                          data_dir,
                                          "nosuch.yaml")
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_index_error_missing_fa(arguments, configuration):
    """
    Test that the rRNA_fasta configuration value being a non-existent
    file causes EXIT_INDEX_ERROR to be returned.

    :param arguments: Python scripts directory, R scripts directory,
    data directory, vignette configuration file
    (pytest fixture defined in this module)
    :type arguments: tuple(str or unicode, str or unicode,
    str or unicode, str or unicode)
    :param configuration: configuration and path to configuration file
    (pytest fixture defined in conftest.py)
    :type configuration: tuple(dict, str or unicode)
    """
    py_scripts, r_scripts, data_dir, _ = arguments
    config, path = configuration
    config["rRNA_fasta"] = "nosuch.fa"
    with open(path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(py_scripts,
                                          r_scripts,
                                          data_dir,
                                          path)
    assert exit_code == prep_riboviz.EXIT_INDEX_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_no_samples_error(arguments, configuration):
    """
    Test that no samples being specified causes
    EXIT_NO_SAMPLES_ERROR to be returned. Indexing is skipped to save
    time.

    :param arguments: Python scripts directory, R scripts directory,
    data directory, vignette configuration file
    (pytest fixture defined in this module)
    :type arguments: tuple(str or unicode, str or unicode,
    str or unicode, str or unicode)
    :param configuration: configuration and path to configuration file
    (pytest fixture defined in conftest.py)
    :type configuration: tuple(dict, str or unicode)
    """
    py_scripts, r_scripts, data_dir, _ = arguments
    config, path = configuration
    config["build_indices"] = False
    config["fq_files"] = []
    with open(path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(py_scripts,
                                          r_scripts,
                                          data_dir,
                                          path)
    assert exit_code == prep_riboviz.EXIT_NO_SAMPLES_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_samples_error_missing_samples(arguments, configuration):
    """
    Test that if all samples are non-existet files then
    EXIT_SAMPLES_ERROR is returned. Indexing is skipped to save time.

    :param arguments: Python scripts directory, R scripts directory,
    data directory, vignette configuration file
    (pytest fixture defined in this module)
    :type arguments: tuple(str or unicode, str or unicode,
    str or unicode, str or unicode)
    :param configuration: configuration and path to configuration file
    (pytest fixture defined in conftest.py)
    :type configuration: tuple(dict, str or unicode)
    """
    py_scripts, r_scripts, data_dir, _ = arguments
    config, path = configuration
    config["build_indices"] = False
    config["fq_files"] = {
        "WT3AT": "nosuch.fastq.gz",
        "WTnone": "nosuch.fastq.gz"
    }
    with open(path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(py_scripts,
                                          r_scripts,
                                          data_dir,
                                          path)
    assert exit_code == prep_riboviz.EXIT_SAMPLES_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_config_error_missing_dir_in(arguments, configuration):
    """
    Test that a missing "dir_in" configuration value causes
    EXIT_CONFIG_ERROR to be returned. Indexing is skipped to save
    time.

    :param arguments: Python scripts directory, R scripts directory,
    data directory, vignette configuration file
    (pytest fixture defined in this module)
    :type arguments: tuple(str or unicode, str or unicode,
    str or unicode, str or unicode)
    :param configuration: configuration and path to configuration file
    (pytest fixture defined in conftest.py)
    :type configuration: tuple(dict, str or unicode)
    """
    py_scripts, r_scripts, data_dir, _ = arguments
    config, path = configuration
    config["build_indices"] = False
    del config["dir_in"]
    with open(path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(py_scripts,
                                          r_scripts,
                                          data_dir,
                                          path)
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code
