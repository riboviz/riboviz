"""
prep_riboviz.py test suite to test error handling and exit codes.
"""
import yaml
import riboviz
import riboviz.test
import riboviz.process_utils
import riboviz.validation
import riboviz.tools
from riboviz.tools import prep_riboviz
from riboviz.test.tools import configuration  # Test fixture


TEST_CONFIG_FILE = riboviz.test.VIGNETTE_CONFIG
"""
YAML configuration used as a template configuration by these tests -
required by configuration test fixture
"""


def test_config_error_missing_config_file():
    """
    Test that a non-existent configuration file causes
    EXIT_CONFIG_ERROR to be returned.
    """
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          "nosuch.yaml")
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_index_error_missing_fa(configuration):
    """
    Test that the rRNA_fasta configuration value being a non-existent
    file causes EXIT_INDEX_ERROR to be returned.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    config["rRNA_fasta"] = "nosuch.fa"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path)
    assert exit_code == prep_riboviz.EXIT_INDEX_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_no_samples_error(configuration):
    """
    Test that no samples being specified causes
    EXIT_NO_SAMPLES_ERROR to be returned. Indexing is skipped to save
    time.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    config["build_indices"] = False
    config["fq_files"] = []
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path)
    assert exit_code == prep_riboviz.EXIT_NO_SAMPLES_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_samples_error_missing_samples(configuration):
    """
    Test that if all samples are non-existet files then
    EXIT_SAMPLES_ERROR is returned. Indexing is skipped to save time.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    config["build_indices"] = False
    config["fq_files"] = {
        "WT3AT": "nosuch.fastq.gz",
        "WTnone": "nosuch.fastq.gz"
    }
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path)
    assert exit_code == prep_riboviz.EXIT_SAMPLES_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_config_error_missing_dir_in(configuration):
    """
    Test that a missing "dir_in" configuration value causes
    EXIT_CONFIG_ERROR to be returned. Indexing is skipped to save
    time.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    config["build_indices"] = False
    del config["dir_in"]
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path)
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code
