"""
prep_riboviz.py test suite to test error handling and exit
codes. prep_riboviz.py is run in "dry-run" mode.
"""
import yaml
import pytest
import riboviz
import riboviz.test
import riboviz.process_utils
import riboviz.validation
import riboviz.tools
from riboviz import params
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
                                          "nosuch.yaml",
                                          True)
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


@pytest.mark.parametrize("index", [params.R_RNA_FASTA_FILE,
                                   params.ORF_FASTA_FILE])
def test_index_error_missing_index_files(configuration, index):
    """
    Test that the rRNA_fasta and orf_fasta configuration value being
    non-existent files causes EXIT_INDEX_ERROR to be returned.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    :param index: index file name configuration parameter
    :type index: str or unicode
    """
    config, config_path = configuration
    config[index] = "nosuch.fa"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_INDEX_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_no_data_error(configuration):
    """
    Test that no samples being specified causes
    EXIT_NO_DATA_ERROR to be returned.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    config[params.FQ_FILES] = []
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_NO_DATA_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_samples_error_missing_samples(configuration):
    """
    Test that if all samples are non-existent files then
    EXIT_DATA_ERROR is returned.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode
    """
    config, config_path = configuration
    config[params.FQ_FILES] = {
        "WT3AT": "nosuch.fastq.gz",
        "WTnone": "nosuch.fastq.gz"
    }
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_DATA_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


@pytest.mark.parametrize("file_config", [params.ORF_GFF_FILE,
                                         params.FEATURES_FILE,
                                         params.T_RNA,
                                         params.CODON_POS])
def test_missing_files_error(configuration, file_config):
    """
    Test that non-existent files being specified for org_gff_file,
    features_file, t_rna and codon_pos then EXIT_DATA_ERROR is
    returned.

    :param configuration: configuration and path to configuration file
    (pytest fixture defined in conftest.py)
    :type configuration: tuple(dict, str or unicode)
    :param file_config: file name configuration parameter
    :type file_config: str or unicode
    """
    config, config_path = configuration
    config[file_config] = "noSuchFile.txt"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_DATA_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_config_error_missing_dir_in(configuration):
    """
    Test that a missing "dir_in" configuration value causes
    EXIT_CONFIG_ERROR to be returned.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    del config[params.INPUT_DIR]
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code
