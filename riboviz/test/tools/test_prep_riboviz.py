"""
:py:mod:`riboviz.tools.prep_riboviz` error handling and exit code tests.

The test suite runs :py:mod:`riboviz.tools.prep_riboviz` in "dry-run"
mode using ``vignette/vignette_config.yaml``
(:py:const:`riboviz.test.VIGNETTE_CONFIG`).
"""
import yaml
import pytest
import riboviz
import riboviz.process_utils
import riboviz.test
import riboviz.tools
from riboviz import params
from riboviz.tools import prep_riboviz
from riboviz.test.tools import configuration  # Test fixture

TEST_CONFIG_FILE = riboviz.test.VIGNETTE_CONFIG
"""
Test file location constant, used by a callback in
:py:func:`riboviz.test.tools.configuration_module`.
"""


def test_missing_config_file():
    """
    Test that using a non-existent configuration file gives the
    expected error code.
    """
    exit_code = prep_riboviz.prep_riboviz(riboviz.R_SCRIPTS,
                                          "nosuch.yaml",
                                          True)
    assert exit_code == prep_riboviz.EXIT_FILE_NOT_FOUND_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_no_fq_files_error(configuration):
    """
    Test that specifying no sample files gives the expected error
    code.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    config[params.FQ_FILES] = []
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_fq_files_multiplex_fq_files_error(configuration):
    """
    Test that specifying both sample files and multiplexed sample
    files gives the expected error code.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    config[params.MULTIPLEX_FQ_FILES] = ["somefile.fastq"]
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_multiplex_fq_files_missing_sample_sheet_error(configuration):
    """
    Test that specifying multiplexed sample files and a non-existent
    sample sheet gives the expected error code.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    del config[params.FQ_FILES]
    config[params.MULTIPLEX_FQ_FILES] = ["somefile.fastq"]
    config[params.SAMPLE_SHEET] = "noSuchFile.tsv"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_FILE_NOT_FOUND_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_missing_fq_files(configuration):
    """
    Test that if sample files are non-existent this gives the expected
    error code.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    config[params.FQ_FILES] = {
        "WT3AT": "nosuch.fastq.gz",
        "WTnone": "nosuch.fastq.gz"
    }
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_PROCESSING_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


def test_missing_dir_in(configuration):
    """
    Test that if the input directory is non-existent this gives the
    expected error code.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, config_path = configuration
    del config[params.INPUT_DIR]
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_CONFIG_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter", [params.RRNA_FASTA_FILE,
                                       params.ORF_FASTA_FILE,
                                       params.ORF_GFF_FILE,
                                       params.FEATURES_FILE,
                                       params.T_RNA_FILE,
                                       params.CODON_POSITIONS_FILE])
def test_missing_files_error(configuration, parameter):
    """
    Test that using non-existent error code for other input files
    gives the expected error code.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param parameter: file name configuration parameter
    :type parameter: str or unicode
    """
    config, config_path = configuration
    config[parameter] = "noSuchFile.txt"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    exit_code = prep_riboviz.prep_riboviz(riboviz.R_SCRIPTS,
                                          config_path,
                                          True)
    assert exit_code == prep_riboviz.EXIT_FILE_NOT_FOUND_ERROR, \
        "prep_riboviz returned with unexpected exit code %d" % exit_code
