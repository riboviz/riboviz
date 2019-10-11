"""
prep_riboviz.py test suite to test adatpor trimming, UMI extraction
and deduplication.

The test suite runs prep_riboviz.py using a copy of
`vignette/vignette_simdata_config.yaml` and the simulated data in
`data/`. It then validates the outputs of the UMI-tools-specific
phases against the expected outputs, also in `data/`.
"""
import os
import pytest
import riboviz
import riboviz.process_utils
import riboviz.test
import riboviz.tools
import riboviz.validation
from riboviz.tools import prep_riboviz
from riboviz.test.tools import configuration_module  # Test fixture


TEST_CONFIG_FILE = riboviz.test.SIMDATA_CONFIG
"""
YAML configuration used as a template configuration by these tests -
required by configuration test fixture
"""


@pytest.fixture(scope="module")
def run_prep_riboviz(configuration_module):
    """
    Fixture to run prep_riboviz.py.

    :param configuration_module: configuration and path to
    configuration file  (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    _, config_path = configuration_module
    exit_code = prep_riboviz.prep_riboviz(riboviz.test.PY_SCRIPTS,
                                          riboviz.test.R_SCRIPTS,
                                          config_path)
    assert exit_code == 0, \
        "prep_riboviz returned non-zero exit code %d" % exit_code


@pytest.mark.usefixtures("run_prep_riboviz")
def test_adaptor_trimming(configuration_module):
    """
    Validate that adaptor trimming, performed by `cutadapt` produces
    the expected results, by comparing the FASTQ file output to a
    pre-calculated one in `data/`.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    expected_output = os.path.join(riboviz.test.DATA_DIR,
                                   "simdata_UMI5and3_4nt.fastq")
    actual_output = os.path.join(config["dir_tmp"],
                                 "simdata5and3_trim.fq")
    riboviz.validation.equal_fastq(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
def test_umi_extract(configuration_module):
    """
    Validate that UMI extraction, performed by `umi_tools extract`
    produces the expected results, by comparing the FASTQ file output
    to a pre-calculated one in `data/`.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    expected_output = os.path.join(riboviz.test.DATA_DIR,
                                   "simdata_extracted_UMI5and3_4nt.fastq")
    actual_output = os.path.join(config["dir_tmp"],
                                 "simdata5and3_extract_trim.fq")
    riboviz.validation.equal_fastq(expected_output, actual_output)
