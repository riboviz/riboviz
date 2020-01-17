"""
riboviz.tools.prep_riboviz test suite to test adaptor trimming,
barcode and UMI extraction, demultpliexing and deduplication.

The test suite runs riboviz.tools.prep_riboviz using a copy of
"vignette/simdata_multiplex_config.yaml" and the simulated data in
"data/simdata/". It then validates the outputs of the adaptor
trimming, barcode and UMI extraction, demultiplexing and deduplication
steps against the expected outputs, also in "data/simdata/".

The simulated data in "data/simdata/" is expected to have been created
using riboviz.tools.create_fastq_simdata.
"""
import os
import pytest
import riboviz
import riboviz.process_utils
import riboviz.test
import riboviz.tools
import riboviz.validation
from riboviz import demultiplex_fastq
from riboviz import params
from riboviz.tools import prep_riboviz
from riboviz.test.tools import configuration_module  # Test fixture
from riboviz.test.tools import run_prep_riboviz  # Test fixture
from riboviz.test.tools.test_prep_riboviz_simdata_umi import check_umi_groups
from riboviz.test.tools.test_prep_riboviz_simdata_umi import check_tpms_collated_tsv


TEST_CONFIG_FILE = riboviz.test.SIMDATA_MULTIPLEX_CONFIG
"""
YAML configuration used as a template configuration by these tests -
required by configuration test fixture
"""


@pytest.mark.usefixtures("run_prep_riboviz")
def test_adaptor_trimming(configuration_module):
    """
    Validate that adaptor trimming, performed by "cutadapt" produces
    the expected results.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR,
        "multiplex_umi_barcode.fastq")
    actual_output = os.path.join(
        config[params.TMP_DIR],
        prep_riboviz.TMP_ADAPTER_TRIM_FQ_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    riboviz.validation.equal_fastq(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
def test_barcode_umi_extract(configuration_module):
    """
    Validate that barcode and UMI extraction, performed by "umi_tools
    extract" produces the expected results.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR,
        "multiplex.fastq")
    actual_output = os.path.join(
        config[params.TMP_DIR],
        prep_riboviz.TMP_UMI_EXTRACT_FQ_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    riboviz.validation.equal_fastq(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
def test_deplex_num_reads(configuration_module):
    """
    Validate that "num_reads.tsv", produced by
    riboviz.demultiplex_fastq has the expected content.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    actual_dir = os.path.join(
        config[params.TMP_DIR],
        prep_riboviz.TMP_DEPLEX_DIR_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    actual_output = os.path.join(actual_dir,
                                 demultiplex_fastq.NUM_READS_FILE)
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR,
        "deplex",
        demultiplex_fastq.NUM_READS_FILE)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.parametrize(
    "fastq", ["Tag0.fastq", "Tag1.fastq", "Tag2.fastq", "Unassigned.fastq"])
@pytest.mark.usefixtures("run_prep_riboviz")
def test_deplex_reads(configuration_module, fastq):
    """
    Validate that ".fastq", produced by
    riboviz.demultiplex_fastq have the expected content.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    :param fastq: FASTQ file
    :type fastq: str or unicode
    """
    config, _ = configuration_module
    actual_dir = os.path.join(
        config[params.TMP_DIR],
        prep_riboviz.TMP_DEPLEX_DIR_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    actual_output = os.path.join(actual_dir, fastq)
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR, "deplex", fastq)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.parametrize("sample_id", ["Tag0", "Tag1", "Tag2"])
@pytest.mark.usefixtures("run_prep_riboviz")
def test_deplex_umi_groups(configuration_module, sample_id):
    """
    Validate the information on UMI groups post-"umi_tools extract",
    for each demultiplexed file, by parsing the ".tsv" file output by
    "umi_tools group".

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID for demultiplexed reads
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    check_umi_groups(config, sample_id, 5)


@pytest.mark.parametrize("sample_id", ["Tag0", "Tag1", "Tag2"])
@pytest.mark.usefixtures("run_prep_riboviz")
def test_deplex_tpms_collated_tsv(configuration_module, sample_id):
    """
    Validate the "TPMs_collated.tsv" file produced by the workflow.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID for demultiplexed reads
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    check_tpms_collated_tsv(config, sample_id, 4)
