"""
:py:mod:`riboviz.tools.prep_riboviz` demultiplexing tests.

The test suite runs :py:mod:`riboviz.tools.prep_riboviz` using
``vignette/simdata_multiplex_config.yaml``
(:py:const:`riboviz.test.SIMDATA_MULTIPLEX_CONFIG`) and simulated data
in ``data/simdata/`` (created by
:py:mod:`riboviz.tools.create_fastq_simdata`).

It then validates the outputs of adaptor trimming, UMI extraction,
demultiplexing and deduplication steps against the expected outputs,
also in ``data/simdata/``. Collated TPMs are also validated.

Each test function is configured with the module-level fixture
:py:func:`riboviz.test.tools.prep_riboviz_fixture` to ensure
that :py:mod:`riboviz.tools.prep_riboviz` is run once before the
test functions are run.
"""
import os
import pytest
import riboviz
import riboviz.test
from riboviz import demultiplex_fastq
from riboviz import fastq
from riboviz import params
from riboviz import utils
from riboviz import workflow_files
from riboviz.test.tools import configuration_module
from riboviz.test.tools import prep_riboviz_fixture
from riboviz.test.tools.test_prep_riboviz_simdata_umi \
    import check_umi_groups
from riboviz.test.tools.test_prep_riboviz_simdata_umi \
    import check_tpms_collated_tsv


TEST_CONFIG_FILE = riboviz.test.SIMDATA_MULTIPLEX_CONFIG
"""
Test file location constant, used by a callback in
:py:func:`riboviz.test.tools.configuration_module`.
"""


@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_adaptor_trimming(configuration_module):
    """
    Test that the results of adaptor trimming are as expected.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR,
        fastq.FASTQ_FORMAT.format("multiplex_umi_barcode"))
    actual_output = os.path.join(
        config[params.TMP_DIR],
        workflow_files.ADAPTER_TRIM_FQ_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    fastq.equal_fastq(expected_output, actual_output)


@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_barcode_umi_extract(configuration_module):
    """
    Test that the results of barcode and UMI extraction are as expected.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR,
        fastq.FASTQ_FORMAT.format("multiplex"))
    actual_output = os.path.join(
        config[params.TMP_DIR],
        workflow_files.UMI_EXTRACT_FQ_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    fastq.equal_fastq(expected_output, actual_output)


@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_deplex_num_reads(configuration_module):
    """
    Test that the number of reads summary, produced during
    demultiplexing, is as expected.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    actual_dir = os.path.join(
        config[params.TMP_DIR],
        workflow_files.DEPLEX_DIR_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    actual_output = os.path.join(actual_dir,
                                 demultiplex_fastq.NUM_READS_FILE)
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR,
        "deplex",
        demultiplex_fastq.NUM_READS_FILE)
    utils.equal_tsv(expected_output, actual_output, na_to_empty_str=True)


@pytest.mark.parametrize(
    "sample_id", ["Tag0", "Tag1", "Tag2", "Unassigned"])
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_deplex_reads(configuration_module, sample_id):
    """
    Test that the FASTQ files output by demultiplexing are as
    expected.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID for demultiplexed reads
    :type sample_id: str or unicode
    """
    # Actual data has a .fq extension.
    actual_file_name = fastq.FQ_FORMAT.format(sample_id)
    # Simulated data has a .fastq extension.
    expected_file_name = fastq.FASTQ_FORMAT.format(sample_id)
    config, _ = configuration_module
    actual_dir = os.path.join(
        config[params.TMP_DIR],
        workflow_files.DEPLEX_DIR_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    actual_output = os.path.join(actual_dir, actual_file_name)
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR, "deplex", expected_file_name)
    fastq.equal_fastq(expected_output, actual_output)


@pytest.mark.parametrize("sample_id", ["Tag3"])
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_deplex_unmatched_samples(configuration_module, sample_id):
    """
    Test that no FASTQ files are output by demultiplexing for
    samples for which there were no matching reads.

    The definition of the simulated data means that Tag3 has no
    matches, as Tag0|1|2 will match any barcodes first. Check
    there is no Tag3-related output file.e

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID for demultiplexed reads
    :type sample_id: str or unicode
    """
    # Actual data has a .fq extension.
    config, _ = configuration_module
    output_dir = os.path.join(
        config[params.TMP_DIR],
        workflow_files.DEPLEX_DIR_FORMAT.format(
            "multiplex_umi_barcode_adaptor"))
    file_name = os.path.join(output_dir, fastq.FQ_FORMAT.format(sample_id))
    assert not os.path.exists(file_name)


@pytest.mark.parametrize("sample_id", ["Tag0", "Tag1", "Tag2"])
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_deplex_umi_groups(configuration_module, sample_id):
    """
    Test that the UMI groups are as expected. See
    :py:func:`riboviz.test.tools.test_prep_riboviz_simdata_umi.check_umi_groups`.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID for demultiplexed reads
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    check_umi_groups(config, sample_id, 5)


@pytest.mark.parametrize("sample_id", ["Tag0", "Tag1", "Tag2"])
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_deplex_tpms_collated_tsv(configuration_module, sample_id):
    """
    Test that the collated TPMs are as expected. See
    :py:func:`riboviz.test.tools.test_prep_riboviz_simdata_umi.check_tpms_collated_tsv`.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID for demultiplexed reads
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    check_tpms_collated_tsv(config, sample_id, 4)
