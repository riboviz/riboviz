"""
Test Nextflow workflow using configuration with UMI extraction and
deduplication.

The test suite runs ``nextflow run prep_riboviz.nf`` using
``vignette/simdata_umi_config.yaml``
(:py:const:`riboviz.test.SIMDATA_UMI_CONFIG`) and simulated data in
``data/simdata/`` (created by
:py:mod:`riboviz.tools.create_fastq_simdata`).

It then validates the outputs of adaptor trimming, UMI extraction,
and deduplication against the expected outputs, also in
``data/simdata/``. Collated TPMs are also validated.

Each test function is configured with the module-level fixture
:py:func:`riboviz.test.nextflow.nextflow_fixture` to ensure
that ``nextflow run prep_riboviz.nf`` is run once before the
test functions are run.
"""
import os
import pytest
import pandas as pd
import riboviz
import riboviz.test
from riboviz import fastq
from riboviz import params
from riboviz import umi_tools
from riboviz import workflow_files
from riboviz import workflow_r
from riboviz.test.nextflow import configuration_module  # pylint: disable=unused-import
from riboviz.test.nextflow import nextflow_fixture  # pylint: disable=unused-import


TEST_CONFIG_FILE = riboviz.test.SIMDATA_UMI_CONFIG
"""
Test file location constant, used by a callback in
:py:func:`riboviz.test.nextflow.configuration_module`.
"""


@pytest.mark.parametrize("sample_id", [riboviz.test.SIMDATA_UMI_SAMPLE])
@pytest.mark.usefixtures("nextflow_fixture")
def test_adaptor_trimming(configuration_module, sample_id):
    """
    Test that the results of adaptor trimming are as expected.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR,
        fastq.FASTQ_FORMAT.format(sample_id + "_umi"))
    actual_output = os.path.join(
        config[params.TMP_DIR],
        sample_id,
        workflow_files.ADAPTER_TRIM_FQ)
    fastq.equal_fastq(expected_output, actual_output)


@pytest.mark.parametrize("sample_id", [riboviz.test.SIMDATA_UMI_SAMPLE])
@pytest.mark.usefixtures("nextflow_fixture")
def test_umi_extract(configuration_module, sample_id):
    """
    Test that the results of UMI extraction are as expected.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    expected_output = os.path.join(
        riboviz.test.SIMDATA_DIR,
        fastq.FASTQ_FORMAT.format(sample_id))
    actual_output = os.path.join(
        config[params.TMP_DIR],
        sample_id,
        workflow_files.UMI_EXTRACT_FQ)
    fastq.equal_fastq(expected_output, actual_output)


def check_umi_groups(tmp_dir, sample_id, num_groups):
    """
    Test that the UMI groups are as expected.

    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    :param sample_id: sample ID
    :type sample_id: str or unicode
    :param num_groups: Expected number of groups
    :type num_groups: int
    """
    groups_tsv = os.path.join(tmp_dir,
                              sample_id,
                              workflow_files.POST_DEDUP_GROUPS_TSV)
    groups = pd.read_csv(groups_tsv, sep="\t")
    assert groups.shape[0] == num_groups, \
        ("Expected %d unique groups but found %d"
         % (num_groups, groups.shape[0]))
    assert (groups[umi_tools.UMI_COUNT] == 1).all(), \
        "Expected each umi_count to be 1"
    assert (groups[umi_tools.FINAL_UMI_COUNT] == 1).all(), \
        "Expected each final_umi_count to be 1"
    # Check group IDs are unique when compared to 1,...,number of
    # groups.
    group_ids = list(groups[umi_tools.UNIQUE_ID])
    group_ids.sort()
    expected_group_ids = list(range(num_groups))
    assert expected_group_ids == group_ids, \
        ("Expected group_ids %s but found %s" % (str(expected_group_ids),
                                                 str(group_ids)))
    # Check each representative read does indeed come from a unique
    # UMI group by parsing the read ID.
    # riboviz.create_fastq_simdata creates read IDs of form:
    # "EWSim-<GROUP>.<MEMBER>-umi<5PRIME>-read<READ>-umi<3PRIME>"
    # where <GROUP> is 1-indexed.
    groups_from_read_ids = [
        int(read_id.split("-")[1].split(".")[0]) - 1
        for read_id in groups[umi_tools.READ_ID]
    ]
    groups_from_read_ids.sort()
    assert groups_from_read_ids == group_ids, \
        ("Reads in read_ids %s are not from unique groups" %
         (str(list(groups[umi_tools.READ_ID]))))


@pytest.mark.parametrize("sample_id", [riboviz.test.SIMDATA_UMI_SAMPLE])
@pytest.mark.usefixtures("nextflow_fixture")
def test_umi_group(configuration_module, sample_id):
    """
    Test that the UMI groups are as expected. See
    :py:func:`check_umi_groups`.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    tmp_dir = config[params.TMP_DIR]
    check_umi_groups(tmp_dir, sample_id, 5)


@pytest.mark.usefixtures("nextflow_fixture")
@pytest.mark.parametrize("sample_id", [riboviz.test.SIMDATA_UMI_SAMPLE])
@pytest.mark.parametrize("stats_file", ["edit_distance.tsv",
                                        "per_umi_per_position.tsv",
                                        "per_umi.tsv"])
def test_umi_stats_files(configuration_module, sample_id, stats_file):
    """
    Test that the UMI deduplication statistics files exist.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID
    :type sample_id: str or unicode
    :param stats_file: statistics file name
    :type stats_file: str or unicode
    """
    config, _ = configuration_module
    stats_file = os.path.join(config[params.TMP_DIR], sample_id,
                              workflow_files.DEDUP_STATS_FORMAT.format(
                                  stats_file))
    assert os.path.exists(stats_file)


def check_tpms_collated_tsv(output_dir, sample_id, expected_num_columns):
    """
    Test that the collated TPMs are as expected.

    :param config: configuration
    :type config: dict
    :param sample_id: sample ID
    :type sample_id: str or unicode
    :param expected_num_columns: Expected number of columns in TSV file
    :type expected_num_columns: int
    """
    tpms_tsv = os.path.join(output_dir, workflow_r.TPMS_COLLATED_TSV)
    tpms = pd.read_csv(tpms_tsv, sep="\t", comment="#")
    num_rows, num_columns = tpms.shape
    assert num_columns == expected_num_columns, \
        "Unexpected number of columns"
    assert num_rows == 68, "Unexpected number of rows"
    columns = list(tpms.columns)
    assert "ORF" in columns, "Missing 'ORF' column"
    assert sample_id in columns, "Missing {} column".format(sample_id)
    yal003w = tpms.loc[tpms["ORF"] == "YAL003W", sample_id]
    assert (yal003w == 607366.8).all(), \
        "{} value for ORF YAL003W is not 607366.8".format(sample_id)
    yal038w = tpms.loc[tpms["ORF"] == "YAL038W", sample_id]
    assert (yal038w == 392633.2).all(), \
        "{} value for ORF YAL038W is not 392633.2".format(sample_id)
    others = tpms[~tpms["ORF"].isin(["YAL003W", "YAL038W"])]
    assert (others[sample_id] == 0).all(), \
        "{} value for non-YAL003W and YAL038W ORFs are not 0".format(sample_id)


@pytest.mark.parametrize("sample_id", [riboviz.test.SIMDATA_UMI_SAMPLE])
@pytest.mark.usefixtures("nextflow_fixture")
def test_tpms_collated_tsv(configuration_module, sample_id):
    """
    Test that the collated TPMs are as expected. See
    :py:func:`check_tpms_collated_tsv`.

    :param configuration_module: temporary configuration and \
    configuration file
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    output_dir = config[params.OUTPUT_DIR]
    check_tpms_collated_tsv(output_dir, sample_id, 2)
