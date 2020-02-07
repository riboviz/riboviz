"""
riboviz.tools.prep_riboviz test suite to test adaptor trimming, UMI
extraction and deduplication.

The test suite runs riboviz.tools.prep_riboviz using a copy of
"vignette/simdata_umi_config.yaml" and the simulated data in
"data/simdata/". It then validates the outputs of the adaptor
trimming, UMI extraction and deduplication steps against the expected
outputs, also in "data/simdata/".

The simulated data in "data/simdata/" is expected to have been created
using riboviz.tools.create_fastq_simdata.
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
from riboviz.tools import prep_riboviz
from riboviz.test.tools import configuration_module  # Test fixture
from riboviz.test.tools import run_prep_riboviz  # Test fixture


TEST_CONFIG_FILE = riboviz.test.SIMDATA_UMI_CONFIG
"""
YAML configuration used as a template configuration by these tests -
required by configuration test fixture
"""


@pytest.mark.parametrize("sample_id", [riboviz.test.SIMDATA_UMI_SAMPLE])
@pytest.mark.usefixtures("run_prep_riboviz")
def test_adaptor_trimming(configuration_module, sample_id):
    """
    Validate that adaptor trimming, performed by "cutadapt" produces
    the expected results.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
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
@pytest.mark.usefixtures("run_prep_riboviz")
def test_umi_extract(configuration_module, sample_id):
    """
    Validate that UMI extraction, performed by "umi_tools extract"
    produces the expected results.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
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


def check_umi_groups(config, sample_id, num_groups):
    """
    Validate the information on UMI groups post-"umi_tools extract",
    by parsing the ".tsv" file output by "umi_tools group".

    :param config: configuration
    :type config: dict
    :param sample_id: sample ID
    :type sample_id: str or unicode
    :param num_groups: Expected number of groups
    :type num_groups: int
    """
    tmp_dir = config[params.TMP_DIR]
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
    # UMI group by parsing the read ID. create_fastq_simdata.py
    # creates read IDs of form:
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
@pytest.mark.usefixtures("run_prep_riboviz")
def test_umi_group(configuration_module, sample_id):
    """
    Validate the information on UMI groups post-"umi_tools extract",
    by parsing the ".tsv" file output by "umi_tools group".

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    check_umi_groups(config, sample_id, 5)


def check_tpms_collated_tsv(config, sample_id, num_columns):
    """
    Validate the "TPMs_collated.tsv" file produced from running the
    workflow.

    :param config: configuration
    :type config: dict
    :param sample_id: sample ID
    :type sample_id: str or unicode
    :param num_columns: Expected number of columns in TSV file
    :type num_columns: int
    """
    output_dir = config[params.OUTPUT_DIR]
    tpms_tsv = os.path.join(output_dir, workflow_r.TPMS_COLLATED_TSV)
    tpms = pd.read_csv(tpms_tsv, sep="\t", comment="#")
    num_rows, num_columns = tpms.shape
    assert num_columns == num_columns, "Unexpected number of columns"
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
@pytest.mark.usefixtures("run_prep_riboviz")
def test_tpms_collated_tsv(configuration_module, sample_id):
    """
    Validate the "TPMs_collated.tsv" file produced from running the
    workflow.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    check_tpms_collated_tsv(config, sample_id, 2)
