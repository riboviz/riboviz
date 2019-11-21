"""
riboviz.tools.prep_riboviz test suite to test adaptor trimming,
barcode and UMI extraction, demultpliexing and deduplication.

The test suite runs riboviz.tools.prep_riboviz using a copy of
"vignette/example_multiplex_config.yaml" and the simulated data in
"data/example/". It then validates the outputs of the adaptor
trimming, barcode and UMI extraction, demultiplexing and deduplication
steps against the expected outputs, also in "data/example/".

The simulated data in "data/example/" is expected to have been created
using riboviz.tools.create_fastq_examples.
"""
import os
import pytest
import pandas as pd
import riboviz
import riboviz.process_utils
import riboviz.test
import riboviz.tools
import riboviz.validation
from riboviz import params
from riboviz.tools import prep_riboviz
from riboviz.tools.demultiplex_fastq import NUM_READS_FILE
from riboviz.test.tools import configuration_module  # Test fixture


TEST_CONFIG_FILE = riboviz.test.EXAMPLE_MULTIPLEX_CONFIG
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
    Validate that adaptor trimming, performed by "cutadapt" produces
    the expected results.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    expected_output = os.path.join(
        riboviz.test.EXAMPLE_DATA_DIR,
        "multiplex_umi_barcode.fastq")
    actual_output = os.path.join(
        config[params.TMP_DIR],
        "multiplex_umi_barcode_adaptor_trim.fq")
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
    expected_output = os.path.join(riboviz.test.EXAMPLE_DATA_DIR,
                                   "multiplex.fastq")
    actual_output = os.path.join(
        config[params.TMP_DIR],
        "multiplex_umi_barcode_adaptor_extract_trim.fq")
    riboviz.validation.equal_fastq(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
def test_deplex_num_reads(configuration_module):
    """
    Validate that "num_reads.tsv", produced by
    riboviz.tools.demultiplex_fastq has the expected content.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    actual_dir = os.path.join(
        config[params.TMP_DIR],
        "multiplex_umi_barcode_adaptor_deplex")
    actual_output = os.path.join(actual_dir, NUM_READS_FILE)
    expected_output = os.path.join(
        riboviz.test.EXAMPLE_DATA_DIR, "deplex", NUM_READS_FILE)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.parametrize(
    "fastq", ["tag0.fastq", "tag1.fastq", "tag2.fastq", "Unassigned.fastq"])
@pytest.mark.usefixtures("run_prep_riboviz")
def test_deplex_reads(configuration_module, fastq):
    """
    Validate that ".fastq", produced by
    riboviz.tools.demultiplex_fastq have the expected content.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    :param fastq: FASTQ file
    :type fastq: str or unicode
    """
    config, _ = configuration_module
    actual_dir = os.path.join(
        config[params.TMP_DIR],
        "multiplex_umi_barcode_adaptor_deplex")
    actual_output = os.path.join(actual_dir, fastq)
    expected_output = os.path.join(
        riboviz.test.EXAMPLE_DATA_DIR, "deplex", fastq)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.parametrize("sample_id", ["tag0", "tag1", "tag2"])
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
    tmp_dir = config[params.TMP_DIR]
    groups_tsv = os.path.join(tmp_dir,
                              sample_id,
                              sample_id + "_" + "post_dedup_groups.tsv")
    groups = pd.read_csv(groups_tsv, sep="\t")
    num_groups = 5
    assert groups.shape[0] == num_groups, \
        ("Expected %d unique groups but found %d"
         % (num_groups, groups.shape[0]))
    assert (groups["umi_count"] == 1).all(), \
        "Expected each umi_count to be 1"
    assert (groups["final_umi_count"] == 1).all(), \
        "Expected each final_umi_count to be 1"
    # Check group IDs are unique when compared to 1,...,number of
    # groups.
    group_ids = list(groups["unique_id"])
    group_ids.sort()
    expected_group_ids = list(range(num_groups))
    assert expected_group_ids == group_ids, \
        ("Expected group_ids %s but found %s" % (str(expected_group_ids),
                                                 str(group_ids)))
    # Check each representative read does indeed come from a unique
    # UMI group by parsing the read ID. create_fastq_examples.py
    # creates read IDs of form:
    # "EWSim-<GROUP>.<MEMBER>-umi<5PRIME>-read<READ>-umi<3PRIME>"
    # where <GROUP> is 1-indexed.
    groups_from_read_ids = [
        int(read_id.split("-")[1].split(".")[0]) - 1
        for read_id in groups["read_id"]
    ]
    groups_from_read_ids.sort()
    assert groups_from_read_ids == group_ids, \
        ("Reads in read_ids %s are not from unique groups" %
         (str(list(groups["read_id"]))))


@pytest.mark.parametrize("sample_id", ["tag0", "tag1", "tag2"])
@pytest.mark.usefixtures("run_prep_riboviz")
def test_deplex_tpms_collated_tsv(configuration_module, sample_id):
    """
    Validate the "TPMs_collated.tsv" file produced by the workflow for
    each demultiplexed file.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    :param sample_id: sample ID for demultiplexed reads
    :type sample_id: str or unicode
    """
    config, _ = configuration_module
    output_dir = config[params.OUTPUT_DIR]
    tpms_tsv = os.path.join(output_dir,
                            sample_id,
                            sample_id + "_" + "TPMs_collated.tsv")
    tpms = pd.read_csv(tpms_tsv, sep="\t")
    num_rows, num_columns = tpms.shape
    assert num_columns == 2, "Unexpected number of columns"
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
