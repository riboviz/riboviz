"""
prep_riboviz.py test suite to test adatpor trimming, UMI extraction
and deduplication.

The test suite runs prep_riboviz.py using a copy of
`vignette/vignette_example_config.yaml` and the simulated data in
`data/example/`. It then validates the outputs of the
UMI-tools-specific phases against the expected outputs, also in
`data/`.

The simulated data in `data/` is expected to have been created using
`create_fastq_examples.py`.
"""
import os
import pytest
import pandas as pd
import riboviz
import riboviz.process_utils
import riboviz.test
import riboviz.tools
import riboviz.validation
from riboviz.tools import prep_riboviz
from riboviz.test.tools import configuration_module  # Test fixture


TEST_CONFIG_FILE = riboviz.test.EXAMPLE_DATA_CONFIG
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
    expected_output = os.path.join(riboviz.test.EXAMPLE_DATA_DIR,
                                   "example_umi5_umi3_umi.fastq")
    actual_output = os.path.join(config["dir_tmp"],
                                 "example_umi5_umi3_trim.fq")
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
    expected_output = os.path.join(riboviz.test.EXAMPLE_DATA_DIR,
                                   "example_umi5_umi3.fastq")
    actual_output = os.path.join(config["dir_tmp"],
                                 "example_umi5_umi3_extract_trim.fq")
    riboviz.validation.equal_fastq(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
def test_umi_group(configuration_module):
    """
    Validate the information on UMI groups post-`umi_tools extract`,
    by parsing the `.tsv` file output by `umi_tools group`.

    :param configuration_module: configuration and path to
    configuration file (pytest fixture)
    :type configuration_module: tuple(dict, str or unicode)
    """
    config, _ = configuration_module
    tmp_dir = config["dir_tmp"]
    groups_tsv = os.path.join(tmp_dir,
                              "example_umi5_umi3_post_dedup_groups.tsv")
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
