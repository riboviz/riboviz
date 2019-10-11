"""
prep_riboviz.py test suite to test adatpor trimming, UMI extraction
and deduplication.

The test suite runs prep_riboviz.py using a copy of
`vignette/vignette_simdata_config.yaml` and the simulated data in
`data/`. It then validates the outputs of the UMI-tools-specific
phases against the expected outputs, also in `data/`.
"""
import os
import shutil
import tempfile
import pytest
import yaml
import riboviz
import riboviz.process_utils
import riboviz.test
import riboviz.tools
import riboviz.validation
from riboviz.tools import prep_riboviz


@pytest.fixture(scope="module")
def riboviz_run():
    """
    Fixture to run prep_riboviz.py.

    * Create temporary index, tmp, and output directories.
    * Create a copy of the simulated data configuration
      file, `vignette/vignette_simdata_config.yaml`, and update
      configuration to use the above directories.
    * Run prep_riboviz.py.
    * Provide configuration and configuration file path to each test
      function.
    * When tests complete, delete temporary directories and other
    * files that were created.

    :return: configuration, path to configuration file
    :rtype: tuple(dict, str or unicode)
    """
    config_yaml = os.path.join(riboviz.test.VIGNETTE_DIR,
                               "vignette_simdata_config.yaml")
    with open(config_yaml, "r") as f:
        config = yaml.load(f, yaml.SafeLoader)
    print(config)
    _, config_path = tempfile.mkstemp(prefix="tmp", suffix=".yaml")
    _, script_path = tempfile.mkstemp(prefix="tmp", suffix=".sh")
    index_dir = tempfile.mkdtemp("tmp_test_prep_riboviz_barcode_umi_index")
    tmp_dir = tempfile.mkdtemp("tmp_test_prep_riboviz_barcode_umi_tmp")
    out_dir = tempfile.mkdtemp("tmp_test_prep_riboviz_barcode_umi_out")
    logs_dir = tempfile.mkdtemp("tmp_test_prep_riboviz_barcode_umi_logs")

    config["dir_index"] = index_dir
    config["dir_tmp"] = tmp_dir
    config["dir_out"] = out_dir
    config["dir_logs"] = logs_dir
    config["cmd_file"] = script_path

    with open(config_path, 'w') as f:
        yaml.dump(config, f)

    exit_code = prep_riboviz.prep_riboviz(
        riboviz.test.PY_SCRIPTS,
        riboviz.test.R_SCRIPTS,
        config_path)
    assert exit_code == 0, \
        "prep_riboviz returned non-zero exit code %d" % exit_code

    yield (config, config_path)
    if os.path.exists(config_path):
        os.remove(config_path)
    if os.path.exists(script_path):
        os.remove(script_path)
    shutil.rmtree(index_dir)
    shutil.rmtree(tmp_dir)
    shutil.rmtree(out_dir)
    shutil.rmtree(logs_dir)


def test_adaptor_trimming(riboviz_run):
    """
    Validate that adaptor trimming, performed by `cutadapt` produces
    the expected results.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    """
    config, _ = riboviz_run
    expected_output = os.path.join(riboviz.test.DATA_DIR,
                                   "simdata_UMI5and3_4nt.fastq")
    actual_output = os.path.join(config["dir_tmp"],
                                 "simdata5and3_trim.fq")
    riboviz.validation.equal_fastq(expected_output, actual_output)


def test_umi_extract(riboviz_run):
    """
    Validate that UMI extraction, performed by `umi_tools extract`
    produces the expected results.

    :param configuration: configuration and path to configuration file
    (pytest fixture)
    :type configuration: tuple(dict, str or unicode)
    """
    config, _ = riboviz_run
    expected_output = os.path.join(riboviz.test.DATA_DIR,
                                   "simdata_extracted_UMI5and3_4nt.fastq")
    actual_output = os.path.join(config["dir_tmp"],
                                 "simdata5and3_extract_trim.fq")
    riboviz.validation.equal_fastq(expected_output, actual_output)
