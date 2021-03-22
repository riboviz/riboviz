"""
Nextflow-related test fixtures.
"""
import os
import os.path
import shutil
import tempfile
import subprocess
import pytest
import yaml
import riboviz
from riboviz import environment
from riboviz import params
from riboviz import test
from riboviz.test import nextflow


@pytest.fixture(name="test_dir", scope='module')
def test_dir_fixture(use_environment_variables):
    """
    Create temporary test directories. The directory structure is of form::

        <test-directory>
          data/
            yeast_codon_pos_i200.RData
            yeast_features.tsv
            yeast_standard_asite_disp_length.txt
            yeast_tRNAs.tsv
          organisms/
            yeast_rRNA_R64-1-1.fa
            yeast_YAL_CDS_w_250utrs.fa
            yeast_YAL_CDS_w_250utrs.gff3
          samples/
            simdata/
              umi5_umi3_umi_adaptor.fastq
              multiplex_umi_barcode_adaptor.fastq
              multiplex_barcodes.tsv

    The files are symbolic links to those in::

        data/
        vignette/input/
        data/simdata/

    :param use_environment_variables: Use environment variables? \
    This is unused but ensures that this fixture is reinvoked for \
    each value of use_environment_variables (fixture, see
    :py:module:`riboviz.test.nextflow.envs.confest`).
    :type use_environment_variables: bool
    :return: base test directory, organisms directory, samples \
    directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    print("---test_dir_fixture---")
    test_dir = tempfile.mkdtemp(prefix="test")
    organisms_dir = os.path.join(test_dir, "organisms")
    os.mkdir(organisms_dir)
    samples_dir = os.path.join(test_dir, "samples")
    os.makedirs(samples_dir)
    simdata_dir = os.path.join(samples_dir, "simdata")
    os.makedirs(simdata_dir)
    data_dir = os.path.join(test_dir, "data")
    os.mkdir(data_dir)
    os.symlink(
        os.path.join(test.VIGNETTE_DIR, "input", "yeast_rRNA_R64-1-1.fa"),
        os.path.join(organisms_dir, "yeast_rRNA_R64-1-1.fa"))
    os.symlink(
        os.path.join(test.VIGNETTE_DIR, "input", "yeast_YAL_CDS_w_250utrs.fa"),
        os.path.join(organisms_dir, "yeast_YAL_CDS_w_250utrs.fa"))
    os.symlink(
        os.path.join(test.VIGNETTE_DIR, "input", "yeast_YAL_CDS_w_250utrs.gff3"),
        os.path.join(organisms_dir, "yeast_YAL_CDS_w_250utrs.gff3"))
    os.symlink(
        os.path.join(test.DATA_DIR, "simdata", "umi5_umi3_umi_adaptor.fastq"),
        os.path.join(simdata_dir, "umi5_umi3_umi_adaptor.fastq"))
    os.symlink(
        os.path.join(test.DATA_DIR, "simdata", "multiplex_umi_barcode_adaptor.fastq"),
        os.path.join(simdata_dir, "multiplex_umi_barcode_adaptor.fastq"))
    os.symlink(
        os.path.join(test.DATA_DIR, "simddata", "multiplex_barcodes.tsv"),
        os.path.join(simdata_dir, "multiplex_barcodes.tsv"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_codon_pos_i200.RData"),
               os.path.join(data_dir, "yeast_codon_pos_i200.RData"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_features.tsv"),
               os.path.join(data_dir, "yeast_features.tsv"))
    os.symlink(
        os.path.join(test.DATA_DIR, "yeast_standard_asite_disp_length.txt"),
        os.path.join(data_dir, "yeast_standard_asite_disp_length.txt"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_tRNAs.tsv"),
               os.path.join(data_dir, "yeast_tRNAs.tsv"))
    print(str(organisms_dir))
    print(str(samples_dir))
    print(str(data_dir))
    yield test_dir, organisms_dir, samples_dir, data_dir
    # shutil.rmtree(test_dir)


@pytest.fixture(name="environment_vars", scope='module')
def environment_vars_fixture(use_environment_variables, test_dir):
    """
    Create a dictionary of environment variables to test directory.
    If ``use_environment_variables`` is ``True`` then return
    dictionary from
    :py:const:`riboviz.params.ENV_RIBOVIZ_SAMPLES`,
    :py:const:`riboviz.params.ENV_RIBOVIZ_ORGANISMS` and
    :py:const:`riboviz.params.ENV_RIBOVIZ_DATA` to complementary
    values in ``test_dir``.
    If ``use_environment_variables`` is ``False`` then return
    empty dictionary.
    :param use_environment_variables: Use environment variables?
    :type use_environment_variables: bool
    :param test_dir: base test directory, organisms directory, samples \
    directory, data directory (fixture)
    :type test_dir: tuple(str or unicode, str or unicode, str or \
    unicode)
    :return: dictionary from environment variables to values
    :rtype: dict
    """
    print("---environment_vars_fixture---")
    print(use_environment_variables)
    print(test_dir)
    if use_environment_variables:
        _, organisms_dir, samples_dir, data_dir = test_dir
        env_vars = {params.ENV_RIBOVIZ_SAMPLES: samples_dir,
                    params.ENV_RIBOVIZ_ORGANISMS: organisms_dir,
                    params.ENV_RIBOVIZ_DATA: data_dir}
    else:
        env_vars = {}
    print(env_vars)
    yield env_vars


@pytest.fixture(name="workflow_config", scope='module')
def workflow_config_fixture(request, environment_vars):
    """
    Get configuration file and customise configuration.

    Load a configuration file and replace any environment variable
    tokens with ``environment_vars``.

    Any test module using fixtures that call this function are
    must define a ``TEST_CONFIG_FILE`` parameter specifying the
    configuration file to be copied.

    :param request: Test module using this test fixture
    :type request: _pytest.fixtures.SubRequest
    :param environment_vars: dictionary from environment variables to \
    values (fixture)
    :type environment_vars: dict
    :return: configuration and configuration file
    :rtype: tuple(dict, str or unicode)
    """
    print("---workflow_config_fixture---")
    print(environment_vars)
    config_file = getattr(request.module, "TEST_CONFIG_FILE")
    print(config_file)
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    environment.update_config_with_env(environment_vars, config)
    print(config)
    yield config, config_file


@pytest.fixture(name="nextflow", scope="module")
def nextflow_fixture(workflow_config, test_dir, environment_vars):
    """
    Run ``nextflow run`` using a given conflguration file.

    :param workflow_config: configuration and configuration file (fixture)
    :rtype: tuple(dict, str or unicode)
    :param test_dir: base test directory, organisms directory, samples \
    directory, data directory (fixture)
    :type test_dir: tuple(str or unicode, str or unicode, str or \
    unicode)
    :param environment_vars: dictionary from environment variables to \
    values (fixture)
    :type environment_vars: dict
    :raise AssertionError: if the exit code is non-zero
    """
    _, config_file = workflow_config
    print("---nextflow_fixture---")
    print(test_dir)
    print(workflow_config)
    print(environment_vars)
    exit_code = nextflow.run_nextflow(config_file, environment_vars)
    exit_code = 0
    assert exit_code == 0, \
        "'nextflow run' returned non-zero exit code %d" % exit_code
