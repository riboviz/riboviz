"""
Nextflow configuration tests. These test the Nextflow workflow
validates configurations that include relative paths, absolute paths,
paths specifying environment variables (when values for the
environment variables are provided) and paths specifying environment
variables (when values for the environment variables are not
provided).

Each test runs ``nextflow run prep_riboviz.nf``.
"""
import yaml
import pytest
from riboviz import environment
from riboviz import params
from riboviz import test
from riboviz.test import nextflow


def tokenize_config(config):
    """
    Tokenize configuration in-place, replacing file and directory
    paths with environment variable tokens.

    The following replacements are done:

    * :py:const.`riboviz.params.ENV_RIBOVIZ_DATA` replaces paths in:
      - :py:const.`riboviz.params.ASITE_DISP_LENGTH_FILE`
      - :py:const.`riboviz.params.CODON_POSITIONS_FILE`
      - :py:const.`riboviz.params.FEATURES_FILE`
      - :py:const.`riboviz.params.T_RNA_FILE`
    * :py:const.`riboviz.params.ENV_RIBOVIZ_ORGANISMS` replaces paths
      in:
      - :py:const.`riboviz.params.ORF_FASTA_FILE`
      - :py:const.`riboviz.params.ORF_GFF_FILE`
      - :py:const.`riboviz.params.RRNA_FASTA_FILE`
    * :py:const.`riboviz.params.ENV_RIBOVIZ_SAMPLES`replaces paths in:
      - :py:const.`riboviz.params.INDEX_DIR`
      - :py:const.`riboviz.params.INPUT_DIR`
      - :py:const.`riboviz.params.OUTPUT_DIR`
      - :py:const.`riboviz.params.TMP_DIR`

    For example:

    * If :py:const:`riboviz.params.INPUT_DIR` has value
      ``data/simdata``, then it is updated to
      ``${RIBOVIZ_SAMPLES}/simdata``.
    * If :py:const.`riboviz.params.ORF_GFF_FILE` has value
      ``vignette/input/yeast_YAL_CDS_w_250utrs.gff3``, then it is
      updated to
      ``${RIBOVIZ_ORGANISMS}/yeast_YAL_CDS_w_250utrs.gff3``.

    See also :py:`func:riboviz.test.customise_path`.

    :param config: Configuration
    :type config: dict
    """
    envs_params = {
        params.ENV_RIBOVIZ_DATA: [params.ASITE_DISP_LENGTH_FILE,
                                  params.CODON_POSITIONS_FILE,
                                  params.FEATURES_FILE,
                                  params.T_RNA_FILE],
        params.ENV_RIBOVIZ_ORGANISMS: [params.ORF_FASTA_FILE,
                                       params.ORF_GFF_FILE,
                                       params.RRNA_FASTA_FILE],
        params.ENV_RIBOVIZ_SAMPLES: [params.INDEX_DIR,
                                     params.INPUT_DIR,
                                     params.OUTPUT_DIR,
                                     params.TMP_DIR]
    }
    for env, parameters in envs_params.items():
        for param in parameters:
            config[param] = test.customise_path(
                environment.ENV_TOKEN_FORMAT.format(env),
                config[param])


def create_vignette_test_dir(directory):
    """
    Create temporary test directories. The directory structure is of
    consistent with that expected by
    :py:const:`ribeviz.test.VIGNETTE_CONFIG`. The
    directory is structured as follows::

        <directory>
          data/
            # Symholic links to <riboviz>/data/ files:
            yeast_codon_pos_i200.RData
            yeast_features.tsv
            yeast_standard_asite_disp_length.txt
            yeast_tRNAs.tsv
          vignette/
            input/
              # Symholic links to <riboviz>/vignette/input/ files:
              yeast_rRNA_R64-1-1.fa
              yeast_YAL_CDS_w_250utrs.fa
              yeast_YAL_CDS_w_250utrs.gff3
              # Symholic links to <riboviz>/vignette/input/ files:
              SRR1042855_s1mi.fastq.gz
              SRR1042864_s1mi.fastq.gz

    The directories returned are as follows:

    * organisms: ``<directory>/vignette/input/``
    * samples: ``<directory>/vignette/``
    * data: ``<directory>/data/``

    :param directory: Directory
    :type directory py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    vignette_dir = directory.mkdir("vignette")
    input_dir = vignette_dir.mkdir("input")
    data_dir = directory.mkdir("data")
    test.symlink_files(input_dir, test.ORGANISM_FILES)
    test.symlink_files(input_dir, test.VIGNETTE_INPUT_FILES)
    test.symlink_files(data_dir, test.DATA_FILES)
    return str(input_dir), str(vignette_dir), str(data_dir)


def create_vignette_token_test_dir(directory):
    """
    Create temporary test directories. The directory structure is of
    consistent with that expected by
    :py:const:`riboviz.test.VIGNETTE_CONFIG` after
    application of :py:func:`tokenize_config` and assuming when
    a workflow is run the user provides values for the environment
    variable (:py:const:`riboviz.params.ENV_DIRS`) corresponding to
    the tokens. The directory is structured as follows::

        <directory>
          data/
            # Symholic links to <riboviz>/data/ files:
            yeast_codon_pos_i200.RData
            yeast_features.tsv
            yeast_standard_asite_disp_length.txt
            yeast_tRNAs.tsv
          organisms/
            # Symholic links to <riboviz>/vignette/input/ files:
            yeast_rRNA_R64-1-1.fa
            yeast_YAL_CDS_w_250utrs.fa
            yeast_YAL_CDS_w_250utrs.gff3
          samples/
            input/
                # Symholic links to <riboviz>/vignette/input/ files:
                SRR1042855_s1mi.fastq.gz
                SRR1042864_s1mi.fastq.gz

    The directories returned are as follows:

    * organisms: ``<directory>/organisms/``
    * samples: ``<directory>/samples/``
    * data: ``<directory>/data/``

    :param directory: Directory
    :type directory py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    organisms_dir = directory.mkdir("organisms")
    samples_dir = directory.mkdir("samples")
    input_dir = samples_dir.mkdir("input")
    data_dir = directory.mkdir("data")
    test.symlink_files(organisms_dir, test.ORGANISM_FILES)
    test.symlink_files(input_dir, test.VIGNETTE_INPUT_FILES)
    test.symlink_files(data_dir, test.DATA_FILES)
    return str(organisms_dir), str(samples_dir), str(data_dir)


def create_vignette_default_token_test_dir(directory):
    """
    Create temporary test directories. The directory structure is of
    consistent with that expected by
    :py:const:`riboviz.test.VIGNETTE_CONFIG` after
    application of :py:func:`tokenize_config` and assuming when
    a workflow is run the user does not provide values for the
    environment variable (:py:const:`riboviz.params.ENV_DIRS`)
    corresponding to the tokens. The directory is structured as
    follows::

        <directory>
          # Symholic links to <riboviz>/data/ files:
          yeast_codon_pos_i200.RData
          yeast_features.tsv
          yeast_standard_asite_disp_length.txt
          yeast_tRNAs.tsv
          # Symholic links to <riboviz>/vignette/input/ files:
          yeast_rRNA_R64-1-1.fa
          yeast_YAL_CDS_w_250utrs.fa
          yeast_YAL_CDS_w_250utrs.gff3
          input/
            # Symholic links to <riboviz>/vignette/input/ files:
            SRR1042855_s1mi.fastq.gz
            SRR1042864_s1mi.fastq.gz

    The directories returned are as follows:

    * organisms: ``<directory>/``
    * samples: ``<directory>/input/``
    * data: ``<directory>/``

    :param directory: Directory
    :type directory py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    input_dir = directory.mkdir("input")
    test.symlink_files(directory, test.ORGANISM_FILES)
    test.symlink_files(input_dir, test.VIGNETTE_INPUT_FILES)
    test.symlink_files(directory, test.DATA_FILES)
    return str(directory), str(input_dir), str(directory)


def create_simdata_test_dir(directory):
    """
    Create temporary test directories. The directory structure is of
    consistent with that expected by
    :py:const:`riboviz.test.SIMDATA_UMI_CONFIG`
    and :py:const:`riboviz.test.SIMDATA_MULTIPLEX_CONFIG`. The
    directory is structured as follows::

        <directory>
          data/
            simdata/
              # Symholic links to <riboviz>/data/simdata/ files:
              umi5_umi3_umi_adaptor.fastq
              multiplex_umi_barcode_adaptor.fastq
              multiplex_barcodes.tsv
            # Symholic links to <riboviz>/data/ files:
            yeast_codon_pos_i200.RData
            yeast_features.tsv
            yeast_standard_asite_disp_length.txt
            yeast_tRNAs.tsv
          vignette/
            input/
              # Symholic links to <riboviz>/vignette/input/ files:
              yeast_rRNA_R64-1-1.fa
              yeast_YAL_CDS_w_250utrs.fa
              yeast_YAL_CDS_w_250utrs.gff3

    The directories returned are as follows:

    * organisms: ``<directory>/vignette/input/``
    * samples: ``<directory>/data/``
    * data: ``<directory>/data/``

    :param directory: Directory
    :type directory py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    input_dir = directory.mkdir("vignette").mkdir("input")
    data_dir = directory.mkdir("data")
    simdata_dir = data_dir.mkdir("simdata")
    test.symlink_files(input_dir, test.ORGANISM_FILES)
    test.symlink_files(simdata_dir, test.SIMDATA_INPUT_FILES)
    test.symlink_files(data_dir, test.DATA_FILES)
    return str(input_dir), str(simdata_dir), str(data_dir)


def create_simdata_token_test_dir(directory):
    """
    Create temporary test directories. The directory structure is of
    consistent with that expected by
    :py:const:`riboviz.test.SIMDATA_UMI_CONFIG`
    and :py:const:`riboviz.test.SIMDATA_MULTIPLEX_CONFIG` after
    application of :py:func:`tokenize_config` and assuming when
    a workflow is run the user provides values for the environment
    variable (:py:const:`riboviz.params.ENV_DIRS`) corresponding to
    the tokens. The directory is structured as follows::

        <directory>
          data/
            # Symholic links to <riboviz>/data/ files:
            yeast_codon_pos_i200.RData
            yeast_features.tsv
            yeast_standard_asite_disp_length.txt
            yeast_tRNAs.tsv
          organisms/
            # Symholic links to <riboviz>/vignette/input/ files:
            yeast_rRNA_R64-1-1.fa
            yeast_YAL_CDS_w_250utrs.fa
            yeast_YAL_CDS_w_250utrs.gff3
          samples/
            simdata/
              # Symholic links to <riboviz>/data/simdata/ files:
              umi5_umi3_umi_adaptor.fastq
              multiplex_umi_barcode_adaptor.fastq
              multiplex_barcodes.tsv

    The directories returned are as follows:

    * organisms: ``<directory>/organisms/``
    * samples: ``<directory>/samples/``
    * data: ``<directory>/data/``

    :param directory: Directory
    :type directory py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    organisms_dir = directory.mkdir("organisms")
    samples_dir = directory.mkdir("samples")
    simdata_dir = samples_dir.mkdir("simdata")
    data_dir = directory.mkdir("data")
    test.symlink_files(organisms_dir, test.ORGANISM_FILES)
    test.symlink_files(simdata_dir, test.SIMDATA_INPUT_FILES)
    test.symlink_files(data_dir, test.DATA_FILES)
    return str(organisms_dir), str(samples_dir), str(data_dir)


def create_simdata_default_token_test_dir(directory):
    """
    Create temporary test directories. The directory structure is of
    consistent with that expected by
    :py:const:`riboviz.test.SIMDATA_UMI_CONFIG`
    and :py:const:`riboviz.test.SIMDATA_MULTIPLEX_CONFIG` after
    application of :py:func:`tokenize_config` and assuming when
    a workflow is run the user does not provide values for the
    environment variable (:py:const:`riboviz.params.ENV_DIRS`)
    corresponding to the tokens. The directory is structured as
    follows::

        <directory>
          # Symholic links to <riboviz>/data/ files:
          yeast_codon_pos_i200.RData
          yeast_features.tsv
          yeast_standard_asite_disp_length.txt
          yeast_tRNAs.tsv
          # Symholic links to <riboviz>/vignette/input/ files:
          yeast_rRNA_R64-1-1.fa
          yeast_YAL_CDS_w_250utrs.fa
          yeast_YAL_CDS_w_250utrs.gff3
          simdata/
            # Symholic links to <riboviz>/data/simdata/ files:
            umi5_umi3_umi_adaptor.fastq
            multiplex_umi_barcode_adaptor.fastq
            multiplex_barcodes.tsv

    The directories returned are as follows:

    * organisms: ``<directory>/``
    * samples: ``<directory>/simdata/``
    * data: ``<directory>/``

    :param directory: Directory
    :type directory py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    simdata_dir = directory.mkdir("simdata")
    test.symlink_files(directory, test.ORGANISM_FILES)
    test.symlink_files(simdata_dir, test.SIMDATA_INPUT_FILES)
    test.symlink_files(directory, test.DATA_FILES)
    return str(directory), str(simdata_dir), str(directory)


@pytest.mark.parametrize(
    "config_file, create_test_dir",
    [(test.VIGNETTE_CONFIG, create_vignette_test_dir),
     (test.SIMDATA_UMI_CONFIG, create_simdata_test_dir),
     (test.SIMDATA_MULTIPLEX_CONFIG, create_simdata_test_dir)])
def test_config_test_dir(tmpdir, config_file, create_test_dir):
    """
    Test that workflow configuration validates in a test
    directory whose structure matches that expected by the workflow
    configuration.
    The workflow is run in a test directory (see
    :py:func:`create_simdata_test_dir`).

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param config_file: Configuration file
    :type config_file: str or unicode
    :param create_test_dir: Callback function to create test directory
    :type create_test_dir: function
    """
    _ = create_test_dir(tmpdir)
    exit_code = nextflow.run_nextflow(config_file,
                                      validate_only=True,
                                      cwd=tmpdir)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize(
    "config_file, create_test_dir",
    [(test.VIGNETTE_CONFIG, create_vignette_token_test_dir),
     (test.SIMDATA_UMI_CONFIG, create_simdata_token_test_dir),
     (test.SIMDATA_MULTIPLEX_CONFIG, create_simdata_token_test_dir)])
def test_environment_token_config(tmpdir, config_file, create_test_dir):
    """
    Test that a workflow configuration with environment variable
    tokens validates when the user provides values for the environment
    variables (:py:const:`riboviz.params.ENV_DIRS`).
    The workflow is run in the default, :py.const:`riboviz.BASE_PATH`,
    directory, with the environment variables specifying paths in a
    test directory, (see :py:func:'create_simdata_token_test_dir`).

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param config_file: Configuration file
    :type config_file: str or unicode
    :param create_test_dir: Callback function to create test directory
    :type create_test_dir: function
    """
    organisms_dir, samples_dir, data_dir = create_test_dir(tmpdir)
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    tokenize_config(config)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    envs = {params.ENV_RIBOVIZ_SAMPLES: samples_dir,
            params.ENV_RIBOVIZ_ORGANISMS: organisms_dir,
            params.ENV_RIBOVIZ_DATA: data_dir}
    exit_code = nextflow.run_nextflow(test_config_file,
                                      envs=envs,
                                      validate_only=True)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize(
    "config_file, create_test_dir",
    [(test.VIGNETTE_CONFIG, create_vignette_default_token_test_dir),
     (test.SIMDATA_UMI_CONFIG, create_simdata_default_token_test_dir),
     (test.SIMDATA_MULTIPLEX_CONFIG, create_simdata_default_token_test_dir)])
def test_default_environment_token_config(tmpdir, config_file,
                                          create_test_dir):
    """
    Test that a workflow configuration with environment variable
    tokens validates when the user does not provide values for the
    environment  variables (:py:const:`riboviz.params.ENV_DIRS`), and
    so the default value
    (:py:const:`riboviz.environment.DEFAULT_ENV_DIR`)
     is used for the environment variables.
    The workflow is run in a test directory (see
    :py:func:'create_simdata_default_token_test_dir`).

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param config_file: Configuration file
    :type config_file: str or unicode
    :param create_test_dir: Callback function to create test directory
    :type create_test_dir: function
    """
    _ = create_test_dir(tmpdir)
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    tokenize_config(config)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    envs = {}
    exit_code = nextflow.run_nextflow(test_config_file,
                                      envs=envs,
                                      validate_only=True,
                                      cwd=tmpdir)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize(
    "config_file, create_test_dir",
    [(test.VIGNETTE_CONFIG, create_vignette_default_token_test_dir),
     (test.SIMDATA_UMI_CONFIG, create_simdata_default_token_test_dir),
     (test.SIMDATA_MULTIPLEX_CONFIG, create_simdata_default_token_test_dir)])
def test_relative_paths_config(tmpdir, config_file, create_test_dir):
    """
    Test that a workflow configuration with relative paths to input
    and output files validates.
    The workflow is run in a test directory (see
    :py:func:'create_simdata_default_token_test_dir`).

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param config_file: Configuration file
    :type config_file: str or unicode
    :param create_test_dir: Callback function to create test directory
    :type create_test_dir: function
    """
    _ = create_test_dir(tmpdir)
    # Tokenise configuration then update configuration with no
    # values for environment variable tokens provided so paths become
    # relative.
    # This is a convenience. Note that no environment variables are
    # passed to run_nextflow below.
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    tokenize_config(config)
    environment.update_config_with_env({}, config)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = nextflow.run_nextflow(test_config_file,
                                      validate_only=True,
                                      cwd=tmpdir)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize(
    "config_file, create_test_dir",
    [(test.VIGNETTE_CONFIG, create_vignette_token_test_dir),
     (test.SIMDATA_UMI_CONFIG, create_simdata_token_test_dir),
     (test.SIMDATA_MULTIPLEX_CONFIG, create_simdata_token_test_dir)])
def test_absolute_paths_config(tmpdir, config_file, create_test_dir):
    """
    Test that a workflow with absolute paths to input and output files
    validates.
    The workflow is run in the default, :py.const:`riboviz.BASE_PATH`,
    directory with the configuration specifying paths in a test
    directory (see :py:func:'create_simdata_token_test_dir`).

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param config_file: Configuration file
    :type config_file: str or unicode
    :param create_test_dir: Callback function to create test directory
    :type create_test_dir: function
    """
    organisms_dir, samples_dir, data_dir = create_test_dir(tmpdir)
    # Tokenise configuration then update configuration with
    # values for environment variable tokens provided so paths become
    # absolute.
    # This is a convenience. Note that no environment variables are
    # passed to run_nextflow below.
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    tokenize_config(config)
    envs = {params.ENV_RIBOVIZ_SAMPLES: samples_dir,
            params.ENV_RIBOVIZ_ORGANISMS: organisms_dir,
            params.ENV_RIBOVIZ_DATA: data_dir}
    environment.update_config_with_env(envs, config)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = nextflow.run_nextflow(test_config_file,
                                      validate_only=True)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize(
    "config_file",
    [test.VIGNETTE_CONFIG,
     test.SIMDATA_UMI_CONFIG,
     test.SIMDATA_MULTIPLEX_CONFIG])
def test_config_riboviz_dir(config_file):
    """
    Test that a workflow configuration validates.
    The workflow is run in the default,
    :py.const:`riboviz.BASE_PATH`, directory.

    :param config_file: Configuration file
    :type config_file: str or unicode
    """
    exit_code = nextflow.run_nextflow(config_file,
                                      validate_only=True)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize(
    "config_file, samples_dir",
    [(test.VIGNETTE_CONFIG, test.VIGNETTE_DIR),
     (test.SIMDATA_UMI_CONFIG, test.DATA_DIR),
     (test.SIMDATA_MULTIPLEX_CONFIG, test.DATA_DIR)])
def test_environment_token_config_riboviz_dir(
        tmpdir, config_file, samples_dir):
    """
    Test that a workflow configuration with environment variable
    tokens validates when the user provides values for the environment
    variables (:py:const:`riboviz.params.ENV_DIRS`).
    The workflow is run in ``tmpdir``, with the configuration, with
    the environment variables specifying paths in the default,
    :py:const:`riboviz.BASE_PATH`, directory.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param config_file: Configuration file
    :type config_file: str or unicode
    :param samples_dir: Samples directory
    :type samples_dir: str or unicode
    """
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    tokenize_config(config)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    envs = {params.ENV_RIBOVIZ_SAMPLES: samples_dir,
            params.ENV_RIBOVIZ_ORGANISMS: test.VIGNETTE_INPUT_DIR,
            params.ENV_RIBOVIZ_DATA: test.DATA_DIR}
    exit_code = nextflow.run_nextflow(test_config_file,
                                      envs=envs,
                                      validate_only=True)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code
