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


def tokenize_config_in_place(config):
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
    application of :py:func:`tokenize_config_in_place` and assuming when
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
    application of :py:func:`tokenize_config_in_place` and assuming when
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
    application of :py:func:`tokenize_config_in_place` and assuming when
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
    application of :py:func:`tokenize_config_in_place` and assuming when
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


TEST_DIR_FUNCTIONS = {
    test.VIGNETTE_CONFIG: create_vignette_test_dir,
    test.SIMDATA_UMI_CONFIG: create_simdata_test_dir,
    test.SIMDATA_MULTIPLEX_CONFIG: create_simdata_test_dir
}
"""
Dictionary mapping configurations to functions that create directories
consistent with that expected by those configurations.
"""
TOKEN_TEST_DIR_FUNCTIONS = {
    test.VIGNETTE_CONFIG: create_vignette_token_test_dir,
    test.SIMDATA_UMI_CONFIG: create_simdata_token_test_dir,
    test.SIMDATA_MULTIPLEX_CONFIG: create_simdata_token_test_dir
}
"""
Dictionary mapping configurations to functions that create directories
consistent with that expected by those configurations after
application of :py:func:`tokenize_config_in_place` and assuming when
a workflow is run the user provides values for the environment
variable (:py:const:`riboviz.params.ENV_DIRS`) corresponding to
the tokens.
"""
DEFAULT_TOKEN_TEST_DIR_FUNCTIONS = {
    test.VIGNETTE_CONFIG: create_vignette_default_token_test_dir,
    test.SIMDATA_UMI_CONFIG: create_simdata_default_token_test_dir,
    test.SIMDATA_MULTIPLEX_CONFIG: create_simdata_default_token_test_dir
}
"""
Dictionary mapping configurations to functiosn that create test
directories consistent with that expected by those configurations
after application of :py:func:`tokenize_config_in_place` and assuming
when a workflow is run the user does not provide values for the
environment variable (:py:const:`riboviz.params.ENV_DIRS`)
corresponding to the tokens.
"""


def get_empty_env_directory_map(organisms_dir, samples_dir, data_dir):
    """
    Get dictionary from environment variables
    (:py:const:`riboviz.params.ENV_DIRS`) to directories.
    Returns an empty dictionary. All parameters are unused.

    :param organisms_dir: Organisms directory
    :type organisms_dir: str or unicode
    :param samples_dir: Samples directory
    :type samples_dir: str or unicode
    :param data_dir: Data directory
    :type data_dir: str or unicode
    :return: Dictionary
    :rtype: dict
    """
    return {}


def get_env_directory_map(organisms_dir, samples_dir, data_dir):
    """
    Get dictionary from environment variables
    (:py:const:`riboviz.params.ENV_DIRS`) to directories.
    Returns:

    * :py:const.`riboviz.params.ENV_RIBOVIZ_SAMPLES` ``samples_dir``.
    * :py:const.`riboviz.params.ENV_RIBOVIZ_ORGANISMS` ``organisms_dir``.
    * :py:const.`riboviz.params.ENV_RIBOVIZ_DATA`: ``data_dir``.

    :param organisms_dir: Organisms directory
    :type organisms_dir: str or unicode
    :param samples_dir: Samples directory
    :type samples_dir: str or unicode
    :param data_dir: Data directory
    :type data_dir: str or unicode
    :return: Dictionary
    :rtype: dict
    """
    return {params.ENV_RIBOVIZ_ORGANISMS: organisms_dir,
            params.ENV_RIBOVIZ_SAMPLES: samples_dir,
            params.ENV_RIBOVIZ_DATA: data_dir}


def get_default_run_dir(directory):
    """
    Get directory within which to run workflow. Returns ``None``
    (for passing to :py:func:`riboviz.test.nextflow.run_nextflow`'s
    ``cwd`` argument. ``directory`` is unused.

    :param directory: Directory
    :type directory py._path.local.LocalPath
    :return: Directory
    :rtype: str or unicode
    """
    return None


def get_given_run_dir(directory):
    """
    Get directory within which to run workflow. Returns ``directory``
    (for passing to :py:func:`riboviz.test.nextflow.run_nextflow`'s
    ``cwd`` argument.

    :param directory: Directory
    :type directory py._path.local.LocalPath
    :return: Directory
    :rtype: str or unicode
    """
    return directory


def identity_config(config, organisms_dir, samples_dir, data_dir):
    """
    Leave configuration unchanged. All parameters are unused.

    :param config: Configuration
    :type config: dict
    :param organisms_dir: Organisms directory
    :type organisms_dir: str or unicode
    :param samples_dir: Samples directory
    :type samples_dir: str or unicode
    :param data_dir: Data directory
    :type data_dir: str or unicode
    """
    pass


def tokenize_config(config, organisms_dir, samples_dir, data_dir):
    """
    Tokenise configuration. See :py:func:`tokenize_config_in_place`.
    Other parameters are unused.

    :param config: Configuration
    :type config: dict
    :param organisms_dir: Organisms directory
    :type organisms_dir: str or unicode
    :param samples_dir: Samples directory
    :type samples_dir: str or unicode
    :param data_dir: Data directory
    :type data_dir: str or unicode
    """
    tokenize_config_in_place(config)


def relative_config(config, organisms_dir, samples_dir, data_dir):
    """
    Make configuration relative, replacing file and directory
    paths with relative paths. Other parameters are unused.

    `""` replaces the relative paths in:

    * :py:const.`riboviz.params.ASITE_DISP_LENGTH_FILE`
    * :py:const.`riboviz.params.CODON_POSITIONS_FILE`
    * :py:const.`riboviz.params.FEATURES_FILE`
    * :py:const.`riboviz.params.T_RNA_FILE`
    * :py:const.`riboviz.params.ORF_FASTA_FILE`
    * :py:const.`riboviz.params.ORF_GFF_FILE`
    * :py:const.`riboviz.params.RRNA_FASTA_FILE`
    * :py:const.`riboviz.params.INDEX_DIR`
    * :py:const.`riboviz.params.INPUT_DIR`
    * :py:const.`riboviz.params.OUTPUT_DIR`
    * :py:const.`riboviz.params.TMP_DIR`

    For example:

    * If :py:const:`riboviz.params.INPUT_DIR` has value
      ``data/simdata``, then it is updated to
      ``simdata``.
    * If :py:const.`riboviz.params.ORF_GFF_FILE` has value
      ``vignette/input/yeast_YAL_CDS_w_250utrs.gff3``, then it is
      updated to ``yeast_YAL_CDS_w_250utrs.gff3``.

    See also :py:`func:riboviz.test.customise_path`.

    :param config: Configuration
    :type config: dict
    :param organisms_dir: Organisms directory
    :type organisms_dir: str or unicode
    :param samples_dir: Samples directory
    :type samples_dir: str or unicode
    :param data_dir: Data directory
    :type data_dir: str or unicode
    """
    parameters = [
        params.ASITE_DISP_LENGTH_FILE,
        params.CODON_POSITIONS_FILE,
        params.FEATURES_FILE,
        params.T_RNA_FILE,
        params.ORF_FASTA_FILE,
        params.ORF_GFF_FILE,
        params.RRNA_FASTA_FILE,
        params.INDEX_DIR,
        params.INPUT_DIR,
        params.OUTPUT_DIR,
        params.TMP_DIR
    ]
    for param in parameters:
        config[param] = test.customise_path("", config[param])


def absolute_config(config, organisms_dir, samples_dir, data_dir):
    """
    Make configuration absolute, replacing file and directory
    paths with absolute paths.

    The following replacements are done:

    * ``<data_dir>`` replaces paths in:
      - :py:const.`riboviz.params.ASITE_DISP_LENGTH_FILE`
      - :py:const.`riboviz.params.CODON_POSITIONS_FILE`
      - :py:const.`riboviz.params.FEATURES_FILE`
      - :py:const.`riboviz.params.T_RNA_FILE`
    * ``<organisms_dir>`` replaces paths in:
      - :py:const.`riboviz.params.ORF_FASTA_FILE`
      - :py:const.`riboviz.params.ORF_GFF_FILE`
      - :py:const.`riboviz.params.RRNA_FASTA_FILE`
    * ``<samples_dir>`` replaces paths in:
      - :py:const.`riboviz.params.INDEX_DIR`
      - :py:const.`riboviz.params.INPUT_DIR`
      - :py:const.`riboviz.params.OUTPUT_DIR`
      - :py:const.`riboviz.params.TMP_DIR`

    For example:

    * If :py:const:`riboviz.params.INPUT_DIR` has value
      ``data/simdata``, then it is updated to
      ``<samples>/simdata``.
    * If :py:const.`riboviz.params.ORF_GFF_FILE` has value
      ``vignette/input/yeast_YAL_CDS_w_250utrs.gff3``, then it is
      updated to ``<organisms>/yeast_YAL_CDS_w_250utrs.gff3``.

    See also :py:`func:riboviz.test.customise_path`.

    :param config: Configuration
    :type config: dict
    """
    envs_params = {
        data_dir: [params.ASITE_DISP_LENGTH_FILE,
                   params.CODON_POSITIONS_FILE,
                   params.FEATURES_FILE,
                   params.T_RNA_FILE],
        organisms_dir: [params.ORF_FASTA_FILE,
                        params.ORF_GFF_FILE,
                        params.RRNA_FASTA_FILE],
        samples_dir: [params.INDEX_DIR,
                      params.INPUT_DIR,
                      params.OUTPUT_DIR,
                      params.TMP_DIR]
    }
    for directory, parameters in envs_params.items():
        for param in parameters:
            config[param] = test.customise_path(directory, config[param])


CONFIG_TEST_CASES = [
    (TEST_DIR_FUNCTIONS, get_empty_env_directory_map,
     identity_config, get_default_run_dir),
    # Tests that workflow configuration validates.
    # The workflow is run in the default,
    # :py.const:`riboviz.BASE_PATH`, directory.
    (TEST_DIR_FUNCTIONS, get_empty_env_directory_map,
     identity_config, get_given_run_dir),
    # Tests that workflow configuration validates in a test directory
    # whose structure matches that expected by the workflow
    # configuration.
    (TOKEN_TEST_DIR_FUNCTIONS, get_env_directory_map,
     tokenize_config, get_default_run_dir),
    # Test that a workflow configuration with environment variable
    # tokens validates when values for the environment variables
    # (:py:const:`riboviz.params.ENV_DIRS`) are provided.
    # The workflow is run in the default, :py.const:`riboviz.BASE_PATH`,
    # directory, with the environment variables specifying paths in a
    # test directory.
    (DEFAULT_TOKEN_TEST_DIR_FUNCTIONS, get_empty_env_directory_map,
     tokenize_config, get_given_run_dir),
    # Test that a workflow configuration with environment variable
    # tokens validates when values for the environment variables
    # (:py:const:`riboviz.params.ENV_DIRS`) are not provided, and
    # so the default value
    # (:py:const:`riboviz.environment.DEFAULT_ENV_DIR`)
    # is used for the environment variable values.
    # The workflow is run in a test directory.
    (DEFAULT_TOKEN_TEST_DIR_FUNCTIONS, get_empty_env_directory_map,
     relative_config, get_given_run_dir),
    # Test that a workflow configuration with relative paths to input
    # and output files validates.
    # The workflow is run in a test directory.
    (TOKEN_TEST_DIR_FUNCTIONS, get_empty_env_directory_map,
     absolute_config, get_default_run_dir)
    # Test that a workflow with absolute paths to input and output files
    # validates.
    # The workflow is run in the default, :py.const:`riboviz.BASE_PATH`,
    # directory with the configuration specifying paths in a test
    # directory.
]
"""
Configuration test cases. Each is a tuple with:
* Dictionary from configuration files to functions that create
  directory structures compatible with these configuration
  files.
* Function to return dictionary from environment variables to paths.
* Function to update a given configuration.
* Function to get directory within which to run workflow.
"""
CONFIG_TEST_CASE_ITEMS = \
    "test_dir_lookup, get_envs, update_config, get_run_dir"
""" Name of each item of a :py:const:`CONFIG_TEST_CASES` tuple. """
CONFIG_TEST_CASE_IDS = [
    "default",
    "default_test_dir",
    "environment_tokens",
    "default_environment_tokens",
    "relative_paths",
    "absolute_paths"
]
""" Identifiers for test cases in :py:const:`CONFIG_TEST_CASES`. """


@pytest.mark.parametrize(CONFIG_TEST_CASE_ITEMS,
                         CONFIG_TEST_CASES,
                         ids=CONFIG_TEST_CASE_IDS)
@pytest.mark.parametrize("config_file",
                         [test.VIGNETTE_CONFIG,
                          test.SIMDATA_UMI_CONFIG,
                          test.SIMDATA_MULTIPLEX_CONFIG])
def test_validate_config(tmpdir, config_file, test_dir_lookup,
                         get_envs, update_config, get_run_dir):
    """
    Run workflow to validate configuration.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param config_file: Configuration file
    :type config_file: str or unicode
    :param test_dir_lookup: Dictionary from configuration files to \
    functions that create directory structures compatible with \
    configuration files.
    :type test_dir_lookup: dict
    :param get_envs: Function to return dictionary from environment \
    variables to paths.
    :type get_ens: function
    :param update_config: Function to update a given configuration.
    :type update_config: function
    :param get_run_dir: Function to get directory within which to run \
    workflow.
    :type get_run_dir: function
    """
    organisms_dir, samples_dir, data_dir = \
        test_dir_lookup[config_file](tmpdir)
    envs = get_envs(organisms_dir, samples_dir, data_dir)
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    update_config(config, organisms_dir, samples_dir, data_dir)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    cwd = get_run_dir(tmpdir)
    exit_code = nextflow.run_nextflow(test_config_file,
                                      validate_only=True,
                                      envs=envs,
                                      cwd=cwd)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code
