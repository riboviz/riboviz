"""
Nextflow configuration tests. These test the Nextflow workflow
validates configurations that include relative paths, absolute paths,
paths specifying environment variables (when values for the
environment variables are provided) and paths specifying environment
variables (when values for the environment variables are not
provided).

Each test runs ``nextflow run prep_riboviz.nf``.
"""
import os
import os.path
import yaml
from riboviz import environment
from riboviz import params
from riboviz import test
from riboviz.test import nextflow


def customise_simdata_config(config):
    config[params.ASITE_DISP_LENGTH_FILE] = "${RIBOVIZ_DATA}/yeast_standard_asite_disp_length.txt"
    config[params.CODON_POSITIONS_FILE] = None
    config[params.FEATURES_FILE] = "${RIBOVIZ_DATA}/yeast_features.tsv"
    config[params.INDEX_DIR] = "${RIBOVIZ_SAMPLES}/index"
    config[params.INPUT_DIR] = "${RIBOVIZ_SAMPLES}/input"
    config[params.OUTPUT_DIR] = "${RIBOVIZ_SAMPLES}/output"
    config[params.ORF_FASTA_FILE] = "${RIBOVIZ_ORGANISMS}/yeast_YAL_CDS_w_250utrs.fa"
    config[params.ORF_GFF_FILE] = "${RIBOVIZ_ORGANISMS}/yeast_YAL_CDS_w_250utrs.gff3"
    config[params.RRNA_FASTA_FILE] = "${RIBOVIZ_ORGANISMS}/yeast_rRNA_R64-1-1.fa"
    config[params.T_RNA_FILE] = None
    config[params.TMP_DIR] = "${RIBOVIZ_SAMPLES}/tmp"


def customise_vignette_config(config):
    config[params.ASITE_DISP_LENGTH_FILE] = "${RIBOVIZ_DATA}/yeast_standard_asite_disp_length.txt"
    config[params.CODON_POSITIONS_FILE] = "${RIBOVIZ_DATA}/yeast_codon_pos_i200.RData"
    config[params.FEATURES_FILE] = "${RIBOVIZ_DATA}/yeast_features.tsv"
    config[params.INDEX_DIR] = "${RIBOVIZ_SAMPLES}/index"
    config[params.INPUT_DIR] = "${RIBOVIZ_SAMPLES}/input"
    config[params.OUTPUT_DIR] = "${RIBOVIZ_SAMPLES}/output"
    config[params.ORF_FASTA_FILE] = "${RIBOVIZ_ORGANISMS}/yeast_YAL_CDS_w_250utrs.fa"
    config[params.ORF_GFF_FILE] = "${RIBOVIZ_ORGANISMS}/yeast_YAL_CDS_w_250utrs.gff3"
    config[params.RRNA_FASTA_FILE] = "${RIBOVIZ_ORGANISMS}/yeast_rRNA_R64-1-1.fa"
    config[params.T_RNA_FILE] = "${RIBOVIZ_DATA}/yeast_tRNAs.tsv"
    config[params.TMP_DIR] = "${RIBOVIZ_SAMPLES}/tmp"


# TODO refector out commonality, pass in subdirectories
# as arguments.
def create_test_dirs_envs(tmpdir):
    """
    Create temporary test directories. The directory structure is of
    form::

        <tmpdir>
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
              # Symholic links to <riboviz>/data/simdata/ files:
              umi5_umi3_umi_adaptor.fastq
              multiplex_umi_barcode_adaptor.fastq
              multiplex_barcodes.tsv

    The directories returned are as follows:

    * organisms: ``<tmpdir>/organisms/``
    * samples: ``<tmpdir>/samples/``
    * data: ``<tmpdir>/data/``

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    print("---create-test-dirs---")
    organisms_dir = tmpdir.mkdir("organisms")
    samples_dir = tmpdir.mkdir("samples")
    input_dir = samples_dir.mkdir("input")
    data_dir = tmpdir.mkdir("data")
    print(str(tmpdir))
    print(str(organisms_dir))
    print(str(samples_dir))
    print(str(data_dir))
    src_vignette_input_dir = os.path.join(test.VIGNETTE_DIR, "input")
    src_simdata_dir = os.path.join(test.DATA_DIR, "simdata")
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_rRNA_R64-1-1.fa"),
        organisms_dir.join("yeast_rRNA_R64-1-1.fa"))
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_YAL_CDS_w_250utrs.fa"),
        organisms_dir.join("yeast_YAL_CDS_w_250utrs.fa"))
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_YAL_CDS_w_250utrs.gff3"),
        organisms_dir.join("yeast_YAL_CDS_w_250utrs.gff3"))
    os.symlink(
        os.path.join(src_simdata_dir, "umi5_umi3_umi_adaptor.fastq"),
        input_dir.join("umi5_umi3_umi_adaptor.fastq"))
    os.symlink(
        os.path.join(src_simdata_dir, "multiplex_umi_barcode_adaptor.fastq"),
        input_dir.join("multiplex_umi_barcode_adaptor.fastq"))
    os.symlink(
        os.path.join(src_simdata_dir, "multiplex_barcodes.tsv"),
        input_dir.join("multiplex_barcodes.tsv"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_codon_pos_i200.RData"),
               data_dir.join("yeast_codon_pos_i200.RData"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_features.tsv"),
               data_dir.join("yeast_features.tsv"))
    os.symlink(
        os.path.join(test.DATA_DIR, "yeast_standard_asite_disp_length.txt"),
        data_dir.join("yeast_standard_asite_disp_length.txt"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_tRNAs.tsv"),
               data_dir.join("yeast_tRNAs.tsv"))
    return str(organisms_dir), str(samples_dir), str(data_dir)


def create_test_dirs_envs_default(tmpdir):
    """
    Create temporary test directories. The directory structure is of
    form::

        <tmpdir>
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
            # Symholic links to <riboviz>/data/simdata/ files:
            umi5_umi3_umi_adaptor.fastq
            multiplex_umi_barcode_adaptor.fastq
            multiplex_barcodes.tsv

    The directories returned are as follows:

    * organisms: ``<tmpdir>/``
    * samples: ``<tmpdir>/input/``
    * data: ``<tmpdir>``

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    print("---create-test-dirs---")
    organisms_dir = tmpdir
    samples_dir = tmpdir.mkdir("input")
    data_dir = tmpdir
    print(str(tmpdir))
    print(str(organisms_dir))
    print(str(samples_dir))
    print(str(data_dir))
    src_vignette_input_dir = os.path.join(test.VIGNETTE_DIR, "input")
    src_simdata_dir = os.path.join(test.DATA_DIR, "simdata")
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_rRNA_R64-1-1.fa"),
        organisms_dir.join("yeast_rRNA_R64-1-1.fa"))
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_YAL_CDS_w_250utrs.fa"),
        organisms_dir.join("yeast_YAL_CDS_w_250utrs.fa"))
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_YAL_CDS_w_250utrs.gff3"),
        organisms_dir.join("yeast_YAL_CDS_w_250utrs.gff3"))

    os.symlink(
        os.path.join(src_simdata_dir, "umi5_umi3_umi_adaptor.fastq"),
        samples_dir.join("umi5_umi3_umi_adaptor.fastq"))
    os.symlink(
        os.path.join(src_simdata_dir, "multiplex_umi_barcode_adaptor.fastq"),
        samples_dir.join("multiplex_umi_barcode_adaptor.fastq"))
    os.symlink(
        os.path.join(src_simdata_dir, "multiplex_barcodes.tsv"),
        samples_dir.join("multiplex_barcodes.tsv"))

    os.symlink(os.path.join(test.DATA_DIR, "yeast_codon_pos_i200.RData"),
               data_dir.join("yeast_codon_pos_i200.RData"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_features.tsv"),
               data_dir.join("yeast_features.tsv"))
    os.symlink(
        os.path.join(test.DATA_DIR, "yeast_standard_asite_disp_length.txt"),
        data_dir.join("yeast_standard_asite_disp_length.txt"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_tRNAs.tsv"),
               data_dir.join("yeast_tRNAs.tsv"))
    return str(organisms_dir), str(samples_dir), str(data_dir)


def create_test_dirs(tmpdir):
    """
    Create temporary test directories. The directory structure is of
    form::

        <tmpdir>
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

    * organisms: ``<tmpdir>/vignette/input/``
    * samples: ``<tmpdir>/data/``
    * data: ``<tmpdir>/data/``

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :return: organisms directory, samples directory, data directory
    :rtype: tuple(str or unicode, str or unicode, str or unicode)
    """
    print("---create-test-dirs---")
    organisms_dir = tmpdir.mkdir("vignette").mkdir("input")
    samples_dir = tmpdir.mkdir("data")
    simdata_dir = samples_dir.mkdir("simdata")
    data_dir = samples_dir
    print(str(tmpdir))
    print(str(organisms_dir))
    print(str(samples_dir))
    print(str(data_dir))
    src_vignette_input_dir = os.path.join(test.VIGNETTE_DIR, "input")
    src_simdata_dir = os.path.join(test.DATA_DIR, "simdata")
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_rRNA_R64-1-1.fa"),
        organisms_dir.join("yeast_rRNA_R64-1-1.fa"))
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_YAL_CDS_w_250utrs.fa"),
        organisms_dir.join("yeast_YAL_CDS_w_250utrs.fa"))
    os.symlink(
        os.path.join(src_vignette_input_dir, "yeast_YAL_CDS_w_250utrs.gff3"),
        organisms_dir.join("yeast_YAL_CDS_w_250utrs.gff3"))
    os.symlink(
        os.path.join(src_simdata_dir, "umi5_umi3_umi_adaptor.fastq"),
        simdata_dir.join("umi5_umi3_umi_adaptor.fastq"))
    os.symlink(
        os.path.join(src_simdata_dir, "multiplex_umi_barcode_adaptor.fastq"),
        simdata_dir.join("multiplex_umi_barcode_adaptor.fastq"))
    os.symlink(
        os.path.join(src_simdata_dir, "multiplex_barcodes.tsv"),
        simdata_dir.join("multiplex_barcodes.tsv"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_codon_pos_i200.RData"),
               data_dir.join("yeast_codon_pos_i200.RData"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_features.tsv"),
               data_dir.join("yeast_features.tsv"))
    os.symlink(
        os.path.join(test.DATA_DIR, "yeast_standard_asite_disp_length.txt"),
        data_dir.join("yeast_standard_asite_disp_length.txt"))
    os.symlink(os.path.join(test.DATA_DIR, "yeast_tRNAs.tsv"),
               data_dir.join("yeast_tRNAs.tsv"))
    return str(organisms_dir), str(samples_dir), str(data_dir)


def test_environment_vars(tmpdir):
    """
    Test that a workflow with environment variables,
    :py:const:`riboviz.params.ENV_DIRS`, validates.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    organisms_dir, samples_dir, data_dir = create_test_dirs_envs(tmpdir)

    with open(test.SIMDATA_UMI_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    customise_simdata_config(config)
    print(config)
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


def test_default_environment_vars(tmpdir):
    """
    Test that a workflow with no environment variables,
    :py:const:`riboviz.params.ENV_DIRS`, validates using the
    default value ``.`` for the environment variables.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    _ = create_test_dirs_envs_default(tmpdir)
    with open(test.SIMDATA_UMI_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    customise_simdata_config(config)
    print(config)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    envs = {}
    exit_code = nextflow.run_nextflow(test_config_file,
                                      envs=envs,
                                      validate_only=True,
                                      cwd=tmpdir)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


def test_relative_paths_template_config(tmpdir):
    """
    Test that a workflow with relative paths to input and output files
    validates.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    _ = create_test_dirs_envs_default(tmpdir)
    with open(test.SIMDATA_UMI_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    customise_simdata_config(config)
    # Make paths relative
    environment.update_config_with_env({}, config)
    print(config)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = nextflow.run_nextflow(test_config_file,
                                      validate_only=True,
                                      cwd=tmpdir)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


def test_relative_paths_default_config(tmpdir):
    """
    Test that a workflow with relative paths to input and output files
    validates.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    # TODO remove and remove cwd and run in current directory?
    _ = create_test_dirs(tmpdir)
    exit_code = nextflow.run_nextflow(test.SIMDATA_UMI_CONFIG,
                                      validate_only=True,
                                      cwd=tmpdir)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


def test_absolute_paths_template_config(tmpdir):
    """
    Test that a workflow with absolute paths to input and output files
    validates.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    organisms_dir, samples_dir, data_dir = create_test_dirs_envs(tmpdir)
    with open(test.SIMDATA_UMI_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    customise_simdata_config(config)
    # Make paths absolute
    envs = {params.ENV_RIBOVIZ_SAMPLES: samples_dir,
            params.ENV_RIBOVIZ_ORGANISMS: organisms_dir,
            params.ENV_RIBOVIZ_DATA: data_dir}
    environment.update_config_with_env(envs, config)
    print(config)
    test_config_file = tmpdir.join("test-config.yaml")
    with open(test_config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = nextflow.run_nextflow(test_config_file,
                                      validate_only=True)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code
