"""
Nextflow validation, error handling and exit code tests.

The test suite runs ``nextflow run prep_riboviz.nf``.
"""
import yaml
import pytest
import riboviz.test
from riboviz import hisat2
from riboviz import params
from riboviz.test.nextflow import run_nextflow


def test_no_sample_multiplex_fq_files(tmpdir):
    """
    Test that missing :py:const:`riboviz.params.FQ_FILES` and
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES` returns a non-zero
    exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.FQ_FILES] = []
    config[params.MULTIPLEX_FQ_FILES] = []
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_both_sample_multiplex_fq_files(tmpdir):
    """
    Test that providing both :py:const:`riboviz.params.FQ_FILES` and
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES` returns a non-zero
    exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        vignette_config = yaml.load(f, yaml.SafeLoader)
    config[params.FQ_FILES] = vignette_config[params.FQ_FILES]
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_dir_in_not_found(tmpdir):
    """
    Test that a non-existent :py:const:`riboviz.params.INPUT_DIR`
    directory raises a non-zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.INPUT_DIR] = "foo"
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_fq_files_not_found(tmpdir):
    """
    Test that non-existent :py:const:`riboviz.params.FQ_FILES` files
    raise a non-zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.FQ_FILES] = {
        "foo1": "foo1.fq", "foo2": "foo2.fq"
    }
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_multiplex_fq_files_not_found(tmpdir):
    """
    Test that non-existent
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    files raise a non-zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.MULTIPLEX_FQ_FILES] = ["foo1.fq", "foo2.fq"]
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_no_sample_sheet(tmpdir):
    """
    Test that providing :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    but no :py:const:`riboviz.params.SAMPLE_SHEET` returns a non-zero
    exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[params.SAMPLE_SHEET]
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_sample_sheet_not_found(tmpdir):
    """
    Test that providing :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    and a non-existent :py:const:`riboviz.params.SAMPLE_SHEET` file
    returns a non-zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.SAMPLE_SHEET] = "foo.tsv"
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [params.INPUT_DIR,
                          params.RRNA_FASTA_FILE,
                          params.ORF_FASTA_FILE,
                          params.ORF_GFF_FILE,
                          params.ADAPTERS,
                          params.ORF_INDEX_PREFIX,
                          params.RRNA_INDEX_PREFIX,
                          params.T_RNA_FILE,
                          params.CODON_POSITIONS_FILE])
def test_no_mandatory_parameter(tmpdir, parameter):
    """
    Test that not providing a mandatory parameter returns a non-zero
    exit code.

    This test also covers the case where if one of
    :py:const:`riboviz.params.T_RNA_FILE` or
    :py:const:`riboviz.params.CODON_POSITIONS_FILE` is provided then
    then the other must be too.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[parameter]
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [params.INPUT_DIR,
                          params.RRNA_FASTA_FILE,
                          params.ORF_FASTA_FILE,
                          params.ORF_GFF_FILE,
                          params.FEATURES_FILE,
                          params.T_RNA_FILE,
                          params.CODON_POSITIONS_FILE,
                          params.ASITE_DISP_LENGTH_FILE])
def test_file_not_found(tmpdir, parameter):
    """
    Test that providing a missing file for a file parameter returns a
    non-zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[parameter] = "foo"
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [params.FEATURES_FILE,
                          params.ASITE_DISP_LENGTH_FILE])
def test_optional_file_not_specified(tmpdir, parameter):
    """
    Test that validating a configuration with missing optional file
    parameters returns a zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[parameter]
    config[params.VALIDATE_ONLY] = True
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_t_rna_codon_positions_not_specified(tmpdir):
    """
    Test that validating a configuration with missing optional file
    parameters t_rna_file and codon_positions_file returns a zero exit
    code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[params.T_RNA_FILE]
    del config[params.CODON_POSITIONS_FILE]
    config[params.VALIDATE_ONLY] = True
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [params.FEATURES_FILE,
                          params.ASITE_DISP_LENGTH_FILE])
def test_optional_file_none(tmpdir, parameter):
    """
    Test that validating a configuration with optional file parameters
    set to 'none' returns a zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[parameter] = None
    config[params.VALIDATE_ONLY] = True
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_t_rna_codon_positions_none(tmpdir):
    """
    Test that validating a configuration with optional file
    parameters t_rna_file and codon_positions_file set to
    'none' returns a zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.T_RNA_FILE] = None
    config[params.CODON_POSITIONS_FILE] = None
    config[params.VALIDATE_ONLY] = True
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [params.T_RNA_FILE,
                          params.CODON_POSITIONS_FILE])
def test_t_rna_codon_positions_either_or(tmpdir, parameter):
    """
    Test that validating a configuration with co-dependent optional
    file parameters - t_rna_file and codon_positions_file - where
    one is defined and the other is not.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[parameter] = None
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_extract_umis_no_umi_regexp(tmpdir):
    """
    Test that if :py:const:`riboviz.params.EXTRACT_UMIS` is `TRUE`
    but no :py:const:`riboviz.params.UMI_REGEXP` returns a non-zero
    exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.SIMDATA_UMI_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[params.UMI_REGEXP]
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [(params.BUFFER, -1),
                          (params.COUNT_THRESHOLD, -1),
                          (params.NUM_PROCESSES, 0),
                          (params.MIN_READ_LENGTH, 0),
                          (params.MAX_READ_LENGTH, 0)], ids=str)
def test_invalid_value(tmpdir, parameter):
    """
    Test that providing invalid values for a parameter returns a
    non-zero exit code.

    :param parameter: Parameter and invalid value
    :type parameter: tuple(str or unicode, int)
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    name, value = parameter
    config[name] = value
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_max_read_length_less_min(tmpdir):
    """
    Test that :py:const:`riboviz.params.MAX_READ_LENGTH` less than
    :py:const:`riboviz.params.MIN_READ_LENGTH` returns a non-zero exit
    code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.MIN_READ_LENGTH] = 10
    config[params.MAX_READ_LENGTH] = 9
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_build_indices_false_no_such_index_dir(tmpdir):
    """
    Test that :py:const:`riboviz.params.BUILD_INDICES` is false
    and no such index directory :py:const:`riboviz.params.INDEX_DIR`
    returns a non-zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.BUILD_INDICES] = False
    config[params.INDEX_DIR] = "NoSuchDirectory"
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_build_indices_false_no_such_orf_index_prefix(tmpdir):
    """
    Test that :py:const:`riboviz.params.BUILD_INDICES` is false
    and no such files with prefixe
    :py:const:`riboviz.params.ORF_INDEX_PREFIX` returns a non-zero
    exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    index_dir = tmpdir.mkdir("index")
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.BUILD_INDICES] = False
    # Cast index_dir to string. If not done, then, when config is saved
    # the value will be
    # !!python/object:py._path.local.LocalPath
    #   strpath: /tmp/.../index
    # Casting to string gives the value /tmp/.../index/ as desired.
    config[params.INDEX_DIR] = str(index_dir)
    # Create empty file for rRNA index prefix, so check for that
    # file will succeed.
    index_file = index_dir.join(hisat2.HT2_FORMAT.format(
        config[params.RRNA_INDEX_PREFIX], 1))
    with open(index_file, 'w') as f:
        pass
    config[params.ORF_INDEX_PREFIX] = "NoSuchOrfPrefix"
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_build_indices_false_no_such_rrna_index_prefix(tmpdir):
    """
    Test that :py:const:`riboviz.params.BUILD_INDICES` is false
    and no such files with prefixe
    :py:const:`riboviz.params.RRNA_INDEX_PREFIX` returns a non-zero
    exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    index_dir = tmpdir.mkdir("index")
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.BUILD_INDICES] = False
    # Cast index_dir to string. If not done, then, when config is saved
    # the value will be
    # !!python/object:py._path.local.LocalPath
    #   strpath: /tmp/.../index
    # Casting to string gives the value /tmp/.../index/ as desired.
    config[params.INDEX_DIR] = str(index_dir)
    # Create empty file for ORF index prefix, so check for that
    # file will succeed.
    index_file = index_dir.join(hisat2.HT2_FORMAT.format(
        config[params.ORF_INDEX_PREFIX], 1))
    with open(index_file, 'w') as f:
        pass
    config[params.RRNA_INDEX_PREFIX] = "NoSuchRrnaPrefix"
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_validate_skip_inputs_dir_in_not_found(tmpdir):
    """
    Test that a non-existent :py:const:`riboviz.params.INPUT_DIR`
    directory #in the presence of both
    :py:const:`riboviz.params.VALIDATE_ONLY` and
    :py:const:`riboviz.params.SKIP_INPUTS` returns a zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.INPUT_DIR] = "foo"
    config[params.VALIDATE_ONLY] = True
    config[params.SKIP_INPUTS] = True
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_validate_skip_inputs_fq_files_not_found(tmpdir):
    """
    Test that non-existent :py:const:`riboviz.params.FQ_FILES` files
    in the presence of both :py:const:`riboviz.params.VALIDATE_ONLY`
    and :py:const:`riboviz.params.SKIP_INPUTS` returns a zero exit
    code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.FQ_FILES] = {
        "foo1": "foo1.fq", "foo2": "foo2.fq"
    }
    config[params.VALIDATE_ONLY] = True
    config[params.SKIP_INPUTS] = True
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_validate_skip_inputs_multiplex_fq_files_not_found(tmpdir):
    """
    Test that non-existent
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES` in the presence of
    both :py:const:`riboviz.params.VALIDATE_ONLY` and
    :py:const:`riboviz.params.SKIP_INPUTS` returns a zero exit code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.MULTIPLEX_FQ_FILES] = ["foo1.fq", "foo2.fq"]
    config[params.VALIDATE_ONLY] = True
    config[params.SKIP_INPUTS] = True
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_validate_skip_inputs_sample_sheet_not_found(tmpdir):
    """
    Test that providing :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    and a non-existent :py:const:`riboviz.params.SAMPLE_SHEET` file
    in the presence of both :py:const:`riboviz.params.VALIDATE_ONLY`
    and :py:const:`riboviz.params.SKIP_INPUTS` returns a zero exit
    code.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.SAMPLE_SHEET] = "foo.tsv"
    config[params.VALIDATE_ONLY] = True
    config[params.SKIP_INPUTS] = True
    config_file = tmpdir.join("config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(config_file, validate_only=True)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_environment_vars(tmpdir):
    """
    Test that a workflow with environment variables,
    :py:const:`riboviz.params.ENV_DIRS`, validates.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir py._path.local.LocalPath
    """
    envs = {env: tmpdir for env in params.ENV_DIRS}
    exit_code = run_nextflow(riboviz.test.VIGNETTE_CONFIG,
                             envs=envs,
                             validate_only=True)
    assert exit_code == 0, "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("env", params.ENV_DIRS)
def test_environment_var_not_found(env):
    """
    Test that a workflow with an environment variable pointing
    to a non-existent path raises a non-zero exit code.

    :param env: Environment variable
    :type env: str or unicode
    """
    exit_code = run_nextflow(riboviz.test.VIGNETTE_CONFIG,
                             envs={env: "noSuchDirectory"},
                             validate_only=True)
    assert exit_code != 0, "Unexpected exit code %d" % exit_code
