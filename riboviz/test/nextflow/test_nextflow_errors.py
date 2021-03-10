"""
Nextflow error handling and exit code tests.

The test suite runs ``nextflow run prep_riboviz.nf``.
"""
import os.path
import shutil
import tempfile
import yaml
import pytest
import riboviz.test
from riboviz import hisat2
from riboviz import params
from riboviz.test.nextflow import run_nextflow


@pytest.fixture(scope="function")
def tmp_file():
    """
    Create a temporary file.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp", suffix=".yaml")
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


@pytest.fixture(scope="function")
def tmp_dir():
    """
    Create a temporary directory.

    :return: directory
    :rtype: str or unicode
    """
    tmp_dir = tempfile.mkdtemp("tmp")
    yield tmp_dir
    shutil.rmtree(tmp_dir)


def test_no_sample_multiplex_fq_files(tmp_file):
    """
    Test that missing :py:const:`riboviz.params.FQ_FILES` and
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES` returns a non-zero
    exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.FQ_FILES] = []
    config[params.MULTIPLEX_FQ_FILES] = []
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_both_sample_multiplex_fq_files(tmp_file):
    """
    Test that providing both :py:const:`riboviz.params.FQ_FILES` and
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES` returns a non-zero
    exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        vignette_config = yaml.load(f, yaml.SafeLoader)
    config[params.FQ_FILES] = vignette_config[params.FQ_FILES]
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_dir_in_not_found(tmp_file):
    """
    Test that a non-existent :py:const:`riboviz.params.INPUT_DIR`
    directory raises a non-zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.INPUT_DIR] = "foo"
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_fq_files_not_found(tmp_file):
    """
    Test that non-existent :py:const:`riboviz.params.FQ_FILES` files
    raise a non-zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.FQ_FILES] = {
        "foo1": "foo1.fq", "foo2": "foo2.fq"
    }
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_multiplex_fq_files_not_found(tmp_file):
    """
    Test that non-existent
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    files raise a non-zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.MULTIPLEX_FQ_FILES] = ["foo1.fq", "foo2.fq"]
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_no_sample_sheet(tmp_file):
    """
    Test that providing :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    but no :py:const:`riboviz.params.SAMPLE_SHEET` returns a non-zero
    exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[params.SAMPLE_SHEET]
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_sample_sheet_not_found(tmp_file):
    """
    Test that providing :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    and a non-existent :py:const:`riboviz.params.SAMPLE_SHEET` file
    returns a non-zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.SAMPLE_SHEET] = "foo.tsv"
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
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
def test_no_mandatory_parameter(tmp_file, parameter):
    """
    Test that not providing a mandatory parameter returns a non-zero
    exit code.

    This test also covers the case where if one of
    :py:const:`riboviz.params.T_RNA_FILE` or
    :py:const:`riboviz.params.CODON_POSITIONS_FILE` is provided then
    then the other must be too.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[parameter]
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
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
def test_file_not_found(tmp_file, parameter):
    """
    Test that providing a missing file for a file parameter returns a
    non-zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[parameter] = "foo"
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [params.FEATURES_FILE,
                          params.ASITE_DISP_LENGTH_FILE])
def test_optional_file_not_specified(tmp_file, parameter):
    """
    Test that validating a configuration with missing optional file
    parameters returns a zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[parameter]
    config[params.VALIDATE_ONLY] = True
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_t_rna_codon_positions_not_specified(tmp_file):
    """
    Test that validating a configuration with missing optional file
    parameters t_rna_file and codon_positions_file returns a zero exit
    code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[params.T_RNA_FILE]
    del config[params.CODON_POSITIONS_FILE]
    config[params.VALIDATE_ONLY] = True
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [params.FEATURES_FILE,
                          params.ASITE_DISP_LENGTH_FILE])
def test_optional_file_none(tmp_file, parameter):
    """
    Test that validating a configuration with optional file parameters
    set to 'none' returns a zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[parameter] = None
    config[params.VALIDATE_ONLY] = True
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_t_rna_codon_positions_none(tmp_file):
    """
    Test that validating a configuration with optional file
    parameters t_rna_file and codon_positions_file set to
    'none' returns a zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.T_RNA_FILE] = None
    config[params.CODON_POSITIONS_FILE] = None
    config[params.VALIDATE_ONLY] = True
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [params.T_RNA_FILE,
                          params.CODON_POSITIONS_FILE])
def test_t_rna_codon_positions_either_or(tmp_file, parameter):
    """
    Test that validating a configuration with co-dependent optional
    file parameters - t_rna_file and codon_positions_file - where
    one is defined and the other is not.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    :param parameter: Parameter
    :type parameter: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[parameter] = None
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_extract_umis_no_umi_regexp(tmp_file):
    """
    Test that if :py:const:`riboviz.params.EXTRACT_UMIS` is `TRUE`
    but no :py:const:`riboviz.params.UMI_REGEXP` returns a non-zero
    exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.SIMDATA_UMI_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    del config[params.UMI_REGEXP]
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


@pytest.mark.parametrize("parameter",
                         [(params.BUFFER, -1),
                          (params.COUNT_THRESHOLD, -1),
                          (params.NUM_PROCESSES, 0),
                          (params.MIN_READ_LENGTH, 0),
                          (params.MAX_READ_LENGTH, 0)], ids=str)
def test_invalid_value(tmp_file, parameter):
    """
    Test that providing invalid values for a parameter returns a
    non-zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    :param parameter: Parameter and invalid value
    :type parameter: tuple(str or unicode, int)
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    name, value = parameter
    config[name] = value
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_max_read_length_less_min(tmp_file):
    """
    Test that :py:const:`riboviz.params.MAX_READ_LENGTH` less than
    :py:const:`riboviz.params.MIN_READ_LENGTH` returns a non-zero exit
    code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.MIN_READ_LENGTH] = 10
    config[params.MAX_READ_LENGTH] = 9
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_build_indices_false_no_such_index_dir(tmp_file):
    """
    Test that :py:const:`riboviz.params.BUILD_INDICES` is false
    and no such index directory :py:const:`riboviz.params.INDEX_DIR`
    returns a non-zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.BUILD_INDICES] = False
    config[params.INDEX_DIR] = "NoSuchDirectory"
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_build_indices_false_no_such_orf_index_prefix(tmp_file, tmp_dir):
    """
    Test that :py:const:`riboviz.params.BUILD_INDICES` is false
    and no such files with prefixe
    :py:const:`riboviz.params.ORF_INDEX_PREFIX` returns a non-zero
    exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    :param tmp_dir: Path to temporary directory
    :type tmp_dir: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.BUILD_INDICES] = False
    config[params.INDEX_DIR] = tmp_dir
    # Create empty file with rRNA index prefix, so check for that
    # file will succeed.
    with open(os.path.join(tmp_dir,
                           hisat2.HT2_FORMAT.format(
                               config[params.RRNA_INDEX_PREFIX],
                               1)), 'w') as f:
        pass
    config[params.ORF_INDEX_PREFIX] = "NoSuchOrfPrefix"
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_build_indices_false_no_such_rrna_index_prefix(tmp_file, tmp_dir):
    """
    Test that :py:const:`riboviz.params.BUILD_INDICES` is false
    and no such files with prefixe
    :py:const:`riboviz.params.RRNA_INDEX_PREFIX` returns a non-zero
    exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    :param tmp_dir: Path to temporary directory
    :type tmp_dir: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.BUILD_INDICES] = False
    config[params.INDEX_DIR] = tmp_dir
    # Create empty file with ORF index prefix, so check for that
    # file will succeed.
    with open(os.path.join(tmp_dir,
                           hisat2.HT2_FORMAT.format(
                               config[params.ORF_INDEX_PREFIX],
                               1)), 'w') as f:
        pass
    config[params.RRNA_INDEX_PREFIX] = "NoSuchRrnaPrefix"
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code != 0, \
        "Unexpected exit code %d" % exit_code


def test_validate_skip_inputs_dir_in_not_found(tmp_file):
    """
    Test that a non-existent :py:const:`riboviz.params.INPUT_DIR`
    directory #in the presence of both
    :py:const:`riboviz.params.VALIDATE_ONLY` and
    :py:const:`riboviz.params.SKIP_INPUTS` returns a zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.INPUT_DIR] = "foo"
    config[params.VALIDATE_ONLY] = True
    config[params.SKIP_INPUTS] = True
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_validate_skip_inputs_fq_files_not_found(tmp_file):
    """
    Test that non-existent :py:const:`riboviz.params.FQ_FILES` files
    in the presence of both :py:const:`riboviz.params.VALIDATE_ONLY`
    and :py:const:`riboviz.params.SKIP_INPUTS` returns a zero exit
    code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.VIGNETTE_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.FQ_FILES] = {
        "foo1": "foo1.fq", "foo2": "foo2.fq"
    }
    config[params.VALIDATE_ONLY] = True
    config[params.SKIP_INPUTS] = True
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_validate_skip_inputs_multiplex_fq_files_not_found(tmp_file):
    """
    Test that non-existent
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES` in the presence of
    both :py:const:`riboviz.params.VALIDATE_ONLY` and
    :py:const:`riboviz.params.SKIP_INPUTS` returns a zero exit code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.MULTIPLEX_FQ_FILES] = ["foo1.fq", "foo2.fq"]
    config[params.VALIDATE_ONLY] = True
    config[params.SKIP_INPUTS] = True
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code


def test_validate_skip_inputs_sample_sheet_not_found(tmp_file):
    """
    Test that providing :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    and a non-existent :py:const:`riboviz.params.SAMPLE_SHEET` file
    in the presence of both :py:const:`riboviz.params.VALIDATE_ONLY`
    and :py:const:`riboviz.params.SKIP_INPUTS` returns a zero exit
    code.

    :param tmp_file: Path to temporary file, to write configuration to
    :type tmp_file: str or unicode
    """
    with open(riboviz.test.SIMDATA_MULTIPLEX_CONFIG, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    config[params.SAMPLE_SHEET] = "foo.tsv"
    config[params.VALIDATE_ONLY] = True
    config[params.SKIP_INPUTS] = True
    with open(tmp_file, 'w') as f:
        yaml.dump(config, f)
    exit_code = run_nextflow(tmp_file)
    assert exit_code == 0, \
        "Unexpected exit code %d" % exit_code
