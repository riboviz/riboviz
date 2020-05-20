"""
Nextflow error handling and exit code tests.

The test suite runs ``nextflow run prep_riboviz.nf``.
"""
import os.path
import tempfile
import yaml
import pytest
import riboviz.test
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


def test_no_sample_multiplex_fq_files(tmp_file):
    """
    Test that missing :py:const:`riboviz.params.FQ_FILES` and
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES` raises a non-zero
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
    :py:const:`riboviz.params.MULTIPLEX_FQ_FILES` raises a non-zero
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


def test_multiplex_fq_files_no_sample_sheet(tmp_file):
    """
    Test that providing :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    but no :py:const:`riboviz.params.SAMPLE_SHEET` raises a non-zero
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


def test_multiplex_fq_files_sample_sheet_not_found(tmp_file):
    """
    Test that providing :py:const:`riboviz.params.MULTIPLEX_FQ_FILES`
    and a non-existent :py:const:`riboviz.params.SAMPLE_SHEET` file
    raises a non-zero exit code.

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
    Test that not providing a mandatory parameter raises a non-zero
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
    Test that providing a missing file for a file parameter raises a
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


def test_extract_umis_no_umi_regexp(tmp_file):
    """
    Test that if :py:const:`riboviz.params.EXTRACT_UMIS` is `TRUE`
    but no :py:const:`riboviz.params.UMI_REGEXP` raises a non-zero
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
                          (params.MAX_READ_LENGTH, 0)])
def test_invalid_value(tmp_file, parameter):
    """
    Test that providing invalid values for a parameter raises a
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
    :py:const:`riboviz.params.MIN_READ_LENGTH` raises a non-zero exit
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
