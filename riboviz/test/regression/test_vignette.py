"""
:py:mod:`riboviz.tools.prep_riboviz` regression test suite.

The test suite accepts three custom command-line parameters:

* ``--expected=<DIRECTORY>``: Directory with expected data files,
  against which files in ``vignette/`` will be checked.
* ``--skip-workflow``: Workflow will not be run prior to checking data
  files.
* `--check-index-tmp`: Check index and temporary files (default is
  that only the output files are checked).

If ``--skip-workflow`` is provided then the module-level fixture
:py:func:`prep_riboviz_fixture` runs
:py:mod:`riboviz.tools.prep_riboviz` regression test suite using the
vignette configuration, ``vignette/vignette-config.yaml``,
(:py:const:`riboviz.test.VIGNETTE_CONFIG`).

The vignette output files (and the index and temporary files, if
``--check-index-tmp`` was provided) in ``vignette/`` are then compared
against those in the directory provided via the ``expected``
parameter.
"""
import os
import shutil
import tempfile
import pytest
import pysam
import riboviz
from riboviz import h5
from riboviz import hisat2
from riboviz import sam_bam
from riboviz import compare_files
from riboviz import workflow_files
from riboviz import workflow_r
from riboviz.tools import prep_riboviz
from riboviz import test


@pytest.fixture(scope="module")
def prep_riboviz_fixture(skip_workflow_fixture):
    """
    Run :py:mod:`riboviz.tools.prep_riboviz` if
    ``skip_workflow_fixture`` is not ``True``.

    :param skip_workflow_fixture: Should workflow not be run?
    :type skip_workflow_fixture: bool
    """
    if not skip_workflow_fixture:
        exit_code = prep_riboviz.prep_riboviz(
            riboviz.R_SCRIPTS,
            test.VIGNETTE_CONFIG)
        assert exit_code == 0, \
            "prep_riboviz returned non-zero exit code %d" % exit_code


@pytest.fixture(scope="function")
def scratch_directory():
    """
    Create a scratch directory.

    :return: directory
    :rtype: str or unicode
    """
    scratch_dir = tempfile.mkdtemp("tmp_scratch")
    yield scratch_dir
    shutil.rmtree(scratch_dir)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("prefix", test.INDEX_PREFIXES)
@pytest.mark.parametrize("index", list(range(1, test.NUM_INDICES)))
def test_index(expected_fixture, prefix, index):
    """
    Test HISAT2 index files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param prefix: File name prefix
    :type prefix: str or unicode
    :param index: File name index
    :type index: int
    """
    file_name = hisat2.HT2_FORMAT.format(prefix, index)
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_INDEX_DIR_NAME,
                     file_name),
        os.path.join(test.VIGNETTE_INDEX_DIR, file_name))


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name", [
    workflow_files.NON_RRNA_FQ,
    workflow_files.ADAPTER_TRIM_FQ,
    workflow_files.UNALIGNED_FQ])
def test_sample_tmp_fq(expected_fixture, sample, file_name):
    """
    Test sample-specific temporary FASTQ files for equality.
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_TMP_DIR_NAME,
                     sample, file_name),
        os.path.join(test.VIGNETTE_TMP_DIR, sample, file_name))


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name", [
    workflow_files.ORF_MAP_CLEAN_SAM,
    workflow_files.ORF_MAP_SAM,
    workflow_files.RRNA_MAP_SAM])
def test_sample_tmp_sam(expected_fixture, scratch_directory, sample,
                        file_name):
    """
    Test sample-specific temporary SAM files for equality. The SAM
    files are sorted into temporary SAM files which are then
    compared. See :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param scratch_directory: scratch files directory
    :type scratch_directory: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    expected_file = os.path.join(
        expected_fixture, test.VIGNETTE_TMP_DIR_NAME, sample, file_name)
    actual_file = os.path.join(
        test.VIGNETTE_TMP_DIR, sample, file_name)
    expected_copy_dir = os.path.join(scratch_directory, "expected")
    os.mkdir(expected_copy_dir)
    actual_copy_dir = os.path.join(scratch_directory, "actual")
    os.mkdir(actual_copy_dir)
    expected_copy_file = os.path.join(expected_copy_dir, file_name)
    actual_copy_file = os.path.join(actual_copy_dir, file_name)
    pysam.sort("-o", expected_copy_file, expected_file)
    pysam.sort("-o", actual_copy_file, actual_file)
    compare_files.compare_files(expected_copy_file, actual_copy_file)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name", [
    workflow_files.TRIM_5P_MISMATCH_TSV])
def test_sample_tmp_tsv(expected_fixture, sample, file_name):
    """
    Test sample-specific temporary TSV files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_TMP_DIR_NAME,
                     sample, file_name),
        os.path.join(test.VIGNETTE_TMP_DIR, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
def test_sample_output_bai(expected_fixture, sample):
    """
    Test sample-specific output BAI files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    """
    file_name = sam_bam.BAI_FORMAT.format(
        sam_bam.BAM_FORMAT.format(sample))
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_OUTPUT_DIR_NAME,
                     sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
def test_sample_output_bam(expected_fixture, sample):
    """
    Test sample-specific output BAM files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    The BAM files are assumed to be sorted by leftmost coordinate
    position.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    """
    file_name = sam_bam.BAM_FORMAT.format(sample)
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_OUTPUT_DIR_NAME,
                     sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name", [
    workflow_files.MINUS_BEDGRAPH,
    workflow_files.PLUS_BEDGRAPH])
def test_sample_output_bedgraph(expected_fixture, sample, file_name):
    """
    Test sample-specific output bedgraph files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_OUTPUT_DIR_NAME,
                     sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
def test_sample_output_h5(expected_fixture, sample):
    """
    Test sample-specific output H5 files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    """
    file_name = h5.H5_FORMAT.format(sample)
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_OUTPUT_DIR_NAME,
                     sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name",
                         [workflow_r.THREE_NT_PERIODICITY_TSV,
                          workflow_r.CODON_RIBODENS_TSV,
                          workflow_r.POS_SP_NT_FREQ_TSV,
                          workflow_r.POS_SP_RPF_NORM_READS_TSV,
                          workflow_r.READ_LENGTHS_TSV,
                          workflow_r.THREE_NT_FRAME_BY_GENE_TSV,
                          workflow_r.TPMS_TSV])
def test_sample_output_tsv(expected_fixture, sample, file_name):
    """
    Test sample-specific output TSV files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_OUTPUT_DIR_NAME,
                     sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name",
                         [workflow_r.THREE_NT_PERIODICITY_PDF,
                          workflow_r.CODON_RIBODENS_PDF,
                          workflow_r.FEATURES_PDF,
                          workflow_r.POS_SP_RPF_NORM_READS_PDF,
                          workflow_r.READ_LENGTHS_PDF,
                          workflow_r.START_CODON_RIBOGRID_BAR_PDF,
                          workflow_r.START_CODON_RIBOGRID_PDF,
                          workflow_r.THREE_NT_FRAME_PROP_BY_GENE_PDF])
def test_sample_output_pdf(expected_fixture, sample, file_name):
    """
    Test sample-specific output PDF files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param sample: sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_OUTPUT_DIR_NAME,
                     sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.TPMS_COLLATED_TSV,
                          workflow_files.READ_COUNTS_FILE])
def test_output_tsv(expected_fixture, file_name):
    """
    Test non-sample-specific output TSV files for equality. See
    :py:func:`riboviz.compare_files.compare_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_files.compare_files(
        os.path.join(expected_fixture, test.VIGNETTE_OUTPUT_DIR_NAME,
                     file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, file_name))
