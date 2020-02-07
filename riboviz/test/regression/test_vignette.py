"""
Vignette regression test suite

The vignette regression test suite optionally runs prep_riboviz.py on
the vignette data, in `vignette/`, then compares the results, in
`vignette/`, to a directory of pre-calculated results.

The tests can be run using pytest:

    pytest riboviz/test/regression/test_vignette.py \
        --expected=<DIRECTORY> \
        [--skip-workflow]

where:

* `--expected=<DIRECTORY>`: directory with expected vignette
   files. This is assumed to have `index/` `tmp/` and `output/`
   directories.
* `--skip-workflow`: request that the `prep_riboviz.py` workflow not
   be run, instead use existing data files in `vignette/` for
   testing.

The directories are assumed to hold the following content:


    index/
      YAL_CDS_w_250.1.ht2
      YAL_CDS_w_250.2.ht2
      YAL_CDS_w_250.3.ht2
      YAL_CDS_w_250.4.ht2
      YAL_CDS_w_250.5.ht2
      YAL_CDS_w_250.6.ht2
      YAL_CDS_w_250.7.ht2
      YAL_CDS_w_250.8.ht2
      yeast_rRNA.1.ht2
      yeast_rRNA.2.ht2
      yeast_rRNA.3.ht2
      yeast_rRNA.4.ht2
      yeast_rRNA.5.ht2
      yeast_rRNA.6.ht2
      yeast_rRNA.7.ht2
      yeast_rRNA.8.ht2
    tmp/
      WT3AT/
        nonrRNA.fq
        orf_map_clean.sam
        orf_map.sam
        rRNA_map.sam
        trim.fq
        trim_5p_mismatch.tsv
        unaligned.fq
      WTnone/
        nonrRNA.fq
        orf_map_clean.sam
        orf_map.sam
        rRNA_map.sam
        trim.fq
        trim_5p_mismatch.tsv
        unaligned.fq
    output/
      TPMs_collated.tsv
      WT3AT/
        3nt_periodicity.pdf
        3nt_periodicity.tsv
        WT3AT.bam
        WT3AT.bam.bai
        3ntframe_bygene.tsv
        3ntframe_propbygene.pdf
        codon_ribodens.pdf
        codon_ribodens.tsv
        features.pdf
        WT3AT.h5
        minus.bedgraph
        plus.bedgraph
        pos_sp_nt_freq.tsv
        pos_sp_rpf_norm_reads.pdf
        pos_sp_rpf_norm_reads.tsv
        read_lengths.pdf
        read_lengths.tsv
        startcodon_ribogridbar.pdf
        startcodon_ribogrid.pdf
        tpms.tsv
      WTnone/
        3nt_periodicity.pdf
        3nt_periodicity.tsv
        WTnone.bam
        WTnone.bam.bai
        3ntframe_bygene.tsv
        3ntframe_propbygene.pdf
        codon_ribodens.pdf
        codon_ribodens.tsv
        features.pdf
        WTnone.h5
        minus.bedgraph
        plus.bedgraph
        pos_sp_nt_freq.tsv
        pos_sp_rpf_norm_reads.pdf
        pos_sp_rpf_norm_reads.tsv
        read_lengths.pdf
        read_lengths.tsv
        startcodon_ribogridbar.pdf
        startcodon_ribogrid.pdf
        tpms.tsv

See riboviz.validation.compare and riboviz.validation functions for
information on the nature of the comparisons for each type of file.
"""
import os
import shutil
import tempfile
import pytest
import pysam
import riboviz
from riboviz import file_names
from riboviz import h5
from riboviz import sam_bam
from riboviz import test
from riboviz import trim_5p_mismatch
from riboviz import validation
from riboviz import workflow
from riboviz import workflow_r
from riboviz.tools import prep_riboviz


INDEX_DIR = "index"
""" Name of index files directory, relative to "expected" directory """
TMP_DIR = "tmp"
""" Name of temporary files directory, relative to "expected" directory """
OUTPUT_DIR = "output"
""" Name of output files directory, relative to "expected" directory """


@pytest.fixture(scope="module")
def run_prep_riboviz(skip_workflow):
    """
    Fixture to optionally run RiboViz workflow on vignette data once
    per module (i.e. once only) then pass vignette directory to
    test methods.

    :param skip_workflow: Should workflow not be run?
    :type skip_workflow: bool
    """
    if not skip_workflow:
        exit_code = prep_riboviz.prep_riboviz(
            riboviz.R_SCRIPTS,
            test.VIGNETTE_CONFIG)
        assert exit_code == 0, \
            "prep_riboviz returned non-zero exit code %d" % exit_code


@pytest.fixture(scope="function")
def scratch_directory():
    """
    Create a temporary directory for any test files, and delete it
    after use.

    :return: directory
    :rtype: str or unicode
    """
    scratch_dir = tempfile.mkdtemp("scratch_test_vignette")
    yield scratch_dir
    shutil.rmtree(scratch_dir)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("prefix", test.INDEX_PREFIXES)
@pytest.mark.parametrize("index",
                         list(range(1, test.NUM_INDICES)))
def test_index(expected, prefix, index):
    """
    Test index/*.ht2 files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param prefix: file name prefix e.g. YAL_CDS_w_250
    :type prefix: str or unicode
    :param index: file name index e.g. 1
    :type index: int
    """
    file_name = "%s.%d.ht2" % (prefix, index)
    validation.compare(
        os.path.join(expected, INDEX_DIR, file_name),
        os.path.join(test.VIGNETTE_INDEX_DIR, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name", [
    file_names.NON_RRNA_FQ,
    file_names.ADAPTER_TRIM_FQ,
    file_names.UNALIGNED_FQ])
def test_sample_tmp_fq(expected, sample, file_name):
    """
    Test tmp/*.fq files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    :param file_name: file name e.g. nonrRNA.fq
    :type file_name: str or unicode
    """
    validation.compare(
        os.path.join(expected, TMP_DIR, sample, file_name),
        os.path.join(test.VIGNETTE_TMP_DIR, sample, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name", [
    file_names.ORF_MAP_CLEAN_SAM,
    file_names.ORF_MAP_SAM,
    file_names.RRNA_MAP_SAM])
def test_sample_tmp_sam(expected, scratch_directory, sample, file_name):
    """
    Test tmp/*.sam files for equality. The SAM files are sorted into
    temporary SAM files which are then compared for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param scratch_directory: scratch files directory for test files
    (pytest fixture defined in this module)
    :type scratch_directory: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    :param file_name: file name e.g. orf_map_clean.sam
    :type file_name: str or unicode
    """
    expected_file = os.path.join(
        expected, TMP_DIR, sample, file_name)
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
    validation.compare(expected_copy_file, actual_copy_file)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name", [
    trim_5p_mismatch.TRIM_5P_MISMATCH_FILE])
def test_sample_tmp_tsv(expected, sample, file_name):
    """
    Test tmp/*tsv files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    :param file_name: file name e.g. trim_5p_mismatch.tsv
    :type file_name: str or unicode
    """
    validation.compare(
        os.path.join(expected, TMP_DIR, sample, file_name),
        os.path.join(test.VIGNETTE_TMP_DIR, sample, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
def test_sample_output_bai(expected, sample):
    """
    Test output/*.bai files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    """
    file_name = sam_bam.BAM_BAI_FORMAT.format(
        sam_bam.BAM_FORMAT.format(sample))
    validation.compare(
        os.path.join(expected, OUTPUT_DIR, sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
def test_sample_output_bam(expected, sample):
    """
    Test output/*.bam files for equality. The BAM files are assumed to
    be sorted by leftmost coordinate position.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    """
    file_name = sam_bam.BAM_FORMAT.format(sample)
    validation.compare(
        os.path.join(expected, OUTPUT_DIR, sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name", [
    file_names.MINUS_BEDGRAPH,
    file_names.PLUS_BEDGRAPH])
def test_sample_output_bedgraph(expected, sample, file_name):
    """
    Test output/*.bedgraph files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    :param content: content e.g. minus
    :type content: str or unicode
    """
    validation.compare(
        os.path.join(expected, OUTPUT_DIR, sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
def test_sample_output_h5(expected, sample):
    """
    Test output/*.h5 files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    """
    file_name = h5.H5_FORMAT.format(sample)
    validation.compare(
        os.path.join(expected, OUTPUT_DIR, sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", test.VIGNETTE_SAMPLES)
@pytest.mark.parametrize("file_name",
                         [workflow_r.THREE_NT_PERIODICITY_TSV,
                          workflow_r.CODON_RIBODENS_TSV,
                          workflow_r.POS_SP_NT_FREQ_TSV,
                          workflow_r.POS_SP_RPF_NORM_READS_TSV,
                          workflow_r.READ_LENGTHS_TSV,
                          workflow_r.THREE_NT_FRAME_BY_GENE_TSV,
                          workflow_r.TPMS_TSV])
def test_sample_output_tsv(expected, sample, file_name):
    """
    Test output/*tsv files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    :param file_name: content e.g. 3nt_periodicity.tsv
    :type file_name: str or unicode
    """
    validation.compare(
        os.path.join(expected, OUTPUT_DIR, sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
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
def test_sample_output_pdf(expected, sample, file_name):
    """
    Test output/*pdf files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    :param file_name: content e.g. 3nt_periodicity.pdf
    :type file_name: str or unicode
    """
    validation.compare(
        os.path.join(expected, OUTPUT_DIR, sample, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, sample, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("file_name",
                         [workflow_r.TPMS_COLLATED_TSV])
def test_output_tsv(expected, file_name):
    """
    Test output/*.tsv files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param file_name: content e.g. TPMs_collated.tsv
    :type file_name: str or unicode
    """
    validation.compare(
        os.path.join(expected, OUTPUT_DIR, file_name),
        os.path.join(test.VIGNETTE_OUTPUT_DIR, file_name))
