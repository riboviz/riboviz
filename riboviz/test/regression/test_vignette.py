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
   files. This is assumed to	have `index/` `tmp/` and `output/`
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
import riboviz.process_utils
import riboviz.test
import riboviz.tools
import riboviz.validation
from riboviz import trim_5p_mismatch
from riboviz import workflow_r
from riboviz.tools import prep_riboviz


@pytest.fixture(scope="module")
def run_prep_riboviz(skip_workflow):
    """
    Fixture to optionally run RiboViz workflow on vignette data once
    per module (i.e. once only) then pass vignette directory to
    test methods.

    :param skip_workflow: Should workflow not be run?
    :type skip_workflow: bool
    """
    print(("Skip workflow:" + str(skip_workflow)))
    if not skip_workflow:
        exit_code = prep_riboviz.prep_riboviz(
            riboviz.R_SCRIPTS,
            riboviz.test.VIGNETTE_CONFIG)
        assert exit_code == 0, \
            "prep_riboviz returned non-zero exit code %d" % exit_code


@pytest.fixture(scope="function")
def tmp_directory():
    """
    Create a temporary directory for any test files, and delete it
    after use.

    :return: directory
    :rtype: str or unicode
    """
    tmpdir = tempfile.mkdtemp("tmp_test_vignette")
    yield tmpdir
    shutil.rmtree(tmpdir)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("index", list(range(1, 9)))
@pytest.mark.parametrize("prefix", ["YAL_CDS_w_250", "yeast_rRNA"])
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
    expected_index = os.path.join(expected, "index")
    actual_index = os.path.join(riboviz.test.VIGNETTE_DIR, "index")
    file_name = "%s.%d.ht2" % (prefix, index)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_index, file_name),
        os.path.join(actual_index, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
@pytest.mark.parametrize("file_name", [
    prep_riboviz.NON_RRNA_FQ,
    prep_riboviz.ADAPTER_TRIM_FQ,
    prep_riboviz.UNALIGNED_FQ])
def test_tmp_fq(expected, sample, file_name):
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
    expected_tmp = os.path.join(
        expected, "tmp", sample, file_name)
    actual_tmp = os.path.join(
        riboviz.test.VIGNETTE_DIR, "tmp", sample, file_name)
    print(file_name)
    riboviz.validation.compare(expected_tmp, actual_tmp)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
@pytest.mark.parametrize("file_name", [
    prep_riboviz.ORF_MAP_CLEAN_SAM,
    prep_riboviz.ORF_MAP_SAM,
    prep_riboviz.RRNA_MAP_SAM])
def test_tmp_sam(expected, tmp_directory, sample, file_name):
    """
    Test tmp/*.sam files for equality. The SAM files are sorted into
    temporary SAM files which are then compared for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param tmp_directory: temporary files directory for test files
    (pytest fixture defined in this module)
    :type tmp_directory: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    :param file_name: file name e.g. orf_map_clean.sam
    :type file_name: str or unicode
    """
    expected_tmp = os.path.join(
        expected, "tmp", sample, file_name)
    actual_tmp = os.path.join(
        riboviz.test.VIGNETTE_DIR, "tmp", sample, file_name)
    expected_tmp_copy_dir = os.path.join(tmp_directory, "expected")
    os.mkdir(expected_tmp_copy_dir)
    actual_tmp_copy_dir = os.path.join(tmp_directory, "actual")
    os.mkdir(actual_tmp_copy_dir)
    print(file_name)
    expected_tmp_copy_file = os.path.join(expected_tmp_copy_dir, file_name)
    actual_tmp_copy_file = os.path.join(actual_tmp_copy_dir, file_name)
    pysam.sort("-o",
               expected_tmp_copy_file,
               expected_tmp)
    pysam.sort("-o",
               actual_tmp_copy_file,
               actual_tmp)
    riboviz.validation.compare(expected_tmp_copy_file, actual_tmp_copy_file)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
@pytest.mark.parametrize("file_name", [
    trim_5p_mismatch.TRIM_5P_MISMATCH_FILE])
def test_tmp_tsv(expected, sample, file_name):
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
    expected_tmp = os.path.join(
        expected, "tmp", sample, file_name)
    actual_tmp = os.path.join(
        riboviz.test.VIGNETTE_DIR, "tmp", sample, file_name)
    print(file_name)
    riboviz.validation.compare(expected_tmp, actual_tmp)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
def test_output_bai(expected, sample):
    """
    Test output/*.bai files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    """
    file_name = prep_riboviz.BAM_BAI_FORMAT.format(sample)
    expected_output = os.path.join(
        expected, "output", sample, file_name)
    actual_output = os.path.join(
        riboviz.test.VIGNETTE_DIR, "output", sample, file_name)
    print(file_name)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
def test_output_bam(expected, sample):
    """
    Test output/*.bam files for equality. The BAM files are assumed to
    be sorted by leftmost coordinate position.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    """
    file_name = prep_riboviz.BAM_FORMAT.format(sample)
    expected_output = os.path.join(
        expected, "output", sample, file_name)
    actual_output = os.path.join(
        riboviz.test.VIGNETTE_DIR, "output", sample, file_name)
    print(file_name)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
@pytest.mark.parametrize("file_name", [
    prep_riboviz.MINUS_BEDGRAPH,
    prep_riboviz.PLUS_BEDGRAPH])
def test_output_bedgraph(expected, sample, file_name):
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
    expected_output = os.path.join(
        expected, "output", sample, file_name)
    actual_output = os.path.join(
        riboviz.test.VIGNETTE_DIR, "output", sample, file_name)
    print(file_name)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
def test_output_h5(expected, sample):
    """
    Test output/*.h5 files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param sample: sample name e.g. WT3AT
    :type sample: str or unicode
    """
    file_name = prep_riboviz.H5_FORMAT.format(sample)
    expected_output = os.path.join(
        expected, "output", sample, file_name)
    actual_output = os.path.join(
        riboviz.test.VIGNETTE_DIR, "output", sample, file_name)
    print(file_name)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
@pytest.mark.parametrize("file_name",
                         [workflow_r.THREE_NT_PERIODICITY_TSV,
                          workflow_r.CODON_RIBODENS_TSV,
                          workflow_r.POS_SP_NT_FREQ_TSV,
                          workflow_r.POS_SP_RPF_NORM_READS_TSV,
                          workflow_r.READ_LENGTHS_TSV,
                          workflow_r.THREE_NT_FRAME_BY_GENE_TSV,
                          workflow_r.TPMS_TSV])
def test_output_tsv(expected, sample, file_name):
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
    expected_output = os.path.join(
        expected, "output", sample, file_name)
    actual_output = os.path.join(
        riboviz.test.VIGNETTE_DIR, "output", sample, file_name)
    print(file_name)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("sample", ["WT3AT", "WTnone"])
@pytest.mark.parametrize("file_name",
                         [workflow_r.THREE_NT_PERIODICITY_PDF,
                          workflow_r.CODON_RIBODENS_PDF,
                          workflow_r.FEATURES_PDF,
                          workflow_r.POS_SP_RPF_NORM_READS_PDF,
                          workflow_r.READ_LENGTHS_PDF,
                          workflow_r.START_CODON_RIBOGRID_BAR_PDF,
                          workflow_r.START_CODON_RIBOGRID_PDF,
                          workflow_r.THREE_NT_FRAME_PROP_BY_GENE_PDF])
def test_output_pdf(expected, sample, file_name):
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
    expected_output = os.path.join(
        expected, "output", sample, file_name)
    actual_output = os.path.join(
        riboviz.test.VIGNETTE_DIR, "output", sample, file_name)
    print(file_name)
    riboviz.validation.compare(expected_output, actual_output)


@pytest.mark.usefixtures("run_prep_riboviz")
def test_output_tpms_collated_tsv(expected):
    """
    Test output/TPMs_collated.tsv files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    """
    file_name = workflow_r.TPMS_COLLATED_TSV
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(riboviz.test.VIGNETTE_DIR, "output")
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))
