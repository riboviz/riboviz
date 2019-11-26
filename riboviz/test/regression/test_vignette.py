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
        WT3AT_nonrRNA.fq
        WT3AT_orf_map_clean.sam
        WT3AT_orf_map.sam
        WT3AT_rRNA_map.sam
        WT3AT_trim.fq
        WT3AT_unaligned.fq
      WTnone/
        WTnone_nonrRNA.fq
        WTnone_orf_map_clean.sam
        WTnone_orf_map.sam
        WTnone_rRNA_map.sam
        WTnone_trim.fq
        WTnone_unaligned.fq
    output/
      TPMs_collated.tsv
      WT3AT/
        WT3AT_3nt_periodicity.pdf
        WT3AT_3nt_periodicity.tsv
        WT3AT.bam
        WT3AT.bam.bai
        WT3AT_3ntframe_bygene.tsv
        WT3AT_3ntframe_propbygene.pdf
        WT3AT_codon_ribodens.pdf
        WT3AT_codon_ribodens.tsv
        WT3AT_features.pdf
        WT3AT.h5
        WT3AT_minus.bedgraph
        WT3AT_plus.bedgraph
        WT3AT_pos_sp_nt_freq.tsv
        WT3AT_pos_sp_rpf_norm_reads.pdf
        WT3AT_pos_sp_rpf_norm_reads.tsv
        WT3AT_read_lengths.pdf
        WT3AT_read_lengths.tsv
        WT3AT_startcodon_ribogridbar.pdf
        WT3AT_startcodon_ribogrid.pdf
        WT3AT_tpms.tsv
      WTnone/
        WTnone_3nt_periodicity.pdf
        WTnone_3nt_periodicity.tsv
        WTnone.bam
        WTnone.bam.bai
        WTnone_3ntframe_bygene.tsv
        WTnone_3ntframe_propbygene.pdf
        WTnone_codon_ribodens.pdf
        WTnone_codon_ribodens.tsv
        WTnone_features.pdf
        WTnone.h5
        WTnone_minus.bedgraph
        WTnone_plus.bedgraph
        WTnone_pos_sp_nt_freq.tsv
        WTnone_pos_sp_rpf_norm_reads.pdf
        WTnone_pos_sp_rpf_norm_reads.tsv
        WTnone_read_lengths.pdf
        WTnone_read_lengths.tsv
        WTnone_startcodon_ribogridbar.pdf
        WTnone_startcodon_ribogrid.pdf
        WTnone_tpms.tsv

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
            riboviz.test.PY_SCRIPTS,
            riboviz.test.R_SCRIPTS,
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
@pytest.mark.parametrize("content", ["nonrRNA", "trim", "unaligned"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_tmp_fq(expected, prefix, content):
    """
    Test tmp/*.fq files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    :param content: content e.g. nonrRNA
    :type content: str or unicode
    """
    expected_tmp = os.path.join(expected, "tmp")
    actual_tmp = os.path.join(riboviz.test.VIGNETTE_DIR, "tmp")
    file_name = "%s_%s.fq" % (prefix, content)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_tmp, prefix, file_name),
        os.path.join(actual_tmp, prefix, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("content", ["orf_map_clean", "orf_map", "rRNA_map"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_tmp_sam(expected, tmp_directory, prefix, content):
    """
    Test tmp/*.sam files for equality. The SAM files are sorted into
    temporary SAM files which are then compared for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param tmp_directory: temporary files directory for test files
    (pytest fixture defined in this module)
    :type tmp_directory: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    :param content: content e.g. orf_map_clean
    :type content: str or unicode
    """
    expected_tmp = os.path.join(expected, "tmp", prefix)
    actual_tmp = os.path.join(riboviz.test.VIGNETTE_DIR, "tmp", prefix)
    expected_tmp_copy = os.path.join(tmp_directory, "expected")
    os.mkdir(expected_tmp_copy)
    actual_tmp_copy = os.path.join(tmp_directory, "actual")
    os.mkdir(actual_tmp_copy)
    file_name = "%s_%s.sam" % (prefix, content)
    print(file_name)
    expected_tmp_file = os.path.join(expected_tmp_copy, file_name)
    actual_tmp_file = os.path.join(actual_tmp_copy, file_name)
    pysam.sort("-o",
               expected_tmp_file,
               os.path.join(expected_tmp, file_name))
    pysam.sort("-o",
               actual_tmp_file,
               os.path.join(actual_tmp, file_name))
    riboviz.validation.compare(
        expected_tmp_file,
        actual_tmp_file)


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_bai(expected, prefix):
    """
    Test output/*.bai files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    """
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(riboviz.test.VIGNETTE_DIR, "output")
    file_name = "%s.bam.bai" % prefix
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, prefix, file_name),
        os.path.join(actual_output, prefix, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_bam(expected, prefix):
    """
    Test output/*.bam files for equality. The BAM files are assumed to
    be sorted by leftmost coordinate position.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    """
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(riboviz.test.VIGNETTE_DIR, "output")
    file_name = "%s.bam" % prefix
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, prefix, file_name),
        os.path.join(actual_output, prefix, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("content", ["minus", "plus"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_bedgraph(expected, prefix, content):
    """
    Test tmp/*.bedgraph files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    :param content: content e.g. minus
    :type content: str or unicode
    """
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(riboviz.test.VIGNETTE_DIR, "output")
    file_name = "%s_%s.bedgraph" % (prefix, content)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, prefix, file_name),
        os.path.join(actual_output, prefix, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_h5(expected, prefix):
    """
    Test output/*.h5 files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    """
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(riboviz.test.VIGNETTE_DIR, "output")
    file_name = "%s.h5" % prefix
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, prefix, file_name),
        os.path.join(actual_output, prefix, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("content",
                         ["3nt_periodicity",
                          "codon_ribodens",
                          "pos_sp_nt_freq",
                          "pos_sp_rpf_norm_reads",
                          "read_lengths",
                          "3ntframe_bygene",
                          "tpms"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_tsv(expected, prefix, content):
    """
    Test output/*tsv files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :type prefix: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :param content: content e.g. 3nt_periodicity
    :type content: str or unicode
    """
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(riboviz.test.VIGNETTE_DIR, "output")
    file_name = "%s_%s.tsv" % (prefix, content)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, prefix, file_name),
        os.path.join(actual_output, prefix, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
@pytest.mark.parametrize("content",
                         ["3nt_periodicity",
                          "codon_ribodens",
                          "features",
                          "pos_sp_rpf_norm_reads",
                          "read_lengths",
                          "startcodon_ribogridbar",
                          "startcodon_ribogrid",
                          "3ntframe_propbygene"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_pdf(expected, prefix, content):
    """
    Test output/*pdf files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    :param content: content e.g. 3nt_periodicity
    :type content: str or unicode
    """
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(riboviz.test.VIGNETTE_DIR, "output")
    file_name = "%s_%s.pdf" % (prefix, content)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, prefix, file_name),
        os.path.join(actual_output, prefix, file_name))


@pytest.mark.usefixtures("run_prep_riboviz")
def test_output_tpms_collated_tsv(expected):
    """
    Test output/TPMs_collated.tsv files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    """
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(riboviz.test.VIGNETTE_DIR, "output")
    file_name = "TPMs_collated.tsv"
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))
