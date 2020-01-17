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
from riboviz.tools import prep_riboviz
from riboviz.trim_5p_mismatch import TRIM_5P_MISMATCH_FILE


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
    file_name = "%s.fq" % content
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
    file_name = "%s.sam" % content
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
@pytest.mark.parametrize("content", [TRIM_5P_MISMATCH_FILE])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_tmp_tsv(expected, prefix, content):
    """
    Test tmp/*tsv files for equality.

    :param expected: expected directory
    (pytest fixture defined in conftest.py)
    :type expected: str or unicode
    :type prefix: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :param content: content e.g. 3nt_periodicity
    :type content: str or unicode
    """
    expected_tmp = os.path.join(expected, "tmp", prefix, content)
    actual_tmp = os.path.join(riboviz.test.VIGNETTE_DIR,
                              "tmp", prefix, content)
    file_name = "%s" % content
    print(file_name)
    riboviz.validation.compare(expected_tmp, actual_tmp)


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
    file_name = "%s.bedgraph" % content
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
    file_name = "%s.tsv" % content
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
    file_name = "%s.pdf" % content
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
