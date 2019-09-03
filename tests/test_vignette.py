"""
Vignette regression test suite

The vignette regression test suite compares two directories, each
assumed to have `index/` `tmp/` and `output/` directories created by
the vignette.

The tests can be run using pytest:

    pytest tests/test_vignette.py --expected=<DIRECTORY> \
                                  [--actual=<DIRECTORY>]

where:

* `--expected=<DIRECTORY>`: directory with expected vignette files.
* `--actual=<DIRECTORY>`: directory to be validated against directory
  with expected vignette files. Default: `vignette/`

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
      WT3AT_nonrRNA.fq
      WT3AT_orf_map_clean.sam
      WT3AT_orf_map.sam
      WT3AT_rRNA_map.sam
      WT3AT_trim.fq
      WT3AT_unaligned.fq
      WTnone_nonrRNA.fq
      WTnone_orf_map_clean.sam
      WTnone_orf_map.sam
      WTnone_rRNA_map.sam
      WTnone_trim.fq
      WTnone_unaligned.fq
    output/
      TPMs_collated.tsv
      WT3AT_3nt_periodicity.pdf
      WT3AT_3nt_periodicity.tsv
      WT3AT.bam
      WT3AT.bam.bai
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
      WT3AT_tpms.tsv
      WTnone_3nt_periodicity.pdf
      WTnone_3nt_periodicity.tsv
      WTnone.bam
      WTnone.bam.bai
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
      WTnone_tpms.tsv
"""
import os
import shutil
import tempfile
import pytest
import pysam
import riboviz
import riboviz.validation


@pytest.fixture(scope="function")
def tmp_directory():
    """
    Create a temporary directory for any test files, and delete it
    after use.
    """
    tmpdir = tempfile.mkdtemp("tmp_test_vignette")
    yield tmpdir
    shutil.rmtree(tmpdir)


@pytest.mark.parametrize("index", range(1, 9))
@pytest.mark.parametrize("prefix", ["YAL_CDS_w_250", "yeast_rRNA"])
def test_index(command_option, prefix, index):
    """
    Test index/*.ht2 files for equality.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    :param prefix: file name prefix e.g. YAL_CDS_w_250
    :type prefix: str or unicode
    :param index: file name index e.g. 1
    :type index: int
    """
    (expected, actual) = command_option
    expected_index = os.path.join(expected, "index")
    actual_index = os.path.join(actual, "index")
    file_name = "%s.%d.ht2" % (prefix, index)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_index, file_name),
        os.path.join(actual_index, file_name))


@pytest.mark.parametrize("content", ["nonrRNA", "trim", "unaligned"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_tmp_fq(command_option, prefix, content):
    """
    Test tmp/*.fq files for equality.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    :param content: content e.g. nonrRNA
    :type content: str or unicode
    """
    (expected, actual) = command_option
    expected_tmp = os.path.join(expected, "tmp")
    actual_tmp = os.path.join(actual, "tmp")
    file_name = "%s_%s.fq" % (prefix, content)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_tmp, file_name),
        os.path.join(actual_tmp, file_name))


@pytest.mark.parametrize("content", ["orf_map_clean", "orf_map", "rRNA_map"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_tmp_sam(command_option, tmp_directory, prefix, content):
    """
    Test tmp/*.sam files for equality. The SAM files are sorted into
    temporary SAM files which are then compared for equality.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    :param tmp_directory: directory for temporary SAM files
    (pytest fixture)
    :type tmp_directory: str or unicode
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    :param content: content e.g. orf_map_clean
    :type content: str or unicode
    """
    (expected, actual) = command_option
    expected_tmp = os.path.join(expected, "tmp")
    actual_tmp = os.path.join(actual, "tmp")
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


@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_bai(command_option, prefix):
    """
    Test output/*.bai files for equality.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    """
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    file_name = "%s.bam.bai" % prefix
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))


@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_bam(command_option, prefix):
    """
    Test output/*.bam files for equality. The BAM files are assumed to
    be sorted by leftmost coordinate position.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    """
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    file_name = "%s.bam" % prefix
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))


@pytest.mark.parametrize("content", ["minus", "plus"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_bedgraph(command_option, prefix, content):
    """
    Test tmp/*.bedgraph files for equality.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    :param content: content e.g. minus
    :type content: str or unicode
    """
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    file_name = "%s_%s.bedgraph" % (prefix, content)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))


@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_h5(command_option, prefix):
    """
    Test output/*.h5 files for equality.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    """
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    file_name = "%s.h5" % prefix
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))


@pytest.mark.parametrize("content",
                         ["3nt_periodicity",
                          "codon_ribodens",
                          "pos_sp_nt_freq",
                          "pos_sp_rpf_norm_reads",
                          "read_lengths",
                          "tpms"])
@pytest.mark.parametrize("prefix", ["WT3AT", "WTnone"])
def test_output_tsv(command_option, prefix, content):
    """
    Test output/*tsv files for equality.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    :param prefix: file name prefix e.g. WT3AT
    :type prefix: str or unicode
    :param content: content e.g. 3nt_periodicity
    :type content: str or unicode
    """
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    file_name = "%s_%s.tsv" % (prefix, content)
    print(file_name)
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))


def test_output_tpms_collated_tsv(command_option):
    """
    Test output/TPMs_collated.tsv files for equality.

    :param command_option: expected directory, actual directory
    (pytest fixture defined in conftest.py)
    :type command_option: tuple(str or unicode, str or unicode)
    """
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    file_name = "TPMs_collated.tsv"
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))
