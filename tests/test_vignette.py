"""
Vignette regression test suite

The vignette (`pyscripts/prepRiboviz.py`) regression test suite
compares two directories, each assumed to have `index/` `tmp/` and
`output/` directories created by the vignette.

The tests can be run using pytest:

    pytest tests/test_vignette.py --expected=<DIR1> [--actual=<DIR2>]

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
import pytest
import riboviz
import riboviz.validation

# TODO replace iteration over files of the same format with parameterized tests, see https://docs.pytest.org/en/latest/parametrize.html and https://docs.pytest.org/en/latest/example/parametrize.html.

# TODO remove test_vignette sample function once tests are completed.
def test_vignette(command_option):
    (expected, actual) = command_option
    print(expected)
    print(actual)


def test_index(command_option):
    (expected, actual) = command_option
    expected_index = os.path.join(expected, "index")
    actual_index = os.path.join(actual, "index")
    for prefix in ["YAL_CDS_w_250", "yeast_rRNA"]:
        for index in range(1, 9):
            file_name = "%s.%d.ht2" % (prefix, index)
            print(file_name)
            riboviz.validation.compare(
                os.path.join(expected_index, file_name),
                os.path.join(actual_index, file_name))


@pytest.mark.skip(reason=".fq files take a long time to check")
def test_tmp_fq(command_option):
    (expected, actual) = command_option
    expected_tmp = os.path.join(expected, "tmp")
    actual_tmp = os.path.join(actual, "tmp")
    for prefix in ["WT3AT", "WTnone"]:
        for content in ["nonrRNA", "trim", "unaligned"]:
            file_name = "%s_%s.fq" % (prefix, content)
            print(file_name)
            riboviz.validation.compare(
                os.path.join(expected_tmp, file_name),
                os.path.join(actual_tmp, file_name))


def test_tmp_sam(command_option):
    (expected, actual) = command_option
    expected_tmp = os.path.join(expected, "tmp")
    actual_tmp = os.path.join(actual, "tmp")
    for prefix in ["WT3AT", "WTnone"]:
        for content in ["orf_map_clean", "orf_map", "rRNA_map"]:
            file_name = "%s_%s.sam" % (prefix, content)
            print(file_name)
            riboviz.validation.compare(
                os.path.join(expected_tmp, file_name),
                os.path.join(actual_tmp, file_name))


def test_output_bai(command_option):
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    for prefix in ["WT3AT", "WTnone"]:
        file_name = "%s.bam.bai" % prefix
        print(file_name)
        riboviz.validation.compare(
            os.path.join(expected_output, file_name),
            os.path.join(actual_output, file_name))


def test_output_bam(command_option):
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    for prefix in ["WT3AT", "WTnone"]:
        file_name = "%s.bam" % prefix
        print(file_name)
        riboviz.validation.compare(
            os.path.join(expected_output, file_name),
            os.path.join(actual_output, file_name))


def test_output_bedgraph(command_option):
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    for prefix in ["WT3AT", "WTnone"]:
        for content in ["minus", "plus"]:
            file_name = "%s_%s.bedgraph" % (prefix, content)
            print(file_name)
            riboviz.validation.compare(
                os.path.join(expected_output, file_name),
                os.path.join(actual_output, file_name))


def test_output_h5(command_option):
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    for prefix in ["WT3AT", "WTnone"]:
        file_name = "%s.h5" % prefix
        print(file_name)
        riboviz.validation.compare(
            os.path.join(expected_output, file_name),
            os.path.join(actual_output, file_name))


def test_output_tsv(command_option):
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    for prefix in ["WT3AT", "WTnone"]:
        for content in ["3nt_periodicity", "codon_ribodens", "pos_sp_nt_freq", "pos_sp_rpf_norm_reads", "read_lengths", "tpms"]:
            file_name = "%s_%s.tsv" % (prefix, content)
            print(file_name)
            riboviz.validation.compare(
                os.path.join(expected_output, file_name),
                os.path.join(actual_output, file_name))


def test_output_tpms_collated_tsv(command_option):
    (expected, actual) = command_option
    expected_output = os.path.join(expected, "output")
    actual_output = os.path.join(actual, "output")
    file_name = "TPMs_collated.tsv"
    riboviz.validation.compare(
        os.path.join(expected_output, file_name),
        os.path.join(actual_output, file_name))
