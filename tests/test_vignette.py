"""
vignette (pyscripts/prepRiboviz.py) regression test suite.

Tests compare a directory's data files against another directory's
data files.

The tests need to be run using pytest and with two command-line
options:

* '--expected=<DIRECTORY>': directory with expected data files.
* '--actual=<DIRECTORY>': directory to be validated against directory
  with expected data files.
"""
import os
import pytest
import riboviz
import riboviz.validation


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
