"""
:py:mod:`riboviz.test_fasta_gff` tests.
"""
import os
import tempfile
import pytest
import pandas as pd
from riboviz import fasta_gff
from riboviz.test import data


TEST_FASTA_FILE = os.path.join(os.path.dirname(data.__file__),
                               "test_fasta_gff_data.fasta")
""" Test FASTA file in :py:mod:`riboviz.test.data`. """
TEST_GFF_FILE = os.path.join(os.path.dirname(data.__file__),
                             "test_fasta_gff_data.gff")
""" Test GFF file in :py:mod:`riboviz.test.data`. """
TEST_GFF_NO_CDS_FILE = os.path.join(os.path.dirname(data.__file__),
                                    "test_fasta_gff_data_no_cds.gff")
""" Test GFF file in :py:mod:`riboviz.test.data` with no CDS. """
TEST_CDS_CODONS = {
    "YAL001C_CDS": ["ATG", "GCC", "CAC", "TGT", "TAA"],
    "YAL002C_CDS": ["ATG", "GTA", "TCA", "GGA", "TAG"],
    "YAL004CSingleCodonCDS_CDS": ["ATG", "AGA", "TGA"],
    "YAL005CMultiCDS_CDS_1": ["ATG", "AGA", "TGA"],
    "YAL005CMultiCDS_CDS_2": ["ATG", "GAT", "TAC", "TAG"],
    "YAL005CMultiCDS_CDS_3": ["ATG", "CCA", "ATT", "TGA"],
    "YAL006CEmptyCDS_CDS": ["ATG", "TGA"],
    "YAL008CNoNameAttr_mRNA_CDS": ["ATG", "GCC", "CAC", "TGT", "TAA"],
    "YAL009CMultiDuplicateCDS_CDS": ["ATG", "AGA", "TGA"],
    "YAL009CMultiDuplicateCDS_CDS.1": ["ATG", "GAT", "TAC", "TAG"],
    "YAL009CMultiDuplicateCDS_CDS.2": ["ATG", "CCA", "ATT", "TGA"]
}
"""
Expected codons for CDS in FASTA file (:py:const:`TEST_FASTA_FILE`)
as defined in GFF file (:py:const:`TEST_GFF_FILE`).
"""


class MockFeature:
    """
    Mock of gffutils.feature.Feature class supporting ``seqid``, ``seq``
    and ``attributes`` attributes.
    """

    def __init__(self, seqid, seq="", attributes={}):
        """
        Constructor.

        :param self: Object reference
        :type self: MockFeature
        :param seqid: Sequence ID
        :type seqid: str or unicode
        :param seq: Sequence
        :type seq: str or unicode
        :param attributes: Attributes
        :type attributes: dict
        """
        self.seqid = seqid
        self.seq = seq
        self.attributes = attributes

    def sequence(self, fasta):
        """
        Mock of.Feature.sequence function, which returns the value of
        the ``seq`` attribute.

        :param self: Object reference
        :type self: MockFeature
        :param fasta: FASTA file (ignored)
        :type fasta: str or unicode
        :return: Sequence
        :rtype: str or unicode
        """
        return self.seq


@pytest.fixture(scope="function")
def tmp_file():
    """
    Create a temporary file.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp")
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


@pytest.mark.parametrize("seq_codons", [
    ("", []),
    ("A", ["A"]),
    ("GA", ["GA"]),
    ("GAT", ["GAT"]),
    ("GATT", ["GAT", "T"]),
    ("ATGAAATAA", ["ATG", "AAA", "TAA"]),
    ("ATGGGGCCCTAG", ["ATG", "GGG", "CCC", "TAG"])])
def test_sequence_to_codons(seq_codons):
    """
    Test :py:func:`riboviz.fasta_gff.sequence_to_codons`
    with valid sequences.

    :param seq_codons: sequence and expected codons
    :type seq_codons: tuple(str or unicode, list(str or unicode))
    """
    sequence, codons = seq_codons
    assert fasta_gff.sequence_to_codons(sequence) == codons,\
        "Unexpected codons"


@pytest.mark.parametrize("features", [
    ("SeqID_mRNA", {}, None),
    ("SeqID_mRNA", {"Name": []}, None),
    ("SeqID_mRNA", {"Name": ["SeqID_CDS"]}, "SeqID_CDS"),
    ("SeqID_mRNA", {"ID": []}, None),
    ("SeqID_mRNA", {"ID": ["SeqID_CDS"]}, "SeqID_CDS")])
def test_get_feature_name(features):
    """
    Test :py:func:`riboviz.fasta_gff.get_feature_name`.

    :param features: sequence ID, attributes, expected feature name
    :type features: tuple(str or unicode, str or unicode, str or unicode)
    """
    seq_id, attributes, feature_name = features
    assert fasta_gff.get_feature_name(
        MockFeature(seq_id, "", attributes)) == feature_name


def test_get_cds_from_fasta():
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_from_fasta` returns
    a sequence.

    :py:class:`MockFeature` is used to mock
    ``gffutils.feature.Feature`` to avoid the need to use a FASTA
    file for this test.
    """
    expected_sequence = "ATGAAATAA"
    feature = MockFeature("SeqID", expected_sequence)
    # Any file name can be used as MockFeature ignores it and returns
    # the sequence given above.
    sequence = fasta_gff.get_cds_from_fasta(feature, "test.fasta")
    assert sequence == expected_sequence


def test_get_cds_from_fasta_bad_length():
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_from_fasta` raises an
    error if given a sequence whose length is not divisible by 3.

    :py:class:`MockFeature` is used to mock
    ``gffutils.feature.Feature`` to avoid the need to use a FASTA
    file for this test.
    """
    feature = MockFeature("SeqID", "ATGATAA")
    with pytest.raises(AssertionError):
        # Any file name can be used as MockFeature ignores it
        # and returns the sequence given above.
        fasta_gff.get_cds_from_fasta(feature, "test.fasta")


def test_get_cds_codons_from_fasta_empty(tmp_file):
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_from_fasta`
    with an empty FASTA file and GFF file
    (:py:const:`TEST_GFF_FILE`).

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    # Use tmp_file as both empty FASTA input file.
    cds_codons = fasta_gff.get_cds_codons_from_fasta(
        tmp_file,
        TEST_GFF_FILE)
    assert cds_codons == {}


def test_get_cds_codons_from_fasta():
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_FILE`) and GFF file
    (:py:const:`TEST_GFF_FILE`).
    """
    cds_codons = fasta_gff.get_cds_codons_from_fasta(
        TEST_FASTA_FILE,
        TEST_GFF_FILE)
    assert cds_codons == TEST_CDS_CODONS


def test_get_cds_codons_from_fasta_exclude_stop_codon():
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_FILE`) and GFF file
    (:py:const:`TEST_GFF_FILE`) where stop codons are to be
    excluded from the codons returned.
    """
    cds_codons = fasta_gff.get_cds_codons_from_fasta(
        TEST_FASTA_FILE,
        TEST_GFF_FILE,
        True)
    cds_codons_minus_stops = {
        name: codons[:-1] for name, codons in TEST_CDS_CODONS.items()
    }
    assert cds_codons == cds_codons_minus_stops


def check_feature_codons_df(feature_codons, df):
    """
    Check contents of given dictionary with codons for features
    those of given DataFrame output by
    :py:func:`riboviz.fasta_gff.feature_codons_to_df`.

    :param feature_codons: Dictionary keyed by feature name and with \
    values that are lists of codons
    :type feature_codons: dict
    :param df: Pandas DataFrame
    :rtype: pandas.core.frame.DataFrame
    """
    num_rows, num_columns = df.shape
    assert num_columns == 3, "Unexpected number of columns"
    for column in [fasta_gff.GENE, fasta_gff.POS, fasta_gff.CODON]:
        assert column in list(df.columns),\
            "Missing column {}".format(column)
    total_codon_length = sum([len(codons)
                              for codons in feature_codons.values()])
    assert num_rows == total_codon_length, \
        "Unexpected number of rows"
    for feature_name, codons in feature_codons.items():
        feature_df = df.loc[df[fasta_gff.GENE] == feature_name]
        assert not feature_df.empty, "Missing row for {}".format(
            feature_name)
        num_rows, _ = feature_df.shape
        assert num_rows == len(codons), \
            "Unexpected number of rows for {}".format(feature_name)
        for pos, codon in zip(feature_df[fasta_gff.POS],
                              feature_df[fasta_gff.CODON]):
            # POS is 1-indexed.
            assert codons[pos - 1] == codon


def test_feature_codons_to_df_empty():
    """
    Test :py:func:`riboviz.fasta_gff.feature_codons_to_df`
    with no values produces an empty data frame.
    """
    df = fasta_gff.feature_codons_to_df({})
    check_feature_codons_df({}, df)


def test_feature_codons_to_df():
    """
    Test :py:func:`riboviz.fasta_gff.feature_codons_to_df`
    with values produces a data frame with the expected columns,
    rows and values.
    """
    feature_codons = {
        "G1": fasta_gff.sequence_to_codons("ATGAAATAA"),
        "G2": fasta_gff.sequence_to_codons("ATGGGGCCCTAG")
    }
    df = fasta_gff.feature_codons_to_df(feature_codons)
    check_feature_codons_df(feature_codons, df)


def test_get_cds_codons_file_empty_fasta(tmp_file):
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_file` with an empty
    FASTA file and GFF file (:py:const:`TEST_GFF_FILE`). A
    header-only TSV file is expected as output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    # Use tmp_file as both empty FASTA input file and TSV output
    # file
    fasta_gff.get_cds_codons_file(tmp_file, TEST_GFF_FILE, tmp_file)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    check_feature_codons_df({}, df)


def test_get_cds_codons_file(tmp_file):
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_file` with
    FASTA file (:py:const:`TEST_FASTA_FILE`) and GFF file
    (:py:const:`TEST_GFF_FILE`) and validate the TSV file output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    fasta_gff.get_cds_codons_file(TEST_FASTA_FILE,
                                  TEST_GFF_FILE,
                                  tmp_file)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    check_feature_codons_df(TEST_CDS_CODONS, df)


def test_get_cds_codons_file_exclude_stop_codon(tmp_file):
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_file` with
    FASTA file (:py:const:`TEST_FASTA_FILE`) and GFF file
    (:py:const:`TEST_GFF_FILE`), where stop codons are to be
    be excluded from the codons returned, and validate the TSV file
    output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    fasta_gff.get_cds_codons_file(TEST_FASTA_FILE,
                                  TEST_GFF_FILE,
                                  tmp_file,
                                  True)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    cds_codons_minus_stops = {
        name: codons[:-1] for name, codons in TEST_CDS_CODONS.items()
    }
    check_feature_codons_df(cds_codons_minus_stops, df)
