"""
:py:mod:`riboviz.test_fasta_gff` tests.
"""
import os
import tempfile
import pytest
import pandas as pd
from riboviz import fasta_gff
from riboviz.test import data


TEST_FASTA_CHECK_FILE = os.path.join(os.path.dirname(data.__file__),
                                     "test_fasta_gff_check.fasta")
""" Test FASTA file in :py:mod:`riboviz.test.data`. """
TEST_GFF_CHECK_FILE = os.path.join(os.path.dirname(data.__file__),
                                   "test_fasta_gff_check.gff")
""" Test GFF file in :py:mod:`riboviz.test.data`. """
TEST_CHECK_ISSUES = [
    ("YAL003CMissingGene_mRNA", "",
     fasta_gff.ISSUE_MISSING_SEQUENCE),
    ("YAL004CNoIDNameAttr_mRNA", "Undefined",
     fasta_gff.ISSUE_NO_ID_NAME),
    ("YAL006CNonUniqueID_mRNA", "YAL005_7CNonUniqueID_CDS",
     fasta_gff.ISSUE_DUPLICATE_FEATURE_ID),
    ("YAL007CNonUniqueID_mRNA", "YAL005_7CNonUniqueID_CDS",
     fasta_gff.ISSUE_DUPLICATE_FEATURE_ID),
    ("YAL008CBadLengthNoStop_mRNA", "YAL008CBadLengthNoStop_CDS",
     fasta_gff.ISSUE_INCOMPLETE),
    ("YAL008CBadLengthNoStop_mRNA", "YAL008CBadLengthNoStop_CDS",
     fasta_gff.ISSUE_NO_STOP),
    ("YAL009C_NoATGStart_mRNA", "YAL09CNoATGStart_CDS",
     fasta_gff.ISSUE_NO_ATG_START),
    ("YAL010C_NoStop_mRNA", "YAL010CNoStop_CDS",
     fasta_gff.ISSUE_NO_STOP),
    ("YAL011C_InternalStop_mRNA", "YAL011CInternalStop_CDS",
     fasta_gff.ISSUE_INTERNAL_STOP),
    ("YAL012C_NoATGStartNoStop_mRNA", "YAL12CNoATGStartNoStop_CDS",
     fasta_gff.ISSUE_NO_ATG_START),
    ("YAL012C_NoATGStartNoStop_mRNA", "YAL12CNoATGStartNoStop_CDS",
     fasta_gff.ISSUE_NO_STOP),
    ("YAL013C_NoATGStartInternalStop_mRNA",
     "YAL13CNoATGStartInternalStop_CDS",
     fasta_gff.ISSUE_NO_ATG_START),
    ("YAL013C_NoATGStartInternalStop_mRNA",
     "YAL13CNoATGStartInternalStop_CDS",
     fasta_gff.ISSUE_INTERNAL_STOP),
    ("YAL014C_NoATGStartInternalStopNoStop_mRNA",
     "YAL14CNoATGStartInternalStopNoStop_CDS",
     fasta_gff.ISSUE_NO_ATG_START),
    ("YAL014C_NoATGStartInternalStopNoStop_mRNA",
     "YAL14CNoATGStartInternalStopNoStop_CDS",
     fasta_gff.ISSUE_NO_STOP),
    ("YAL014C_NoATGStartInternalStopNoStop_mRNA",
     "YAL14CNoATGStartInternalStopNoStop_CDS",
     fasta_gff.ISSUE_INTERNAL_STOP),
    ("YAL015CMultiCDS_mRNA", "", fasta_gff.ISSUE_MULTIPLE_CDS)
]
"""
Expected issues (sequence ID, feature ID, issue) for checking FASTA
file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
(:py:const:`TEST_GFF_CHECK_FILE`).
"""

TEST_FASTA_CODONS_FILE = os.path.join(os.path.dirname(data.__file__),
                                      "test_fasta_gff_codons.fasta")
""" Test FASTA file in :py:mod:`riboviz.test.data`. """
TEST_GFF_CODONS_FILE = os.path.join(os.path.dirname(data.__file__),
                                    "test_fasta_gff_codons.gff")
""" Test GFF file in :py:mod:`riboviz.test.data`. """
TEST_GFF_NO_CDS_FILE = os.path.join(
    os.path.dirname(data.__file__),
    "test_fasta_gff_no_codons.gff")
""" Test GFF file in :py:mod:`riboviz.test.data` with no CDS. """
TEST_CDS_CODONS = {
    "YAL001C_CDS": ["ATG", "GCC", "CAC", "TGT", "TAA"],
    "YAL002C_CDS": ["ATG", "GTA", "TCA", "GGA", "TAG"],
    "YAL004CSingleCodonCDS_CDS": ["ATG", "AGA", "TGA"],
    "YAL005CMultiCDS_CDS_1": ["ATG", "AGA", "TGA"],
    "YAL005CMultiCDS_CDS_2": ["ATG", "GAT", "TAC", "TAG"],
    "YAL005CMultiCDS_CDS_3": ["ATG", "CCA", "ATT", "TGA"],
    "YAL006CEmptyCDS_CDS": ["ATG", "TGA"],
    "YAL008CNoIdNameAttr_mRNA_CDS": ["ATG", "GCC", "CAC", "TGT", "TAA"],
    "YAL009CMultiDuplicateCDS_CDS": ["ATG", "AGA", "TGA"],
    "YAL009CMultiDuplicateCDS_CDS.1": ["ATG", "GAT", "TAC", "TAG"],
    "YAL009CMultiDuplicateCDS_CDS.2": ["ATG", "CCA", "ATT", "TGA"]
}
"""
Expected codons for CDS in FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`)
as defined in GFF file (:py:const:`TEST_GFF_CODONS_FILE`).
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


def test_get_fasta_gff_cds_issues():
    """
    Test :py:func:`riboviz.fasta_gff.get_fasta_gff_cds_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and check all issues match
    expected issues in :py:const:`TEST_CHECK_ISSUES`).
    """
    issues = fasta_gff.get_fasta_gff_cds_issues(
        TEST_FASTA_CHECK_FILE,
        TEST_GFF_CHECK_FILE)
    for issue in issues:
        assert issue in TEST_CHECK_ISSUES


def check_fasta_gff_issues_df(issues, df):
    """
    Check contents of given list of tuples with issues held within
    a DataFrame output by
    :py:func:`riboviz.fasta_gff.fasta_gff_issues_to_df`.

    :param issues: List of issues for sequences and features.
    :type issues: list(tuple(str or unicode, str or unicode, \
    str or unicode))
    :param df: Pandas DataFrame
    :rtype: pandas.core.frame.DataFrame
    """
    num_rows, num_columns = df.shape
    assert num_columns == 3, "Unexpected number of columns"
    for column in [fasta_gff.SEQUENCE, fasta_gff.FEATURE, fasta_gff.ISSUE]:
        assert column in list(df.columns),\
            "Missing column {}".format(column)
    assert num_rows == len(issues), \
        "Unexpected number of rows"
    for sequence, feature, issue in issues:
        issue_df = df[(df[fasta_gff.SEQUENCE] == sequence) &
                      (df[fasta_gff.FEATURE] == feature) &
                      (df[fasta_gff.ISSUE] == issue)]
        assert not issue_df.empty, "Missing row for {}".format(issue)


def test_fasta_gff_issues_to_df():
    """
    Test :py:func:`riboviz.fasta_gff.get_fasta_gff_cds_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and check all issues match
    expected issues in :py:const:`TEST_CHECK_ISSUES`).
    """
    df = fasta_gff.fasta_gff_issues_to_df(TEST_CHECK_ISSUES)
    check_fasta_gff_issues_df(TEST_CHECK_ISSUES, df)

# TODO
# def check_fasta_gff(fasta, gff, feature_issues_file, delimiter="\t"):
# TODO


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
def test_get_feature_id(features):
    """
    Test :py:func:`riboviz.fasta_gff.get_feature_id`.

    :param features: sequence ID, attributes, expected feature ID
    :type features: tuple(str or unicode, str or unicode, str or unicode)
    """
    seq_id, attributes, feature_id = features
    assert fasta_gff.get_feature_id(
        MockFeature(seq_id, "", attributes)) == feature_id


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
    (:py:const:`TEST_GFF_CODONS_FILE`).

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    # Use tmp_file as both empty FASTA input file.
    cds_codons = fasta_gff.get_cds_codons_from_fasta(
        tmp_file,
        TEST_GFF_CODONS_FILE)
    assert cds_codons == {}


def test_get_cds_codons_from_fasta_no_cds(tmp_file):
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_NO_CDS_FILE`) which defines no CDS.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    cds_codons = fasta_gff.get_cds_codons_from_fasta(
        TEST_FASTA_CODONS_FILE,
        TEST_GFF_NO_CDS_FILE)
    assert cds_codons == {}


def test_get_cds_codons_from_fasta():
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`).
    """
    cds_codons = fasta_gff.get_cds_codons_from_fasta(
        TEST_FASTA_CODONS_FILE,
        TEST_GFF_CODONS_FILE)
    assert cds_codons == TEST_CDS_CODONS


def test_get_cds_codons_from_fasta_exclude_stop_codon():
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`) where stop codons are to be
    excluded from the codons returned.
    """
    cds_codons = fasta_gff.get_cds_codons_from_fasta(
        TEST_FASTA_CODONS_FILE,
        TEST_GFF_CODONS_FILE,
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
    FASTA file and GFF file (:py:const:`TEST_GFF_CODONS_FILE`). A
    header-only TSV file is expected as output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    # Use tmp_file as both empty FASTA input file and TSV output
    # file
    fasta_gff.get_cds_codons_file(tmp_file, TEST_GFF_CODONS_FILE, tmp_file)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    check_feature_codons_df({}, df)


def test_get_cds_codons_file(tmp_file):
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_file` with
    FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`) and validate the TSV file output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    fasta_gff.get_cds_codons_file(TEST_FASTA_CODONS_FILE,
                                  TEST_GFF_CODONS_FILE,
                                  tmp_file)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    check_feature_codons_df(TEST_CDS_CODONS, df)


def test_get_cds_codons_file_exclude_stop_codon(tmp_file):
    """
    Test :py:func:`riboviz.fasta_gff.get_cds_codons_file` with
    FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`), where stop codons are to be
    be excluded from the codons returned, and validate the TSV file
    output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    fasta_gff.get_cds_codons_file(TEST_FASTA_CODONS_FILE,
                                  TEST_GFF_CODONS_FILE,
                                  tmp_file,
                                  True)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    cds_codons_minus_stops = {
        name: codons[:-1] for name, codons in TEST_CDS_CODONS.items()
    }
    check_feature_codons_df(cds_codons_minus_stops, df)
