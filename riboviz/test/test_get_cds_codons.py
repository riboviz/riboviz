"""
:py:mod:`riboviz.get_cds_codons` tests.
"""
import os
import tempfile
import pytest
import pandas as pd
from riboviz import get_cds_codons
from riboviz.test import data

TEST_FASTA_CODONS_FILE = os.path.join(os.path.dirname(data.__file__),
                                      "test_get_cds_codons.fasta")
""" Test FASTA file in :py:mod:`riboviz.test.data`. """
TEST_GFF_CODONS_FILE = os.path.join(os.path.dirname(data.__file__),
                                    "test_get_cds_codons.gff")
""" Test GFF file in :py:mod:`riboviz.test.data`. """
TEST_GFF_NO_CDS_FILE = os.path.join(
    os.path.dirname(data.__file__),
    "test_get_cds_codons_no_codons.gff")
""" Test GFF file in :py:mod:`riboviz.test.data` with no CDS. """
TEST_CDS_CODONS = {
    "YAL001C_CDS": ["ATG", "GCC", "CAC", "TGT", "TAA"],
    "YAL002C_CDS": ["ATG", "GTA", "TCA", "GGA", "TAG"],
    "YAL004CSingleCodonCDS_CDS": ["ATG", "AGA", "TGA"],
    "YAL005CMultiCDS_CDS_1": ["ATG", "AGA", "TGA"],
    "YAL005CMultiCDS_CDS_2": ["ATG", "GAT", "TAC", "TAG"],
    "YAL005CMultiCDS_CDS_3": ["ATG", "CCA", "ATT", "TGA"],
    "YAL006CEmptyCDS_CDS": ["ATG", "TGA"],
    "YAL008CNoIdNameAttrCDS_mRNA_CDS": ["ATG", "GCC", "CAC", "TGT", "TAA"],
    "YAL009CIdNameAttrIdCDS_CDS": ["ATG", "AGA", "TGA"],
    "YAL010CMultiDuplicateCDS_CDS": ["ATG", "AGA", "TGA"],
    "YAL010CMultiDuplicateCDS_CDS.1": ["ATG", "GAT", "TAC", "TAG"],
    "YAL010CMultiDuplicateCDS_CDS.2": ["ATG", "CCA", "ATT", "TGA"]
}
"""
Expected codons for CDS in FASTA file
(:py:const:`TEST_FASTA_CODONS_FILE`) as defined in GFF file
(:py:const:`TEST_GFF_CODONS_FILE`).
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

    def sequence(self, _):
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
    Test :py:func:`riboviz.get_cds_codons.sequence_to_codons`
    with valid sequences.

    :param seq_codons: sequence and expected codons
    :type seq_codons: tuple(str or unicode, list(str or unicode))
    """
    sequence, codons = seq_codons
    assert get_cds_codons.sequence_to_codons(sequence) == codons,\
        "Unexpected codons"


@pytest.mark.parametrize("feature", [
    ("SeqID_mRNA", {}, False, None),
    ("SeqID_mRNA", {"Name": []}, False, None),
    ("SeqID_mRNA", {"Name": ["SeqName_CDS"]}, False, "SeqName_CDS"),
    ("SeqID_mRNA", {"ID": []}, False, None),
    ("SeqID_mRNA", {"ID": ["SeqID_CDS"]}, False, "SeqID_CDS"),
    ("SeqID_mRNA", {"ID": ["SeqID_CDS"], "Name": ["SeqName_CDS"]},
     False, "SeqID_CDS"),
    ("SeqID_mRNA", {}, True, None),
    ("SeqID_mRNA", {"Name": []}, True, None),
    ("SeqID_mRNA", {"Name": ["SeqName_CDS"]}, True, "SeqName_CDS"),
    ("SeqID_mRNA", {"ID": []}, True, None),
    ("SeqID_mRNA", {"ID": ["SeqID_CDS"]}, True, "SeqID_CDS"),
    ("SeqID_mRNA", {"ID": ["SeqID_CDS"], "Name": ["SeqName_CDS"]},
     True, "SeqName_CDS")])
def test_get_feature_id(feature):
    """
    Test :py:func:`riboviz.get_cds_codons.get_feature_id`.

    :param feature: sequence ID, attributes, use 'ID' as name if both \
    'ID' and 'Name' are defined, expected feature ID
    :type feature: tuple(str or unicode, str or unicode, bool, \
    str or unicode)
    """
    seq_id, attributes, report_name, feature_id = feature
    assert get_cds_codons.get_feature_id(
        MockFeature(seq_id, "", attributes), report_name) == feature_id


def test_get_cds_from_fasta():
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_from_fasta` returns
    a sequence.

    :py:class:`MockFeature` is used to mock
    ``gffutils.feature.Feature`` to avoid the need to use a FASTA
    file for this test.
    """
    expected_sequence = "ATGAAATAA"
    feature = MockFeature("SeqID", expected_sequence)
    # Any file name can be used as MockFeature ignores it and returns
    # the sequence given above.
    sequence = get_cds_codons.get_cds_from_fasta(feature, "test.fasta")
    assert sequence == expected_sequence


def test_get_cds_from_fasta_bad_length():
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_from_fasta` raises an
    error if given a sequence whose length is not divisible by 3.

    :py:class:`MockFeature` is used to mock
    ``gffutils.feature.Feature`` to avoid the need to use a FASTA
    file for this test.
    """
    feature = MockFeature("SeqID", "ATGATAA")
    with pytest.raises(AssertionError):
        # Any file name can be used as MockFeature ignores it
        # and returns the sequence given above.
        get_cds_codons.get_cds_from_fasta(feature, "test.fasta")


def test_get_cds_codons_from_fasta_empty(tmp_file):
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_from_fasta`
    with an empty FASTA file and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`).

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    # Use tmp_file as empty FASTA input file.
    cds_codons = get_cds_codons.get_cds_codons_from_fasta(
        tmp_file,
        TEST_GFF_CODONS_FILE)
    assert cds_codons == {}


def test_get_cds_codons_from_fasta_no_cds():
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_NO_CDS_FILE`) which defines no CDS.
    """
    cds_codons = get_cds_codons.get_cds_codons_from_fasta(
        TEST_FASTA_CODONS_FILE,
        TEST_GFF_NO_CDS_FILE)
    assert cds_codons == {}


def test_get_cds_codons_from_fasta():
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`).
    """
    cds_codons = get_cds_codons.get_cds_codons_from_fasta(
        TEST_FASTA_CODONS_FILE,
        TEST_GFF_CODONS_FILE)
    assert cds_codons == TEST_CDS_CODONS


def test_get_cds_codons_from_fasta_report_name_true():
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`) and ``report_name=True``.
    """
    cds_codons = get_cds_codons.get_cds_codons_from_fasta(
        TEST_FASTA_CODONS_FILE,
        TEST_GFF_CODONS_FILE,
        report_name=True)
    # Update TEST_CDS_CODONS with the expected result when
    # report_name_True.
    test_cds_codons = TEST_CDS_CODONS.copy()
    codons = test_cds_codons["YAL009CIdNameAttrIdCDS_CDS"]
    del test_cds_codons["YAL009CIdNameAttrIdCDS_CDS"]
    test_cds_codons["YAL009CIdNameAttrNameCDS_CDS"] = codons
    assert cds_codons == test_cds_codons


def test_get_cds_codons_from_fasta_exclude_stop_codons():
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_from_fasta`
    with FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`) where stop codons are to be
    excluded from the codons returned.
    """
    cds_codons = get_cds_codons.get_cds_codons_from_fasta(
        TEST_FASTA_CODONS_FILE,
        TEST_GFF_CODONS_FILE,
        exclude_stop_codons=True)
    cds_codons_minus_stops = {
        name: codons[:-1] for name, codons in TEST_CDS_CODONS.items()
    }
    assert cds_codons == cds_codons_minus_stops


def check_feature_codons_df(feature_codons, df):
    """
    Check contents of given dictionary with codons for features
    those of given DataFrame output by
    :py:func:`riboviz.get_cds_codons.feature_codons_to_df`.

    :param feature_codons: Dictionary keyed by feature name and with \
    values that are lists of codons
    :type feature_codons: dict
    :param df: Pandas DataFrame
    :rtype: pandas.core.frame.DataFrame
    """
    num_rows, num_columns = df.shape
    assert num_columns == 3, "Unexpected number of columns"
    for column in [get_cds_codons.GENE, get_cds_codons.POS, get_cds_codons.CODON]:
        assert column in list(df.columns),\
            "Missing column {}".format(column)
    total_codon_length = sum([len(codons)
                              for codons in feature_codons.values()])
    assert num_rows == total_codon_length, \
        "Unexpected number of rows"
    for feature_name, codons in feature_codons.items():
        feature_df = df.loc[df[get_cds_codons.GENE] == feature_name]
        assert not feature_df.empty, "Missing row for {}".format(
            feature_name)
        num_rows, _ = feature_df.shape
        assert num_rows == len(codons), \
            "Unexpected number of rows for {}".format(feature_name)
        for pos, codon in zip(feature_df[get_cds_codons.POS],
                              feature_df[get_cds_codons.CODON]):
            # POS is 1-indexed.
            assert codons[pos - 1] == codon


def test_feature_codons_to_df_empty():
    """
    Test :py:func:`riboviz.get_cds_codons.feature_codons_to_df`
    with no values produces an empty data frame.
    """
    df = get_cds_codons.feature_codons_to_df({})
    check_feature_codons_df({}, df)


def test_feature_codons_to_df():
    """
    Test :py:func:`riboviz.get_cds_codons.feature_codons_to_df`
    with values produces a data frame with the expected columns,
    rows and values.
    """
    feature_codons = {
        "G1": get_cds_codons.sequence_to_codons("ATGAAATAA"),
        "G2": get_cds_codons.sequence_to_codons("ATGGGGCCCTAG")
    }
    df = get_cds_codons.feature_codons_to_df(feature_codons)
    check_feature_codons_df(feature_codons, df)


def test_get_cds_codons_file_empty_fasta(tmp_file):
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_file` with an empty
    FASTA file and GFF file (:py:const:`TEST_GFF_CODONS_FILE`). A
    header-only TSV file is expected as output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    # Use tmp_file as both empty FASTA input file and TSV output
    # file
    get_cds_codons.get_cds_codons_file(tmp_file, TEST_GFF_CODONS_FILE, tmp_file)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    check_feature_codons_df({}, df)


def test_get_cds_codons_file(tmp_file):
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_file` with
    FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`) and validate the TSV file output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    get_cds_codons.get_cds_codons_file(TEST_FASTA_CODONS_FILE,
                                       TEST_GFF_CODONS_FILE,
                                       tmp_file)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    check_feature_codons_df(TEST_CDS_CODONS, df)


def test_get_cds_codons_file_cds_format(tmp_file):
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_file` with
    FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`) and a custom
    CDS feature name format and validate the TSV file output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    get_cds_codons.get_cds_codons_file(TEST_FASTA_CODONS_FILE,
                                       TEST_GFF_CODONS_FILE,
                                       tmp_file,
                                       cds_feature_format="{}-Custom")
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    # Delete the test entry that uses the default CDS_FEATURE_FORMAT.
    test_cds_codons = TEST_CDS_CODONS.copy()
    codons = test_cds_codons["YAL008CNoIdNameAttrCDS_mRNA_CDS"]
    del test_cds_codons["YAL008CNoIdNameAttrCDS_mRNA_CDS"]
    # Add the expected result for the custom entry.
    test_cds_codons["YAL008CNoIdNameAttrCDS_mRNA-Custom"] = codons
    check_feature_codons_df(test_cds_codons, df)


def test_get_cds_codons_file_exclude_stop_codons(tmp_file):
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_file` with
    FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`), where stop codons are to be
    be excluded from the codons returned, and validate the TSV file
    output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    get_cds_codons.get_cds_codons_file(TEST_FASTA_CODONS_FILE,
                                       TEST_GFF_CODONS_FILE,
                                       tmp_file,
                                       exclude_stop_codons=True)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    cds_codons_minus_stops = {
        name: codons[:-1] for name, codons in TEST_CDS_CODONS.items()
    }
    check_feature_codons_df(cds_codons_minus_stops, df)


def test_get_cds_codons_file_report_name_true(tmp_file):
    """
    Test :py:func:`riboviz.get_cds_codons.get_cds_codons_file` with
    FASTA file (:py:const:`TEST_FASTA_CODONS_FILE`) and GFF file
    (:py:const:`TEST_GFF_CODONS_FILE`) and ``report_name=True``.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    get_cds_codons.get_cds_codons_file(TEST_FASTA_CODONS_FILE,
                                       TEST_GFF_CODONS_FILE,
                                       tmp_file,
                                       report_name=True)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    # Update TEST_CDS_CODONS with the expected result when
    # report_name=True.
    test_cds_codons = TEST_CDS_CODONS.copy()
    codons = test_cds_codons["YAL009CIdNameAttrIdCDS_CDS"]
    del test_cds_codons["YAL009CIdNameAttrIdCDS_CDS"]
    test_cds_codons["YAL009CIdNameAttrNameCDS_CDS"] = codons
    check_feature_codons_df(test_cds_codons, df)
