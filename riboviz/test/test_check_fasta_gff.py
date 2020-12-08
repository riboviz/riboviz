"""
:py:mod:`riboviz.check_fasta_gff` tests.
"""
import os
import tempfile
import pytest
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from riboviz import check_fasta_gff
from riboviz.fasta_gff import CDS_FEATURE_FORMAT
from riboviz.test import data


TEST_FASTA_CHECK_FILE = os.path.join(os.path.dirname(data.__file__),
                                     "test_check_fasta_gff.fasta")
""" Test FASTA file in :py:mod:`riboviz.test.data`. """
TEST_GFF_CHECK_FILE = os.path.join(os.path.dirname(data.__file__),
                                   "test_check_fasta_gff.gff")
""" Test GFF file in :py:mod:`riboviz.test.data`. """
TEST_NO_ID_NAME_ATTR_PREFIX = "YAL004CNoIDNameAttr_mRNA"
""" Prefix for name of feature in test issue. """
TEST_NO_ID_NAME_ATTR_ISSUE = (
    "YAL004CNoIDNameAttr_mRNA",
    CDS_FEATURE_FORMAT.format(TEST_NO_ID_NAME_ATTR_PREFIX),
    check_fasta_gff.NO_ID_NAME, None)
"""
Test issue with feature name derived from
:py:const:`TEST_NO_ID_NAME_ATTR_PREFIX` and
:py:const:`riboviz.fasta_gff.CDS_FEATURE_FORMAT`.
"""
TEST_NO_ATG_START_ID_NAME_ATTR_ID = "YAL010CNoATGStartID_CDS"
""" Name of feature in test issue, derived from feature ID attribute. """
TEST_NO_ATG_START_ID_NAME_ATTR_NAME = "YAL010CNoATGStartName_CDS"
""" Name of feature in test issue, derived from feature Name attribute. """
TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE = (
    "YAL010CNoATGStartIDNameAttr_mRNA",
    TEST_NO_ATG_START_ID_NAME_ATTR_ID,
    check_fasta_gff.NO_START_CODON, "AAG")
""" Test issue. """
TEST_CHECK_GFF_ISSUES = [
    TEST_NO_ID_NAME_ATTR_ISSUE,
    ("YAL005CNonUniqueID_mRNA", "YAL005_7CNonUniqueID_CDS",
     check_fasta_gff.DUPLICATE_FEATURE_ID, None),
    ("YAL006CNonUniqueID_mRNA", "YAL005_7CNonUniqueID_CDS",
     check_fasta_gff.DUPLICATE_FEATURE_ID, None),
    ("YAL007CNonUniqueID_mRNA", "YAL005_7CNonUniqueID_CDS",
     check_fasta_gff.DUPLICATE_FEATURE_ID, None),
    (check_fasta_gff.WILDCARD, "YAL005_7CNonUniqueID_CDS",
     check_fasta_gff.DUPLICATE_FEATURE_IDS, 3),
    ("YAL016CMultiCDS_mRNA", check_fasta_gff.WILDCARD,
     check_fasta_gff.MULTIPLE_CDS, 3),
]
"""
Expected GFF-specific issues (sequence ID, feature ID, issue type,
issue data) for GFF file (:py:const:`TEST_GFF_CHECK_FILE`) only.
"""
TEST_CHECK_FASTA_ISSUES = [
    ("YAL003CMissingSequence_mRNA", check_fasta_gff.NOT_APPLICABLE,
     check_fasta_gff.SEQUENCE_NOT_IN_FASTA, None),
    ("YAL008CBadLengthNoStop_mRNA", "YAL008CBadLengthNoStop_CDS",
     check_fasta_gff.INCOMPLETE_FEATURE, None),
    ("YAL008CBadLengthNoStop_mRNA", "YAL008CBadLengthNoStop_CDS",
     check_fasta_gff.NO_STOP_CODON, "TAN"),
    ("YAL009CNoATGStart_mRNA", "YAL009CNoATGStart_CDS",
     check_fasta_gff.NO_START_CODON, "AAG"),
    TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE,
    ("YAL011CNoStop_mRNA", "YAL011CNoStop_CDS",
     check_fasta_gff.NO_STOP_CODON, "AAG"),
    ("YAL012CInternalStop_mRNA", "YAL012CInternalStop_CDS",
     check_fasta_gff.INTERNAL_STOP_CODON, None),
    ("YAL013CNoATGStartNoStop_mRNA", "YAL013CNoATGStartNoStop_CDS",
     check_fasta_gff.NO_START_CODON, "AAG"),
    ("YAL013CNoATGStartNoStop_mRNA", "YAL013CNoATGStartNoStop_CDS",
     check_fasta_gff.NO_STOP_CODON, "AAG"),
    ("YAL014CNoATGStartInternalStop_mRNA",
     "YAL014CNoATGStartInternalStop_CDS",
     check_fasta_gff.NO_START_CODON, "AAG"),
    ("YAL014CNoATGStartInternalStop_mRNA",
     "YAL014CNoATGStartInternalStop_CDS",
     check_fasta_gff.INTERNAL_STOP_CODON, None),
    ("YAL015CNoATGStartInternalStopNoStop_mRNA",
     "YAL015CNoATGStartInternalStopNoStop_CDS",
     check_fasta_gff.NO_START_CODON, "AAG"),
    ("YAL015CNoATGStartInternalStopNoStop_mRNA",
     "YAL015CNoATGStartInternalStopNoStop_CDS",
     check_fasta_gff.NO_STOP_CODON, "AAG"),
    ("YAL015CNoATGStartInternalStopNoStop_mRNA",
     "YAL015CNoATGStartInternalStopNoStop_CDS",
     check_fasta_gff.INTERNAL_STOP_CODON, None),
    ("YAL017CMissingGFF_mRNA", check_fasta_gff.NOT_APPLICABLE,
     check_fasta_gff.SEQUENCE_NOT_IN_GFF, None),
    ("YAL018CMissingGFF_mRNA", check_fasta_gff.NOT_APPLICABLE,
     check_fasta_gff.SEQUENCE_NOT_IN_GFF, None),
    ("YAL019CMissingSequence_mRNA", check_fasta_gff.NOT_APPLICABLE,
     check_fasta_gff.SEQUENCE_NOT_IN_FASTA, None),
]
"""
Expected FASTA-specific issues (sequence ID, feature ID, issue type,
issue data) for checking FASTA file
(:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
(:py:const:`TEST_GFF_CHECK_FILE`).
"""
TEST_CHECK_ISSUES = TEST_CHECK_GFF_ISSUES + TEST_CHECK_FASTA_ISSUES
"""
All expected issues (sequence ID, feature ID, issue) for checking
FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
(:py:const:`TEST_GFF_CHECK_FILE`).
"""


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


def test_get_fasta_sequence_ids(tmp_file):
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_sequence_ids`
    with a FASTA file and check that the expected sequence IDs
    are returned.
    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    fasta_seq_ids = ["SeqID_{}".format(i) for i in range(0, 5)]
    with open(tmp_file, "w") as f:
        for seq_id in fasta_seq_ids:
            record = SeqRecord(Seq("GATTACCA"),
                               id=seq_id,
                               description="test")
            SeqIO.write(record, f, "fasta")
    seq_ids = check_fasta_gff.get_fasta_sequence_ids(tmp_file)
    assert set(fasta_seq_ids) == seq_ids, "Unexpected set of sequence IDs"


def test_get_fasta_sequence_ids_empty_fasta(tmp_file):
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_sequence_ids`
    with an empty FASTA file and check that an empty set is
    returned.
    """
    seq_ids = check_fasta_gff.get_fasta_sequence_ids(tmp_file)
    assert not seq_ids, "Expected empty set of sequence IDs"


def test_get_fasta_gff_cds_issues():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_gff_cds_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and check all issues match
    expected issues in :py:const:`TEST_CHECK_ISSUES`).
    """
    issues = check_fasta_gff.get_fasta_gff_cds_issues(
        TEST_FASTA_CHECK_FILE,
        TEST_GFF_CHECK_FILE)
    for issue in issues:
        assert issue in TEST_CHECK_ISSUES


def test_get_fasta_gff_cds_issues_feature_format():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_gff_cds_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and a custom feature name.
    """
    issues = check_fasta_gff.get_fasta_gff_cds_issues(
        TEST_FASTA_CHECK_FILE,
        TEST_GFF_CHECK_FILE,
        feature_format="{}-Custom")
    # Update TEST_CHECK_ISSUES with the expected result when
    # custom cds_feature_format is used.
    test_check_issues = TEST_CHECK_ISSUES.copy()
    test_check_issues.remove(TEST_NO_ID_NAME_ATTR_ISSUE)
    test_check_issues.append(
        (TEST_NO_ID_NAME_ATTR_ISSUE[0],
         "{}-Custom".format(TEST_NO_ID_NAME_ATTR_PREFIX),
         TEST_NO_ID_NAME_ATTR_ISSUE[2],
         TEST_NO_ID_NAME_ATTR_ISSUE[3]))
    for issue in issues:
        assert issue in test_check_issues


def test_get_fasta_gff_cds_issues_use_feature_name_true():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_gff_cds_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and ``use_feature_name=True``.
    """
    issues = check_fasta_gff.get_fasta_gff_cds_issues(
        TEST_FASTA_CHECK_FILE,
        TEST_GFF_CHECK_FILE,
        use_feature_name=True)
    # Update TEST_CHECK_ISSUES with the expected result when
    # use_feature_name=True.
    test_check_issues = TEST_CHECK_ISSUES.copy()
    test_check_issues.remove(TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE)
    test_check_issues.append(
        (TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[0],
         TEST_NO_ATG_START_ID_NAME_ATTR_NAME,
         TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[2],
         TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[3]))
    for issue in issues:
        assert issue in test_check_issues


def test_get_fasta_gff_cds_issues_empty_fasta(tmp_file):
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_gff_cds_issues`
    with an empty FASTA file and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`)  and check all issues match
    expected issues in :py:const:`TEST_CHECK_GFF_ISSUES`).
    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    issues = check_fasta_gff.get_fasta_gff_cds_issues(
        tmp_file,
        TEST_GFF_CHECK_FILE)
    for issue in issues:
        assert issue in TEST_CHECK_GFF_ISSUES


def check_fasta_gff_issues_df(issues, df):
    """
    Check contents of given list of tuples with issues held within
    a DataFrame output by
    :py:func:`riboviz.check_fasta_gff.fasta_gff_issues_to_df`.

    :param issues: List of issues for sequences and features.
    :type issues: list(tuple(str or unicode, str or unicode, \
    str or unicode))
    :param df: Pandas DataFrame
    :rtype: pandas.core.frame.DataFrame
    """
    num_rows, num_columns = df.shape
    assert num_columns == 4, "Unexpected number of columns"
    for column in [check_fasta_gff.SEQUENCE,
                   check_fasta_gff.FEATURE,
                   check_fasta_gff.ISSUE_TYPE,
                   check_fasta_gff.ISSUE_DATA]:
        assert column in list(df.columns),\
            "Missing column {}".format(column)
    assert num_rows == len(issues), \
        "Unexpected number of rows"
    for sequence, feature, issue_type, issue_data in issues:
        issue_df = df[(df[check_fasta_gff.SEQUENCE] == sequence) &
                      (df[check_fasta_gff.FEATURE] == feature) &
                      (df[check_fasta_gff.ISSUE_TYPE] == issue_type)]
        assert not issue_df.empty, \
            "Expected 1 matching row for {} {} {}".format(
                sequence, feature, issue_type)
        if issue_data is None:
            expected_data = ""
        else:
            expected_data = issue_data
        df_issue_data = issue_df.iloc[0][check_fasta_gff.ISSUE_DATA]
        if type(expected_data) is int:
            df_issue_data = int(df_issue_data)
        elif type(expected_data) is float:
            df_issue_data = float(df_issue_data)
        else:
            df_issue_data = str(df_issue_data)
        assert df_issue_data == expected_data, \
            "Unexpected data ({}) for {} {} {}".format(
                df_issue_data, sequence, feature, issue_type)


def test_fasta_gff_issues_to_df():
    """
    Test :py:func:`riboviz.check_fasta_gff.fasta_gff_issues_to_df`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and check all issues match
    expected issues in :py:const:`TEST_CHECK_ISSUES`).
    """
    df = check_fasta_gff.fasta_gff_issues_to_df(TEST_CHECK_ISSUES)
    df = df.fillna("")  # Force None to "" not "nan"
    check_fasta_gff_issues_df(TEST_CHECK_ISSUES, df)


def test_fasta_gff_issues_to_df_empty():
    """
    Test :py:func:`riboviz.check_fasta_gff.fasta_gff_issues_to_df`
    with no values produces an empty data frame.
    """
    df = check_fasta_gff.fasta_gff_issues_to_df([])
    check_fasta_gff_issues_df([], df)


def test_check_fasta_gff_empty_fasta(tmp_file):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with an
    empty FASTA file and GFF file (:py:const:`TEST_GFF_CHECK_FILE`)
    and validate the TSV file output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    # Use tmp_file as both empty FASTA input file and TSV output
    # file
    check_fasta_gff.check_fasta_gff(tmp_file, TEST_GFF_CHECK_FILE, tmp_file)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    df = df.fillna("")  # Force None to "" not "nan"
    check_fasta_gff_issues_df(TEST_CHECK_GFF_ISSUES, df)


def test_check_fasta_gff(tmp_file):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with
    FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and validate the TSV file
    output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                    TEST_GFF_CHECK_FILE,
                                    tmp_file)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    df = df.fillna("")  # Force None to "" not "nan"
    check_fasta_gff_issues_df(TEST_CHECK_ISSUES, df)


def test_check_fasta_gff_feature_format(tmp_file):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with
    FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and a custom feature name
    format and validate the TSV file output.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                    TEST_GFF_CHECK_FILE,
                                    tmp_file,
                                    feature_format="{}-Custom")
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    df = df.fillna("")  # Force None to "" not "nan"
    # Update TEST_CHECK_ISSUES with the expected result when
    # custom cds_feature_format is used.
    test_check_issues = TEST_CHECK_ISSUES.copy()
    test_check_issues.remove(TEST_NO_ID_NAME_ATTR_ISSUE)
    test_check_issues.append(
        (TEST_NO_ID_NAME_ATTR_ISSUE[0],
         "{}-Custom".format(TEST_NO_ID_NAME_ATTR_PREFIX),
         TEST_NO_ID_NAME_ATTR_ISSUE[2],
         TEST_NO_ID_NAME_ATTR_ISSUE[3]))
    check_fasta_gff_issues_df(test_check_issues, df)


def test_check_fasta_gff_use_feature_name_true(tmp_file):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with
    FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and ``use_feature_name=True``.

    :param tmp_file: Temporary file
    :type tmp_file: str or unicode
    """
    check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                    TEST_GFF_CHECK_FILE,
                                    tmp_file,
                                    use_feature_name=True)
    df = pd.read_csv(tmp_file, delimiter="\t", comment="#")
    df = df.fillna("")  # Force None to "" not "nan"
    # Update TEST_CHECK_ISSUES with the expected result when
    # custom cds_feature_format is used.
    test_check_issues = TEST_CHECK_ISSUES.copy()
    test_check_issues.remove(TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE)
    test_check_issues.append(
        (TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[0],
         TEST_NO_ATG_START_ID_NAME_ATTR_NAME,
         TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[2],
         TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[3]))
    check_fasta_gff_issues_df(test_check_issues, df)
