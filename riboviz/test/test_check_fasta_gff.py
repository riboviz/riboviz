"""
:py:mod:`riboviz.check_fasta_gff` tests.
"""
import os
import pytest
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import FastaIndexingError
from riboviz import check_fasta_gff
from riboviz.fasta_gff import START_CODON
from riboviz.fasta_gff import CDS_FEATURE_FORMAT
from riboviz.test import data


TEST_FASTA_CHECK_FILE = os.path.join(os.path.dirname(data.__file__),
                                     "test_check_fasta_gff.fasta")
""" Test FASTA file in :py:mod:`riboviz.test.data`. """
TEST_GFF_CHECK_FILE = os.path.join(os.path.dirname(data.__file__),
                                   "test_check_fasta_gff.gff")
""" Test GFF file in :py:mod:`riboviz.test.data`. """
TEST_NUM_SEQUENCES = 17
""" Number of sequences in :py:const:`TEST_FASTA_CHECK_FILE`. """
TEST_NUM_FEATURES = 19
""" Number of features in :py:const:`TEST_GFF_CHECK_FILE`. """
TEST_NUM_CDS_FEATURES = 19
""" Number of CDS features in :py:const:`TEST_GFF_CHECK_FILE`. """
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
TEST_START_CODONS = [START_CODON, "AAG"]
""" List of start codons in test data. """
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


def test_get_fasta_sequence_ids(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_sequence_ids`
    with a FASTA file and check that the expected sequence IDs
    are returned.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    fasta_file = tmpdir.join("fasta.fasta")
    fasta_seq_ids = ["SeqID_{}".format(i) for i in range(0, 5)]
    with open(fasta_file, 'w') as f:
        for seq_id in fasta_seq_ids:
            record = SeqRecord(Seq("GATTACCA"),
                               id=seq_id,
                               description="test")
            SeqIO.write(record, f, "fasta")
    seq_ids = check_fasta_gff.get_fasta_sequence_ids(fasta_file)
    assert set(fasta_seq_ids) == seq_ids, "Unexpected set of sequence IDs"


def test_get_fasta_sequence_ids_no_such_fasta_file():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_sequence_ids`
    with a non-existent FASTA raises an exception.
    """
    with pytest.raises(FileNotFoundError):
        check_fasta_gff.get_fasta_sequence_ids("nosuch.fasta")


def test_get_fasta_sequence_ids_empty_fasta_file(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.get_fasta_sequence_ids`
    with an empty FASTA file and check that an empty set is
    returned.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    fasta_file = tmpdir.join("fasta.fasta")
    open(fasta_file, 'a').close()
    seq_ids = check_fasta_gff.get_fasta_sequence_ids(fasta_file)
    assert not seq_ids, "Expected empty set of sequence IDs"


def test_get_issues_no_such_fasta_file():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_issues`
    with an empty FASTA file and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) raises an exception.
    """
    with pytest.raises(FileNotFoundError):
        check_fasta_gff.get_issues("nosuch.fasta", TEST_GFF_CHECK_FILE)


def test_get_issues_no_such_gff_file():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and
    a non-existent GFF file raises an exception.
    """
    with pytest.raises(FileNotFoundError):
        check_fasta_gff.get_issues(TEST_FASTA_CHECK_FILE, "nosuch.gff")


def test_get_issues_empty_fasta_file(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.get_issues`
    with an empty FASTA file and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) raises an exception.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    # Cast to str to avoid "'LocalPath' object is not subscriptable"
    # error from gff library.
    fasta_file = str(tmpdir.join("fasta.fasta"))
    open(fasta_file, 'a').close()
    with pytest.raises(FastaIndexingError):
        check_fasta_gff.get_issues(fasta_file, TEST_GFF_CHECK_FILE)


def test_get_issues_empty_gff_file(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.get_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and
    an empty GFF file raises an exception.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    # Cast to str to avoid "'LocalPath' object is not subscriptable"
    # error from gff library.
    gff_file = str(tmpdir.join("gff.gff"))
    open(gff_file, 'a').close()
    with pytest.raises(ValueError):
        check_fasta_gff.get_issues(TEST_FASTA_CHECK_FILE, gff_file)


def test_get_issues():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and check all issues match
    expected issues in :py:const:`TEST_CHECK_ISSUES`).
    """
    num_sequences, num_features, num_cds_features, issues = \
        check_fasta_gff.get_issues(TEST_FASTA_CHECK_FILE,
                                   TEST_GFF_CHECK_FILE)
    assert num_sequences == TEST_NUM_SEQUENCES, \
        "Unexpected number of sequences"
    assert num_features == TEST_NUM_FEATURES, \
        "Unexpected number of features"
    assert num_cds_features == TEST_NUM_CDS_FEATURES, \
        "Unexpected number of CDS features"
    for issue in issues:
        assert issue in TEST_CHECK_ISSUES


def test_get_issues_use_feature_name_true():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and ``use_feature_name=True``.
    """
    num_sequences, num_features, num_cds_features, issues = \
        check_fasta_gff.get_issues(TEST_FASTA_CHECK_FILE,
                                   TEST_GFF_CHECK_FILE,
                                   use_feature_name=True)
    assert num_sequences == TEST_NUM_SEQUENCES, \
        "Unexpected number of sequences"
    assert num_features == TEST_NUM_FEATURES, \
        "Unexpected number of features"
    assert num_cds_features == TEST_NUM_CDS_FEATURES, \
        "Unexpected number of CDS features"
    # Create expected results when use_feature_name=True.
    test_check_issues = TEST_CHECK_ISSUES.copy()
    test_check_issues.remove(TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE)
    test_check_issues.append(
        (TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[0],
         TEST_NO_ATG_START_ID_NAME_ATTR_NAME,
         TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[2],
         TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[3]))
    for issue in issues:
        assert issue in test_check_issues


def test_get_issues_feature_format():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and custom
    ``cds_feature_format``.
    """
    cds_feature_format = "{}-Custom"
    num_sequences, num_features, num_cds_features, issues = \
        check_fasta_gff.get_issues(TEST_FASTA_CHECK_FILE,
                                   TEST_GFF_CHECK_FILE,
                                   feature_format=cds_feature_format)
    assert num_sequences == TEST_NUM_SEQUENCES, \
        "Unexpected number of sequences"
    assert num_features == TEST_NUM_FEATURES, \
        "Unexpected number of features"
    assert num_cds_features == TEST_NUM_CDS_FEATURES, \
        "Unexpected number of CDS features"
    # Create expected results for custom cds_feature_format.
    test_check_issues = TEST_CHECK_ISSUES.copy()
    test_check_issues.remove(TEST_NO_ID_NAME_ATTR_ISSUE)
    test_check_issues.append(
        (TEST_NO_ID_NAME_ATTR_ISSUE[0],
         cds_feature_format.format(TEST_NO_ID_NAME_ATTR_PREFIX),
         TEST_NO_ID_NAME_ATTR_ISSUE[2],
         TEST_NO_ID_NAME_ATTR_ISSUE[3]))
    for issue in issues:
        assert issue in test_check_issues


def test_get_issues_start_codons():
    """
    Test :py:func:`riboviz.check_fasta_gff.get_issues`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and custom ``start_codons``.
    """
    num_sequences, num_features, num_cds_features, issues = \
        check_fasta_gff.get_issues(TEST_FASTA_CHECK_FILE,
                                   TEST_GFF_CHECK_FILE,
                                   start_codons=TEST_START_CODONS)
    assert num_sequences == TEST_NUM_SEQUENCES, \
        "Unexpected number of sequences"
    assert num_features == TEST_NUM_FEATURES, \
        "Unexpected number of features"
    assert num_cds_features == TEST_NUM_CDS_FEATURES, \
        "Unexpected number of CDS features"
    # TEST_START_CODONS includes all start codons in test data so
    # expect no NO_START_CODON issues in results.
    test_check_issues = [(s, f, i, d)
                         for (s, f, i, d) in TEST_CHECK_ISSUES
                         if i != check_fasta_gff.NO_START_CODON]
    for issue in issues:
        assert issue in test_check_issues


def check_fasta_gff_issues_csv(issues, csv_file):
    """
    Check contents of given list of tuples with issues held within
    a CSV file written by
    :py:func:`riboviz.check_fasta_gff.write_issues_to_csv`.

    :param issues: List of issues for sequences and features.
    :type issues: list(tuple(str or unicode, str or unicode, \
    str or unicode))
    :param csv_file: CSV file
    :type csv_file: str or unicode
    """
    df = pd.read_csv(csv_file, delimiter="\t", comment="#")
    df = df.fillna("")  # Force None to "" not "nan"
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
        if isinstance(expected_data, int):
            df_issue_data = int(df_issue_data)
        elif isinstance(expected_data, float):
            df_issue_data = float(df_issue_data)
        else:
            df_issue_data = str(df_issue_data)
        assert df_issue_data == expected_data, \
            "Unexpected data ({}) for {} {} {}".format(
                df_issue_data, sequence, feature, issue_type)


def test_write_issues_to_csv(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.write_issues_to_csv`
    produces a CSV file with the expected columns, rows and values.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    issues_file = tmpdir.join("issues.tsv")
    check_fasta_gff.write_issues_to_csv(
        TEST_CHECK_ISSUES, issues_file)
    check_fasta_gff_issues_csv(TEST_CHECK_ISSUES, issues_file)


def test_write_issues_to_csv_empty(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.write_issues_to_csv`
    with no values produces a header-only CSV file.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    issues_file = tmpdir.join("issues.tsv")
    check_fasta_gff.write_issues_to_csv([], issues_file)
    check_fasta_gff_issues_csv([], issues_file)


def test_count_issues_empty():
    """
    Test :py:func:`riboviz.check_fasta_gff.count_issues`
    with no values produces a list of issues with 0 counts
    for each issue.
    """
    counts = check_fasta_gff.count_issues({})
    assert len(counts) == len(check_fasta_gff.ISSUE_TYPES)
    for issue in check_fasta_gff.ISSUE_TYPES:
        assert (issue, 0) in counts


def test_count_issues():
    """
    Test :py:func:`riboviz.check_fasta_gff.count_issues`
    with values produces the expected list of issues and
    counts.
    """
    issues = {
        ("s1", "f1", check_fasta_gff.NO_START_CODON, None),
        ("s2", "f2", check_fasta_gff.NO_STOP_CODON, None),
        ("s3", "f3", check_fasta_gff.NO_STOP_CODON, None),
        ("s4", "f4", check_fasta_gff.NO_START_CODON, None),
        ("s5", "f5", check_fasta_gff.NO_START_CODON, None)
    }
    counts = check_fasta_gff.count_issues(issues)
    assert len(counts) == len(check_fasta_gff.ISSUE_TYPES)

    zero_issues = list(check_fasta_gff.ISSUE_TYPES)
    zero_issues.remove(check_fasta_gff.NO_START_CODON)
    zero_issues.remove(check_fasta_gff.NO_STOP_CODON)
    for issue in zero_issues:
        assert (issue, 0) in counts
    # counts should be in sorted order.
    assert counts[0] == (check_fasta_gff.NO_START_CODON, 3)
    assert counts[1] == (check_fasta_gff.NO_STOP_CODON, 2)


def test_check_fasta_gff_no_such_fasta_file(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff`
    with an empty FASTA file and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) raises an exception.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    issues_file = tmpdir.join("issues.tsv")
    with pytest.raises(FileNotFoundError):
        check_fasta_gff.check_fasta_gff("nosuch.fasta",
                                        TEST_GFF_CHECK_FILE,
                                        issues_file)


def test_check_fasta_gff_no_such_gff_file(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and
    a non-existent GFF file raises an exception.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    issues_file = tmpdir.join("issues.tsv")
    with pytest.raises(FileNotFoundError):
        check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                        "nosuch.gff",
                                        issues_file)


def test_check_fasta_gff_empty_fasta_file(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff`
    with an empty FASTA file and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) raises an exception.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    # Cast to str to avoid "'LocalPath' object is not subscriptable"
    # error from gff library.
    fasta_file = str(tmpdir.join("fasta.fasta"))
    open(fasta_file, 'a').close()
    issues_file = tmpdir.join("issues.tsv")
    with pytest.raises(FastaIndexingError):
        check_fasta_gff.check_fasta_gff(fasta_file,
                                        TEST_GFF_CHECK_FILE,
                                        issues_file)


def test_check_fasta_gff_empty_gff_file(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff`
    with FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and
    an empty GFF file raises an exception.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    # Cast to str to avoid "'LocalPath' object is not subscriptable"
    # error from gff library.
    gff_file = str(tmpdir.join("gff.gff"))
    open(gff_file, 'a').close()
    issues_file = tmpdir.join("issues.tsv")
    with pytest.raises(ValueError):
        check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                        gff_file,
                                        issues_file)


def test_check_fasta_gff(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with
    FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and validate the TSV file
    output.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    issues_file = tmpdir.join("issues.tsv")
    check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                    TEST_GFF_CHECK_FILE,
                                    issues_file)
    check_fasta_gff_issues_csv(TEST_CHECK_ISSUES, issues_file)


def test_check_fasta_gff_stdout(tmpdir, capsys):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with
    FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and validate the standard
    output.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    :param capsys: Capture standard output (pytest built-in fixture)
    :type capsys: _pytest.capture.CaptureFixture
    """
    issues_file = tmpdir.join("issues.tsv")
    check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                    TEST_GFF_CHECK_FILE,
                                    issues_file,
                                    is_verbose=True)
    lines = capsys.readouterr()[0].split("\n")
    expected_config = [
        "Configuration:",
        "{}\t{}".format(check_fasta_gff.FASTA_FILE, TEST_FASTA_CHECK_FILE),
        "{}\t{}".format(check_fasta_gff.GFF_FILE, TEST_GFF_CHECK_FILE),
        "{}\t{}".format(check_fasta_gff.START_CODONS, [START_CODON]),
        ""]
    assert lines[0:len(expected_config)] == expected_config,\
        "Unexpected Configuration header"
    lines = lines[len(expected_config):]
    expected_metadata = [
        "Metadata:",
        "{}\t{}".format(check_fasta_gff.NUM_SEQUENCES, TEST_NUM_SEQUENCES),
        "{}\t{}".format(check_fasta_gff.NUM_FEATURES, TEST_NUM_FEATURES),
        "{}\t{}".format(check_fasta_gff.NUM_CDS_FEATURES,
                        TEST_NUM_CDS_FEATURES),
        ""]
    assert lines[0:len(expected_metadata)] == expected_metadata,\
        "Unexpected Metadata header"
    lines = lines[len(expected_metadata):]
    expected_summary = ["Issue summary:", "{}\t{}".format("Issue", "Count")]
    assert lines[0:len(expected_summary)] == expected_summary,\
        "Unexpected Issue summary header"
    lines = lines[len(expected_summary):]
    counts = check_fasta_gff.count_issues(TEST_CHECK_ISSUES)
    expected_counts = [
        "{}\t{}".format(issue, count) for issue, count in counts]
    actual_counts = lines[0:len(expected_counts)]
    # Sort so test passes if same counts are reported but in
    # different order.
    expected_counts.sort()
    actual_counts.sort()
    assert actual_counts == expected_counts, "Unexpected counts"
    lines = lines[len(expected_counts):]
    line = lines.pop(0)
    assert line == "", "Missing blank line"
    line = lines.pop(0)
    assert line == "Issue details:", "Unexpected Issue details header"
    expected_details = [
        check_fasta_gff.ISSUE_FORMATS[issue_type].format(
            sequence=seq_id,
            feature=feature_id,
            data=issue_data)
        for seq_id, feature_id, issue_type, issue_data in TEST_CHECK_ISSUES]
    actual_details = lines[0:len(expected_details)]
    # Sort so test passes if same details are reported but in
    # different order.
    expected_details.sort()
    actual_details.sort()
    assert actual_details == expected_details, "Unexpected details"
    lines = lines[len(expected_details):]
    assert lines == [""], "Missing blank line"


def test_check_fasta_gff_use_feature_name_true(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with
    FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and ``use_feature_name=True``
    and validate the TSV file output.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    issues_file = tmpdir.join("issues.tsv")
    check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                    TEST_GFF_CHECK_FILE,
                                    issues_file,
                                    use_feature_name=True)
    # Create expected results when use_feature_name=True.
    test_check_issues = TEST_CHECK_ISSUES.copy()
    test_check_issues.remove(TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE)
    test_check_issues.append(
        (TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[0],
         TEST_NO_ATG_START_ID_NAME_ATTR_NAME,
         TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[2],
         TEST_NO_ATG_START_ID_NAME_ATTR_ISSUE[3]))
    check_fasta_gff_issues_csv(test_check_issues, issues_file)


def test_check_fasta_gff_feature_format(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with
    FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and custom
    ``cds_feature_format`` and validate the TSV file output.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    issues_file = tmpdir.join("issues.tsv")
    cds_feature_format = "{}-Custom"
    check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                    TEST_GFF_CHECK_FILE,
                                    issues_file,
                                    feature_format=cds_feature_format)
    # Create expected results for custom cds_feature_format.
    test_check_issues = TEST_CHECK_ISSUES.copy()
    test_check_issues.remove(TEST_NO_ID_NAME_ATTR_ISSUE)
    test_check_issues.append(
        (TEST_NO_ID_NAME_ATTR_ISSUE[0],
         cds_feature_format.format(TEST_NO_ID_NAME_ATTR_PREFIX),
         TEST_NO_ID_NAME_ATTR_ISSUE[2],
         TEST_NO_ID_NAME_ATTR_ISSUE[3]))
    check_fasta_gff_issues_csv(test_check_issues, issues_file)


def test_check_fasta_gff_start_codons(tmpdir):
    """
    Test :py:func:`riboviz.check_fasta_gff.check_fasta_gff` with
    FASTA file (:py:const:`TEST_FASTA_CHECK_FILE`) and GFF file
    (:py:const:`TEST_GFF_CHECK_FILE`) and custom
    ``start_codons`` and validate the TSV file output.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    """
    issues_file = tmpdir.join("issues.tsv")
    check_fasta_gff.check_fasta_gff(TEST_FASTA_CHECK_FILE,
                                    TEST_GFF_CHECK_FILE,
                                    issues_file,
                                    start_codons=TEST_START_CODONS)
    # TEST_START_CODONS includes all start codons in test data so
    # expect no NO_START_CODON issues in results.
    test_check_issues = [(s, f, i, d)
                         for (s, f, i, d) in TEST_CHECK_ISSUES
                         if i != check_fasta_gff.NO_START_CODON]
    check_fasta_gff_issues_csv(test_check_issues, issues_file)
