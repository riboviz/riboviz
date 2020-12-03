"""
Functions to check FASTA and GFF files for compatibility.
"""
import bisect
import warnings
from Bio.Seq import Seq
import gffutils
import pandas as pd
from riboviz.fasta_gff import CDS_FEATURE_FORMAT
from riboviz import provenance

SEQUENCE = "Sequence"
""" FASTA-GFF compatibility column name (sequence ID). """
FEATURE = "Feature"
""" FASTA-GFF compatibility column name (feature ID). """
ISSUE_TYPE = "Issue"
""" FASTA-GFF compatibility column name (issue type). """
ISSUE_DATA = "Data"
""" FASTA-GFF compatibility column name (issue data). """
ISSUE_INCOMPLETE = "Incomplete"
""" FASTA-GFF compatibility issue column value. """
ISSUE_NO_ATG_START = "NoATGStart"
""" FASTA-GFF compatibility issue column value. """
ISSUE_NO_STOP = "NoStop"
""" FASTA-GFF compatibility issue column value. """
ISSUE_INTERNAL_STOP = "InternalStop"
""" FASTA-GFF compatibility issue column value. """
ISSUE_NO_ID_NAME = "NoIdName"
""" FASTA-GFF compatibility issue column value. """
ISSUE_DUPLICATE_FEATURE_ID = "DuplicateFeatureId"
""" FASTA-GFF compatibility issue column value. """
ISSUE_MULTIPLE_CDS = "MultipleCDS"
""" FASTA-GFF compatibility issue column value. """
ISSUE_MISSING_SEQUENCE = "MissingSequence"
""" FASTA-GFF compatibility issue column value. """
ISSUE_FEATURE_FORMATS = {
    ISSUE_INCOMPLETE: "Sequence {} feature {} has length not divisible by 3",
    ISSUE_NO_ATG_START: "Sequence {} feature {} doesn't start with ATG",
    ISSUE_NO_STOP: "Sequence {} feature {} doesn't stop at end",
    ISSUE_INTERNAL_STOP: "Sequence {} feature {} has internal STOP",
    ISSUE_NO_ID_NAME: "Sequence {} feature {} has no 'ID' or 'Name' attribute",
    ISSUE_DUPLICATE_FEATURE_ID: "Sequence {} has non-unique 'ID' attribute {}"
}
""" Format strings for printing FASTA-GFF compatibility issues. """
ISSUE_SEQUENCE_FORMATS = {
    ISSUE_MULTIPLE_CDS: "Sequence {} has multiple CDS",
    ISSUE_MISSING_SEQUENCE: "Sequence {} missing in FASTA file"
}
""" Format strings for printing FASTA-GFF compatibility issues. """


def get_fasta_gff_cds_issues(fasta, gff, feature_format=CDS_FEATURE_FORMAT):
    """
    Check FASTA and GFF files for compatibility and return a list of
    issues for relating to coding sequences, ``CDS``, features. A list
    of tuples of form (sequence ID, feature ID ('' if not applicable
    to the issue), issue type, issue data) is returned.

    The feature ID is one of:

    * ``''`` if the issue relates to the sequence, not the feature.
    * Value of ``ID`` attribute for feature, if defined.
    * Value of ``Name`` attribute for feature, if defined, and if
      ``ID`` is undefined.
    * Sequence ID formatted using ``feature_format`` if both ``ID``
      and ``Name`` are undefined.

    Issue types are one of:

    * :py:const:`ISSUE_INCOMPLETE`: The CDS has length not divisible by
      3.
    * :py:const:`ISSUE_NO_ATG_START`: The beginning of a CDS does not
      have a start codon (``ATG``, translates to ``M``)
    * :py:const:`ISSUE_NO_STOP`: The end codon of the CDS is not a
      stop codon (``TAG``, ``TGA``, ``TAA``, translates to ``*``).
    * :py:const:`ISSUE_INTERNAL_STOP`: There are stop codons internal
      to the CDS.
    * :py:const:`ISSUE_NO_ID_NAME`: The CDS has no ``ID`` or ``Name``
      attribute.
    * :py:const:`ISSUE_DUPLICATE_FEATURE_ID`: The CDS has non-unique
      ``ID`` attribute.
    * :py:const:`ISSUE_MULTIPLE_CDS`: Sequence has multiple CDS.
    * :py:const:`ISSUE_MISSING_SEQUENCE`: Sequence missing in FASTA
      file.

    Issue data is supplementary data relating to the issue e.g. for
    :py:const:`ISSUE_DUPLICATE_FEATURE_ID` the number of features with
    the same feature ID. This is a dictionary or ``None``.

    Some unusual genes (e.g. frame shifts) might have these issues.

    If a feature has no issues then it has no entry in the
    issues dictionary.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :param feature_format: Feature name format for features which \
    do not define ``ID``  or ``Name`` attributes. This format is \
    applied to the sequence ID to create a feature name.
    :type feature_format: str or unicode
    :return: List of issues for sequences and features.
    :rtype: list(tuple(str or unicode, str or unicode, str or \
    unicode, dict or None))
    """
    gffdb = gffutils.create_db(gff,
                               dbfn='gff.db',
                               force=True,
                               keep_order=True,
                               merge_strategy='merge',
                               sort_attribute_values=True)
    issues = []
    # Track IDs of features encountered. Each ID must be unique within
    # a GFF file. See http://gmod.org/wiki/GFF3.
    feature_ids = {}
    sequence_features = {}
    for feature in gffdb.features_of_type('CDS'):
        feature_id = None
        if "ID" in feature.attributes:
            feature_id = feature.attributes["ID"][0].strip()
            if feature_id in feature_ids:
                feature_ids[feature_id] += 1
            else:
                feature_ids[feature_id] = 1
        feature_name = None
        if "Name" in feature.attributes:
            feature_name = feature.attributes["Name"][0].strip()
        if feature_id:
            feature_id_name = feature_id
        elif feature_name:
            feature_id_name = feature_name
        else:
            feature_id = feature_format.format(feature.seqid)
            issues.append((feature.seqid, feature_id,
                           ISSUE_NO_ID_NAME, None))
        try:
            sequence = feature.sequence(fasta)
        except KeyError as e:
            issues.append((feature.seqid, '', ISSUE_MISSING_SEQUENCE,
                           None))
        except Exception as e:
            # Log and continue with other CDSs. A typical exception
            # that can be thrown by gffutils.Feature.sequence is
            # KeyError. This can arise if the GFF file contains
            # information on a sequence that is not in the FASTA
            # file.
            warnings.warn(str(e))
            continue

        seq_len_remainder = len(sequence) % 3
        if seq_len_remainder != 0:
            issues.append((feature.seqid, feature_id_name,
                           ISSUE_INCOMPLETE, None))
            sequence += ("N" * (3 - seq_len_remainder))

        translation = Seq(sequence).translate()
        if translation[0] != "M":
            issues.append((feature.seqid, feature_id_name,
                           ISSUE_NO_ATG_START, None))

        if translation[-1] != "*":
            issues.append((feature.seqid, feature_id_name,
                           ISSUE_NO_STOP, None))

        if any([L == "*" for L in translation[:-1]]):
            issues.append((feature.seqid, feature_id_name,
                           ISSUE_INTERNAL_STOP, None))

        if feature.seqid not in sequence_features:
            sequence_features[feature.seqid] = 0
        sequence_features[feature.seqid] += 1
    for feature_id, count in feature_ids.items():
        if count > 1:
            issues.append(('', feature_id, ISSUE_DUPLICATE_FEATURE_ID,
                           {'count': count}))

    for sequence, num_features in list(sequence_features.items()):
        if num_features > 1:
            # Insert issue into list so entries remain grouped
            # by sequence.
            bisect.insort(issues, (sequence, '', ISSUE_MULTIPLE_CDS, None))
    return issues


def fasta_gff_issues_to_df(issues):
    """
    Given dictionary of the issues for features, keyed by feature
    name, return a Pandas data frame with the issues. The data frame
    is sorted by sequence ID.

    The data frame has columns:

    * :py:const:`SEQUENCE`: sequence ID.
    * :py:const:`FEATURE`: feature ID.
    * :py:const:`ISSUE_TYPE`: issue type.
    * :py:const:`ISSUE_DATA`: issue data or ``None``.

    :param issues: List of issues for sequences and features.
    :type issues: list(tuple(str or unicode, str or unicode, \
    str or unicode, dict))
    :return: data frame
    :rtype: pandas.core.frame.DataFrame
    """
    issues_list = []
    for (sequence_id, feature_id, issue_type, issue_data) in issues:
        issues_list.append({SEQUENCE: sequence_id,
                            FEATURE: feature_id,
                            ISSUE_TYPE: issue_type,
                            ISSUE_DATA: issue_data})
    # Create empty DataFrame so if issues and
    # issues_list are empty we still have an empty DataFrame
    # with the column names.
    df = pd.DataFrame(columns=[SEQUENCE, FEATURE, ISSUE_TYPE, ISSUE_DATA])
    df = df.append(pd.DataFrame(issues_list),
                   ignore_index=True)
    return df


def check_fasta_gff(fasta, gff, issues_file,
                    feature_format=CDS_FEATURE_FORMAT,
                    delimiter="\t"):
    """
    Check FASTA and GFF files for compatibility and both print and
    save a list of issues for each coding sequence, ``CDS``,
    feature.

    A tab-separated values file of the issues for each CDS, keyed by
    sequence ID and CDS feature name, is saved.

    See :py:func:`get_fasta_gff_cds_issues` for feature IDs and
    issues.

    See :py:func:`fasta_gff_issues_to_df` for tab-separated values
    file columns.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :param fasta: FASTA file
    :param issues_file: Feature issues file
    :type issues_file: str or unicode
    :param feature_format: Feature name format for features which \
    do not define ``ID``  or ``Name`` attributes. This format is \
    applied to the sequence ID to create a feature name.
    :type feature_format: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    issues = get_fasta_gff_cds_issues(fasta, gff, feature_format)
    issues_df = fasta_gff_issues_to_df(issues)
    for _, row in issues_df.iterrows():
        issue = row[ISSUE_TYPE]
        if issue in ISSUE_FEATURE_FORMATS:
            print(ISSUE_FEATURE_FORMATS[issue].format(row[SEQUENCE],
                                                      row[FEATURE]))
        elif issue in ISSUE_SEQUENCE_FORMATS:
            print(ISSUE_SEQUENCE_FORMATS[issue].format(row[SEQUENCE]))
    provenance.write_provenance_header(__file__, issues_file)
    issues_df[list(issues_df.columns)].to_csv(
        issues_file,
        mode='a',
        sep=delimiter,
        index=False)
