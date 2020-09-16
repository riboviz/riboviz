"""
FASTA and GFF-related functions.
"""
import bisect
import warnings
from Bio.Seq import Seq
import gffutils
import pandas as pd
from riboviz import provenance


GENE = "Gene"
""" Codon positions column name (gene name). """
POS = "Pos"
"""
Codon positions column name (codon position in coding sequence,
1-indexed by codon).
"""
CODON = "Codon"
""" Codon positions column name (codon). """
CDS_FEATURE_FORMAT = "{}_CDS"
"""
CDS feature name format for CDS features which do not define ``ID``
or ``Name`` attributes.
"""

SEQUENCE = "Sequence"
""" FASTA-GFF compatibility column name (sequence ID). """
FEATURE = "Feature"
""" FASTA-GFF compatibility column name (feature ID). """
UNDEFINED_FEATURE = "Undefined"
"""
FASTA-GFF compatibility feature column value for features with
undefined names
"""
ISSUE = "Issue"
""" FASTA-GFF compatibility column name (issue). """
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


def get_fasta_gff_cds_issues(fasta, gff):
    """
    Check FASTA and GFF files for compatibility and return a list of
    issues for relating to coding sequences, ``CDS``, features. A list
    of tuples of form (sequence ID, feature ID ('' if not applicable
    to the issue), issue) is returned.

    The feature ID is one of:

    * ``''`` if the issue relates to the sequence, not the feature.
    * Value of ``ID`` attribute for feature, if defined.
    * Value of ``Name`` attribute for feature, if defined, and if
      ``ID`` is undefined.
    * ``Undefined`` otherwise.

    Issues are one of:

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

    Some unusual genes (e.g. frame shifts) might have these issues.

    If a feature has no issues then it has a no entry in the
    issues dictionary.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :return: List of issues for sequences and features.
    :rtype: list(tuple(str or unicode, str or unicode, str or \
    unicode))
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
    feature_ids = set()
    sequence_features = {}
    for feature in gffdb.features_of_type('CDS'):
        feature_id = None
        if "ID" in feature.attributes:
            feature_id = feature.attributes["ID"][0].strip()
            if feature_id in feature_ids:
                issues.append((feature.seqid,
                               feature_id,
                               ISSUE_DUPLICATE_FEATURE_ID))
            else:
                feature_ids.add(feature_id)

        feature_name = None
        if "Name" in feature.attributes:
            feature_name = feature.attributes["Name"][0].strip()
        if not feature_id and not feature_name:
            issues.append((feature.seqid,
                           UNDEFINED_FEATURE,
                           ISSUE_NO_ID_NAME))
        if feature_id:
            feature_id_name = feature_id
        elif feature_name:
            feature_id_name = feature_name
        else:
            feature_id_name = UNDEFINED_FEATURE

        try:
            sequence = feature.sequence(fasta)
        except KeyError as e:
            issues.append((feature.seqid, '', ISSUE_MISSING_SEQUENCE))
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
                           ISSUE_INCOMPLETE))
            sequence += ("N" * (3 - seq_len_remainder))

        translation = Seq(sequence).translate()
        if translation[0] != "M":
            issues.append((feature.seqid, feature_id_name,
                           ISSUE_NO_ATG_START))

        if translation[-1] != "*":
            issues.append((feature.seqid, feature_id_name,
                           ISSUE_NO_STOP))

        if any([L == "*" for L in translation[:-1]]):
            issues.append((feature.seqid, feature_id_name,
                           ISSUE_INTERNAL_STOP))

        if feature.seqid not in sequence_features:
            sequence_features[feature.seqid] = 0
        sequence_features[feature.seqid] += 1

    for sequence, num_features in list(sequence_features.items()):
        if num_features > 1:
            # Insert issue into list so entries remain grouped
            # by sequence.
            bisect.insort(issues, (sequence, '', ISSUE_MULTIPLE_CDS))
    return issues


def fasta_gff_issues_to_df(issues):
    """
    Given dictionary of the issues for features, keyed by feature
    name, return a Pandas data frame with the issues. The data frame
    is sorted by sequence ID.

    The data frame has columns:

    * :py:const:`SEQUENCE`: sequence ID.
    * :py:const:`FEATURE`: feature ID.
    * :py:const:`ISSUE`: issue.

    :param issues: List of issues for sequences and features.
    :type issues: list(tuple(str or unicode, str or unicode, \
    str or unicode))
    :return: data frame
    :rtype: pandas.core.frame.DataFrame
    """
    issues_list = []
    for (sequence_id, feature_id, issue) in issues:
        issues_list.append({SEQUENCE: sequence_id,
                            FEATURE: feature_id,
                            ISSUE: issue})
    # Create empty DataFrame so if issues and
    # issues_list are empty we still have an empty DataFrame
    # with the column names.
    df = pd.DataFrame(columns=[SEQUENCE, FEATURE, ISSUE])
    df = df.append(pd.DataFrame(issues_list),
                   ignore_index=True)
    return df


def check_fasta_gff(fasta, gff, issues_file, delimiter="\t"):
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
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    issues = get_fasta_gff_cds_issues(fasta, gff)
    issues_df = fasta_gff_issues_to_df(issues)
    for _, row in issues_df.iterrows():
        issue = row[ISSUE]
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


def sequence_to_codons(sequence):
    """
    Given a sequence, split into a list of codons.

    :param sequence: Sequence
    :type sequence: str or unicode
    :return: list of codons
    :rtype: list(str or unicode)
    """
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    # Validate split was done correctly.
    assert "".join(codons) == sequence
    return codons


def get_feature_id(feature):
    """
    Get the name of a GFF feature. If there is no feature name (no
    ``ID`` or ``Name``) attribute defined in the feature then ``None``
    is returned. If both are defined then ``ID`` is returned by
    preference.

    :param feature: GFF feature
    :type feature: gffutils.feature.Feature
    :return: Feature name or ``None``
    :rtype: str or unicode
    """
    name = None
    if "ID" in feature.attributes:
        name_attr = feature.attributes["ID"]
    elif "Name" in feature.attributes:
        name_attr = feature.attributes["Name"]
    else:
        name_attr = []
    if name_attr != []:
        name = name_attr[0].strip()
    return name


def get_cds_from_fasta(feature, fasta):
    """
    Get the sequence of a CDS from a FASTA file using a GFF feature
    for the CDS.

    :param feature: GFF feature for the CDS
    :type feature: gffutils.feature.Feature
    :param fasta: FASTA file
    :type fasta: str or unicode
    :return: sequence
    :rtype: str or unicode
    :raises AssertionError: If sequence has length not divisible by 3
    :raises Exception: Exceptions specific to \
    gffutils.Feature.sequence (these are undocumented in the \
    gffutils documentation - a typical exception that can be thrown \
    is ``KeyError`` which can arise if the GFF file contains \
    information on a sequence that is not in the FASTA file).
    """
    sequence = feature.sequence(fasta)
    assert (len(sequence) % 3) == 0, \
        "Feature {} ({}) has length not divisible by 3".format(
            feature.seqid, sequence)
    return sequence


def get_cds_codons_from_fasta(fasta,
                              gff,
                              exclude_stop_codons=False,
                              feature_format=CDS_FEATURE_FORMAT):
    """
    Using CDS entries within a GFF file, get the codons in each coding
    sequence in the complementary FASTA file.

    A dictionary of the codons for each CDS, keyed by CDS feature
    name, is returned.

    CDSs whose sequences don't have a length divisible by 3 are
    ignored.

    If there is no feature name (no ``ID`` or ``Name``) attribute
    defined in the feature the the ``feature_format`` is used to
    format the sequence ID into a feature name. If both ``ID`` and
    ``Name`` are defined then ``ID`` is used.

    If two or more CDSs for the same sequence have the same feature
    name then the first CDS for that sequence has a feature name, as
    defined above. Subsequent CDSs for that sequence have the feature
    name with with the suffix ``.1``, ``.2`` etc. appended.

    See also :py:func:`get_cds_from_fasta` and
    :py:func:`sequence_to_codons`.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :param exclude_stop_codons: Should stop codons be excluded from \
    codons returned?
    :type exclude_stop_codons: bool
    :param feature_format: Feature name format
    :type feature_format: str or unicode
    :return: Codons for each coding sequence, keyed by feature name
    :rtype: dict(str or unicode -> list(str or unicode))
    :raises Exception: Exceptions specific to gffutils.create_db \
    (these are undocumented in the gffutils documentation)
    """
    gffdb = gffutils.create_db(gff,
                               dbfn='gff.db',
                               force=True,
                               keep_order=True,
                               merge_strategy='merge',
                               sort_attribute_values=True)
    cds_codons = {}
    same_feature_id_count = 0
    for feature in gffdb.features_of_type('CDS'):
        try:
            sequence = get_cds_from_fasta(feature, fasta)
            codons = sequence_to_codons(sequence)
        except Exception as e:
            # Log and continue with other CDSs.
            warnings.warn(str(e))
            continue
        feature_id = get_feature_id(feature)
        if feature_id is None:
            feature_id = feature_format.format(feature.seqid)
        if feature_id not in cds_codons:
            same_feature_id_count = 0
        else:
            same_feature_id_count += 1
            feature_id = "{}.{}".format(feature_id,
                                        same_feature_id_count)
        if exclude_stop_codons:
            codons = codons[:-1]
        cds_codons[feature_id] = codons
    return cds_codons


def feature_codons_to_df(feature_codons):
    """
    Given dictionary of the codons for features, keyed by feature
    name, return a Pandas data frame with the codons, also
    including the position of each codon in its enclosing sequence.

    The data frame has columns:

    * :py:const:`GENE`: feature name.
    * :py:const:`CODON`: codon.
    * :py:const:`POS`: codon position in coding sequence (1-indexed).

    :param feature_codons: Codons for each feature, keyed by feature \
    name
    :type feature_codons: dict(str or unicode -> list(str or unicode))
    :return: data frame
    :rtype: pandas.core.frame.DataFrame
    """
    feature_codons_list = []
    num_feature_codons = 0
    for feature_id, codons in list(feature_codons.items()):
        num_feature_codons += len(codons)
        for pos, codon in zip(range(0, len(codons)), codons):
            feature_codons_list.append({GENE: feature_id,
                                        CODON: codon,
                                        POS: pos + 1})
    # Create empty DataFrame so if feature_codons and
    # feature_codons_list are empty we still have an empty DataFrame
    # with the column names.
    df = pd.DataFrame(columns=[GENE, CODON, POS])
    df = df.append(pd.DataFrame(feature_codons_list),
                   ignore_index=True)
    # Validate number of codons
    assert num_feature_codons == df.shape[0]
    return df


def get_cds_codons_file(fasta,
                        gff,
                        cds_codons_file,
                        exclude_stop_codons=False,
                        feature_format=CDS_FEATURE_FORMAT,
                        delimiter="\t"):
    """
    Using CDS entries within a GFF file, get the codons in each coding
    sequence in the complementary FASTA file.

    A tab-separated values file of the codons for each CDS, keyed by
    CDS feature name, is saved.

    See :py:func:`get_cds_codons_from_fasta`.

    See :py:func:`feature_codons_to_df` for tab-separated values
    file columns.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :param cds_codons_file: Coding sequence codons file
    :type cds_codons_file: str or unicode
    :param exclude_stop_codons: Should stop codons be excluded from \
    codons returned?
    :type exclude_stop_codons: bool
    :param feature_format: Feature name format
    :type feature_format: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    cds_codons = get_cds_codons_from_fasta(fasta,
                                           gff,
                                           exclude_stop_codons,
                                           feature_format)
    cds_codons_df = feature_codons_to_df(cds_codons)
    provenance.write_provenance_header(__file__, cds_codons_file)
    cds_codons_df[list(cds_codons_df.columns)].to_csv(cds_codons_file,
                                                      mode='a',
                                                      sep=delimiter,
                                                      index=False)
