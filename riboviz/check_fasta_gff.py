"""
Functions to check FASTA and GFF files for compatibility.
"""
import bisect
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
import pandas as pd
from riboviz.fasta_gff import CDS_FEATURE_FORMAT
from riboviz.get_cds_codons import sequence_to_codons
from riboviz import provenance

SEQUENCE = "Sequence"
""" FASTA-GFF compatibility column name (sequence ID). """
FEATURE = "Feature"
""" FASTA-GFF compatibility column name (feature ID). """
SEQUENCE_WILDCARD = "*"
""" Sequence ID value for issue that applies to multiple sequences. """
ISSUE_TYPE = "Issue"
""" FASTA-GFF compatibility column name (issue type). """
ISSUE_DATA = "Data"
""" FASTA-GFF compatibility column name (issue data). """
INCOMPLETE_FEATURE = "IncompleteFeature"
""" FASTA-GFF compatibility issue column value. """
NO_ATG_START_CODON = "NoATGStartCodon"
""" FASTA-GFF compatibility issue column value. """
NO_STOP_CODON = "NoStopCodon"
""" FASTA-GFF compatibility issue column value. """
INTERNAL_STOP_CODON = "InternalStopCodon"
""" FASTA-GFF compatibility issue column value. """
NO_ID_NAME = "NoIdName"
""" FASTA-GFF compatibility issue column value. """
DUPLICATE_FEATURE_ID = "DuplicateFeatureId"
""" FASTA-GFF compatibility issue column value. """
DUPLICATE_FEATURE_IDS = "DuplicateFeatureIds"
""" FASTA-GFF compatibility issue column value. """
MULTIPLE_CDS = "MultipleCDS"
""" FASTA-GFF compatibility issue column value. """
SEQUENCE_NOT_IN_FASTA = "SequenceNotInFASTA"
""" FASTA-GFF compatibility issue column value. """
SEQUENCE_NOT_IN_GFF = "SequenceNotInGFF"
""" FASTA-GFF compatibility issue column value. """
ISSUE_FORMATS = {
    INCOMPLETE_FEATURE: "Sequence {sequence} feature {feature} has length not divisible by 3",
    NO_ATG_START_CODON: "Sequence {sequence} feature {feature} doesn't start with ATG",
    NO_STOP_CODON: "Sequence {sequence} feature {feature} doesn't stop at end",
    INTERNAL_STOP_CODON: "Sequence {sequence} feature {feature} has internal STOP",
    NO_ID_NAME: "Sequence {sequence} feature {feature} has no 'ID' or 'Name' attribute",
    DUPLICATE_FEATURE_ID: "Sequence {sequence} has non-unique 'ID' attribute {feature}",
    MULTIPLE_CDS: "Sequence {sequence} has multiple CDS",
    SEQUENCE_NOT_IN_FASTA: "Sequence {sequence} in GFF file is not in FASTA file",
    SEQUENCE_NOT_IN_GFF: "Sequence {sequence} in FASTA file is not in GFF file",
    DUPLICATE_FEATURE_IDS: "Non-unique 'ID' attribute {feature} ({data} occurrences)"
}
""" Format strings for printing compatibility issues. """


def get_fasta_sequence_ids(fasta):
    """
    Get a list of the unique IDs of sequences in a FASTA file.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :return: Unique sequence IDs
    :rtype: set(str or unicode)
    """
    seq_ids = set()
    with open(fasta, "r") as f:
        # 'fasta' is https://biopython.org/wiki/SeqIO file type.
        for record in SeqIO.parse(f, "fasta"):
            seq_ids.add(record.id)
    return seq_ids


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

    * :py:const:`INCOMPLETE_FEATURE`: The CDS has length not
      divisible by 3.
    * :py:const:`NO_ATG_START_CODON`: The beginning of a CDS
      does not have a start codon (``ATG``, translates to ``M``)
    * :py:const:`NO_STOP_CODON`: The end codon of the CDS is not
      a stop codon (``TAG``, ``TGA``, ``TAA``, translates to ``*``).
    * :py:const:`INTERNAL_STOP_CODON`: There are stop codons
      internal to the CDS.
    * :py:const:`NO_ID_NAME`: The CDS has no ``ID`` or ``Name``
      attribute.
    * :py:const:`DUPLICATE_FEATURE_ID`: The CDS has non-unique
      ``ID`` attribute.
    * :py:const:`MULTIPLE_CDS`: Sequence has multiple CDS.
    * :py:const:`SEQUENCE_NOT_IN_FASTA`: Sequence in GFF
      file is not in FASTA file.
    * :py:const:`SEQUENCE_NOT_IN_GFF`: Sequence in FASTA
      file is not in GFF file.
    * :py:const:`DUPLICATE_FEATURE_IDS`: CDSs have non-unique
      ``ID`` attributes. This summarises the count of all CDSs that
     share a common ``ID`` attribute.

    Issue data is supplementary data relating to the issue e.g. for
    :py:const:`DUPLICATE_FEATURE_ID` the number of features with
    the same feature ID. This is a string or int or float or ``None``.

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
    :return: List of unique sequence IDs in GFF file and list \
    of issues for sequences and features.
    :rtype: list(tuple(str or unicode, str or unicode, str or unicode, *))
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
        if feature.seqid not in sequence_features:
            sequence_features[feature.seqid] = 0
        sequence_features[feature.seqid] += 1

        feature_id = None
        if "ID" in feature.attributes:
            feature_id = feature.attributes["ID"][0].strip()
            if feature_id in feature_ids:
                feature_ids[feature_id].append(feature.seqid)
            else:
                feature_ids[feature_id] = [feature.seqid]
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
                           NO_ID_NAME, None))
        try:
            sequence = feature.sequence(fasta)
        except KeyError as e:
            issues.append((feature.seqid, '', SEQUENCE_NOT_IN_FASTA,
                           None))
            continue
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
                           INCOMPLETE_FEATURE, None))
            sequence += ("N" * (3 - seq_len_remainder))

        sequence_codons = sequence_to_codons(sequence)
        if sequence_codons[0] != "ATG":
            issues.append((feature.seqid, feature_id_name,
                           NO_ATG_START_CODON, None))
        if not sequence_codons[-1] in ["TAA", "TAG", "TGA"]:
            issues.append((feature.seqid, feature_id_name,
                           NO_STOP_CODON, None))
        if any([codon in ["TAA", "TAG", "TGA"]
                for codon in sequence_codons[:-1]]):
            issues.append((feature.seqid, feature_id_name,
                           INTERNAL_STOP_CODON, None))

    for sequence, num_features in list(sequence_features.items()):
        if num_features > 1:
            # Insert issue, keeping entries grouped by sequence ID.
            bisect.insort(issues, (sequence, '', MULTIPLE_CDS, None))

    for feature_id, seq_ids in feature_ids.items():
        if len(seq_ids) > 1:
            issues.append((SEQUENCE_WILDCARD,
                           feature_id,
                           DUPLICATE_FEATURE_IDS,
                           len(seq_ids)))
            for seq_id in seq_ids:
                # Insert issue, keeping entries grouped by sequence ID.
                bisect.insort(issues, (seq_id, feature_id,
                                       DUPLICATE_FEATURE_ID, None))

    gff_seq_ids = set(sequence_features.keys())
    fasta_seq_ids = get_fasta_sequence_ids(fasta)
    fasta_only_seq_ids = fasta_seq_ids - gff_seq_ids
    for seq_id in fasta_only_seq_ids:
        issues.append((seq_id, '', SEQUENCE_NOT_IN_GFF, None))
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
    str or unicode, *))
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
        if issue in ISSUE_FORMATS:
            print(ISSUE_FORMATS[issue].format(
                sequence=row[SEQUENCE],
                feature=row[FEATURE],
                data=row[ISSUE_DATA]))
    provenance.write_provenance_header(__file__, issues_file)
    issues_df[list(issues_df.columns)].to_csv(
        issues_file,
        mode='a',
        sep=delimiter,
        index=False)
