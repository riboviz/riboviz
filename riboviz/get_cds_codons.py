"""
Functions use CDS entries within a GFF file to get the codons from
each coding sequence in a complementary FASTA file.
"""
import warnings
import gffutils
import pandas as pd
from riboviz import provenance
from riboviz.fasta_gff import CDS_FEATURE_FORMAT


GENE = "Gene"
""" Codon positions column name (gene name). """
POS = "Pos"
"""
Codon positions column name (codon position in coding sequence,
1-indexed by codon).
"""
CODON = "Codon"
""" Codon positions column name (codon). """


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
                              cds_feature_format=CDS_FEATURE_FORMAT):
    """
    Using CDS entries within a GFF file, get the codons in each coding
    sequence in the complementary FASTA file.

    A dictionary of the codons for each CDS, keyed by CDS feature
    name, is returned.

    CDSs whose sequences don't have a length divisible by 3 are
    ignored.

    If there is no feature name (no ``ID`` or ``Name``) attribute
    defined in the feature the the ``cds_feature_format`` is used to
    format the sequence ID into a feature name.

    If both ``ID`` and ``Name`` are defined then ``ID`` is used.

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
    :param cds_feature_format: CDS feature name format for CDS \
    features which do not define ``ID``  or ``Name`` attributes. \
    This format is applied to the sequence ID to create a \
    feature name.
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
            feature_id = cds_feature_format.format(feature.seqid)
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
                        cds_feature_format=CDS_FEATURE_FORMAT,
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
    :param cds_feature_format: CDS feature name format for CDS \
    features which do not define ``ID``  or ``Name`` attributes. \
    This format is applied to the sequence ID to create a \
    feature name.
    :type cds_feature_format: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    cds_codons = get_cds_codons_from_fasta(fasta,
                                           gff,
                                           exclude_stop_codons,
                                           cds_feature_format)
    cds_codons_df = feature_codons_to_df(cds_codons)
    provenance.write_provenance_header(__file__, cds_codons_file)
    cds_codons_df[list(cds_codons_df.columns)].to_csv(cds_codons_file,
                                                      mode='a',
                                                      sep=delimiter,
                                                      index=False)
