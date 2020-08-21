"""
FASTA and GFF-related functions.
"""
import warnings
from Bio.Seq import Seq
import gffutils
import pandas as pd
from riboviz import provenance


GENE = "Gene"
""" Column name (gene name). """
POS = "Pos"
"""
Column name (codon position in coding sequence, 1-indexed by codon).
"""
CODON = "Codon"
""" Column name (codon). """


def check_fasta_gff(fasta, gff):
    """
    Check FASTA and GFF files for compatibility. Check that:

    * The beginning of every CDS is a start codon (ATG; translates to
      ``M``).
    * The stop of every CDS is a stop codon (TAG, TGA, TAA; translates
      to ``*``)
    * There are no stop codons internal to the CDS.

    Some unusual genes (e.g. frameshifts) might not have this.

    Information is currently printed to standard output.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    """
    print(("Checking fasta file " + fasta))
    print(("with gff file " + gff))
    gffdb = gffutils.create_db(gff,
                               dbfn='test.db',
                               force=True,
                               keep_order=True,
                               merge_strategy='merge',
                               sort_attribute_values=True)
    for cds_coord in gffdb.features_of_type('CDS'):
        try:
            cds_seq = cds_coord.sequence(fasta)
        except Exception as e:
            # Log and continue with other CDSs. A typical exception
            # that can be thrown by gffutils.Feature.sequence is
            # KeyError. This can arise if the GFF file contains
            # information on a sequence that is not in the FASTA
            # file.
            warnings.warn(str(e))
            continue
        cds_len_remainder = len(cds_seq) % 3
        if cds_len_remainder != 0:
            warnings.warn(
                cds_coord.seqid + " has length not divisible by 3")
            cds_seq += ("N" * (3 - cds_len_remainder))

        cds_trans = Seq(cds_seq).translate()

        if cds_trans[0] != "M":
            print((cds_coord.seqid + " doesn't start with ATG."))
        if cds_trans[-1] != "*":
            print((cds_coord.seqid + " doesn't stop at end."))
        if any([L == "*" for L in cds_trans[:-1]]):
            print((cds_coord.seqid + " has internal STOP."))


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


def get_feature_name(feature, feature_format):
    """
    Get the name of a GFF feature.

    If there is no feature name (no ``Name`` or ``ID``) attribute
    defined in the feature the the ``feature_format`` is used to
    format the sequence ID into a feature name. If both ``Name`` and
    ``ID`` are defined then ``Name`` is used.

    :param feature: GFF feature
    :type feature: gffutils.feature.Feature
    :param feature_format: Feature name format
    :type feature_format: str or unicode
    :return: Feature name
    :rtype: str or unicode
    """
    if "Name" in feature.attributes:
        name_attr = feature.attributes["Name"]
    elif "ID" in feature.attributes:
        name_attr = feature.attributes["ID"]
    else:
        name_attr = []
    if name_attr != []:
        name = name_attr[0].strip()
    else:
        name = feature_format.format(feature.seqid)
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
    :raises Exception: Exceptions specific to
    gffutils.Feature.sequence (these are undocumented in the gffutils
    documentation). A typical exception that can be thrown is
    KeyError. This can arise if the GFF file contains information on a
    sequence that is not in the FASTA file.
    """
    sequence = feature.sequence(fasta)
    assert (len(sequence) % 3) == 0, \
        "Feature {} ({}) has length not divisible by 3".format(
            feature.seqid, sequence)
    return sequence


def get_cds_codons_from_fasta(fasta,
                              gff,
                              feature_format="{}_CDS"):
    """
    Using CDS entries within a GFF file, get the codons in each coding
    sequence in the complementary FASTA file.

    A dictionary of the codons for each CDS, keyed by CDS feature
    name, is returned.

    CDSs whose sequences don't have a length divisible by 3 are
    ignored.

    If there is no feature name (no ``Name`` or ``ID``) attribute
    defined in the feature the the ``feature_format`` is used to
    format the sequence ID into a feature name. If both ``Name`` and
    ``ID`` are defined then ``Name`` is used.

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
    :param feature_format: Feature name format
    :type feature_format: str or unicode
    :return: Codons for each coding sequence, keyed by feature name
    :rtype: dict(str or unicode -> list(str or unicode))
    :raises Exception: Exceptions specific to gffutils.create_db
    (these are undocumented in the gffutils documentation)
    """
    gffdb = gffutils.create_db(gff,
                               dbfn='gff.db',
                               force=True,
                               keep_order=True,
                               merge_strategy='merge',
                               sort_attribute_values=True)
    cds_codons = {}
    same_feature_name_count = 0
    for feature in gffdb.features_of_type('CDS'):
        try:
            sequence = get_cds_from_fasta(feature, fasta)
            codons = sequence_to_codons(sequence)
        except Exception as e:
            # Log and continue with other CDSs.
            warnings.warn(str(e))
            continue
        feature_name = get_feature_name(feature, feature_format)
        if feature_name not in cds_codons:
            same_feature_name_count = 0
        else:
            same_feature_name_count += 1
            feature_name = "{}.{}".format(feature_name,
                                          same_feature_name_count)
        cds_codons[feature_name] = codons
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

    :param feature_codons: Codons for each feature, keyed by feature
    name
    :type feature_codons: dict(str or unicode -> list(str or unicode))
    :return: data frame
    :rtype: pandas.core.frame.DataFrame
    """
    feature_codons_list = []
    num_feature_codons = 0
    for feature_name, codons in list(feature_codons.items()):
        num_feature_codons += len(codons)
        for pos, codon in zip(range(0, len(codons)), codons):
            feature_codons_list.append({GENE: feature_name,
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


def extract_cds_codons(fasta,
                       gff,
                       cds_codons_file,
                       delimiter="\t"):
    """
    Using CDS entries within a GFF file, get the codons in each coding
    sequence in the complementary FASTA file.

    A tab-separated values file of the codons for each CDS, keyed by
    CDS feature name, is saved.

    The tab-separated values file has columns:

    * :py:const:`GENE`: feature name.
    * :py:const:`CODON`: codon.
    * :py:const:`POS`: codon position in coding sequence (1-indexed).

    See also :py:func:`get_all_cds_codons_from_fasta` and
    :py:func:`feature_codons_to_df`.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :params cds_codons_file: Coding sequence codons file
    :type cds_codons_file: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    cds_codons = get_cds_codons_from_fasta(fasta, gff)
    cds_codons_df = feature_codons_to_df(cds_codons)
    provenance.write_provenance_header(__file__, cds_codons_file)
    cds_codons_df[list(cds_codons_df.columns)].to_csv(cds_codons_file,
                                                      mode='a',
                                                      sep=delimiter,
                                                      index=False)
