"""
Functions use CDS entries within a GFF file to get the codons from
each coding sequence in a complementary FASTA file.
"""
import csv
import os
import warnings
import gffutils
from pyfaidx import Fasta
from riboviz import provenance
from riboviz.fasta_gff import CDS_FEATURE_FORMAT


GENE = "Gene"
""" Codon positions column name (gene name). """
POS = "PosCodon"
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


def get_feature_id(feature, use_feature_name=False):
    """
    Get the name of a GFF feature. If there is no feature name (no
    ``ID`` or ``Name``) attribute defined in the feature then ``None``
    is returned.

    :param feature: GFF feature
    :type feature: gffutils.feature.Feature
    :param use_feature_name: If a feature defines both ``ID`` and ``Name`` \
    attributes then use ``Name`` as its identifier, otherwise use ``ID``.
    :type use_feature_name: bool
    :return: Feature name or ``None``
    :rtype: str or unicode
    """
    name = None
    id_attr = None
    name_attr = None
    if "ID" in feature.attributes:
        id_attr = feature.attributes["ID"]
    if "Name" in feature.attributes:
        name_attr = feature.attributes["Name"]
    if id_attr and name_attr:
        if use_feature_name:
            name = name_attr
        else:
            name = id_attr
    elif id_attr:
        name = id_attr
    elif name_attr:
        name = name_attr
    if name:
        name = name[0].strip()
    return name


def get_cds_from_fasta(feature, fasta):
    """
    Get the sequence of a CDS from a FASTA file using a GFF feature
    for the CDS.

    :param feature: GFF feature for the CDS
    :type feature: gffutils.feature.Feature
    :param fasta: FASTA genes
    :type fasta: pyfaidx.Fasta
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
                              cds_feature_format=CDS_FEATURE_FORMAT,
                              use_feature_name=False):
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
    :param use_feature_name: If a feature defines both ``ID`` and \
    ``Name`` attributes then use ``Name`` in reporting, otherwise use \
    ``ID``.
    :type use_feature_name: bool
    :return: Codons for each coding sequence, keyed by feature name
    :rtype: dict(str or unicode -> list(str or unicode))
    :raises pyfaidx.FastaIndexingError: If the FASTA file has badly \
    formatted sequences
    :raises FileNotFoundError: If the FASTA or GFF files \
    cannot be found
    :raises ValueError: If GFF file is empty
    :raises Exception: Exceptions specific to gffutils.create_db \
    (these are undocumented in the gffutils documentation)
    """
    for f in [fasta, gff]:
        if not os.path.exists(f) or (not os.path.isfile(f)):
            raise FileNotFoundError(f)
    try:
        gffdb = gffutils.create_db(gff,
                                   dbfn='gff.db',
                                   force=True,
                                   keep_order=True,
                                   merge_strategy='merge',
                                   sort_attribute_values=True)
    except ValueError as e:
        # Wrap and rethrow exception so file name is included
        raise ValueError("{} ({})".format(e, gff)) from e
    cds_codons = {}
    same_feature_id_count = 0
    fasta_genes = Fasta(fasta)
    for feature in gffdb.features_of_type('CDS'):
        try:
            sequence = get_cds_from_fasta(feature, fasta_genes)
            codons = sequence_to_codons(sequence)
        except KeyError as e:  # Missing sequence.
            warnings.warn(str(e))
            continue
        except AssertionError as e:  # Sequence length not divisible by 3.
            warnings.warn(str(e))
            continue
        feature_id = get_feature_id(feature, use_feature_name)
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


def write_feature_codons_to_csv(feature_codons, csv_file, delimiter="\t"):
    """
    Write a dictionary of the codons for features, keyed by feature
    name into a CSV file, including a header.

    The CSV file has columns:

    * :py:const:`GENE`: feature name.
    * :py:const:`POS`: codon position in coding sequence (1-indexed).
    * :py:const:`CODON`: codon.

    :param feature_codons: Codons for each feature, keyed by feature \
    name
    :type feature_codons: dict(str or unicode -> list(str or unicode))
    :param csv_file: CSV file name
    :type csv_file: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    provenance.write_provenance_header(__file__, csv_file)
    with open(csv_file, "a") as f:
        writer = csv.writer(f, delimiter=delimiter, lineterminator='\n')
        writer.writerow([GENE, POS, CODON])
        for feature_id, codons in list(feature_codons.items()):
            for pos, codon in zip(range(0, len(codons)), codons):
                writer.writerow([feature_id, pos+1, codon])


def get_cds_codons_file(fasta,
                        gff,
                        cds_codons_file,
                        exclude_stop_codons=False,
                        cds_feature_format=CDS_FEATURE_FORMAT,
                        use_feature_name=False,
                        delimiter="\t"):
    """
    Using CDS entries within a GFF file, get the codons in each coding
    sequence in the complementary FASTA file.

    A tab-separated values file of the codons for each CDS, keyed by
    CDS feature name, is saved.

    See :py:func:`get_cds_codons_from_fasta`.

    See :py:func:`write_feature_codons_to_csv` for tab-separated values
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
    :param use_feature_name: If a feature defines both ``ID`` and \
    ``Name`` attributes then use ``Name`` in reporting, otherwise use \
    ``ID``.
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    :raises FileNotFoundError: If the FASTA or GFF files \
    cannot be found
    :raises pyfaidx.FastaIndexingError: If the FASTA file has badly \
    formatted sequences
    :raises ValueError: If GFF file is empty
    :raises Exception: Exceptions specific to gffutils.create_db \
    (these are undocumented in the gffutils documentation)
    """
    cds_codons = get_cds_codons_from_fasta(fasta,
                                           gff,
                                           exclude_stop_codons,
                                           cds_feature_format,
                                           use_feature_name)
    write_feature_codons_to_csv(cds_codons, cds_codons_file, delimiter)
