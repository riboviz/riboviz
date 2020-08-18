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
            # that can be thrown by get_cds_codons, from
            # gffutils.Feature.sequence is KeyError. This can arise if
            # the GFF file contains information on a sequence that is
            # not in the FASTA file.
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


def extract_cds(gff):
    """
    Extract information on genes and coding sequences from a GFF
    file. A dictionary indexed by gene ID is returned. For each
    gene ID is a list of tuples of form (CDS start position, CDS end
    position).

    :param gff: GFF file
    :type gff: str or unicode
    :return: coding sequence information
    :rtype: dict(str or unicode -> list(tuple(int, int)))
    """
    db = gffutils.create_db(gff,
                            'gff.db',
                            merge_strategy="create_unique",
                            keep_order=True,
                            force=True)
    db = gffutils.FeatureDB('gff.db')
    cds = {}
    for c in db.features_of_type("CDS"):
        if c.seqid not in cds:
            cds[c.seqid] = []
        cds[c.seqid].append((c.start, c.end))
    return cds


def get_cds(cds_coord, fasta):
    """
    Given a GFF feature for a CDS extract the coding
    sequence from the given FASTA file.
    If the length of coding sequence is not a factor of
    3 then it is padded with ``N``.

    :param cds_coord: GFF feature for a CDS
    :type cds_coord: gffutils.feature.Feature
    :param fasta: FASTA file
    :type fasta: str or unicode
    :return: sequence
    :rtype: str or unicode
    :raises Exception: Exceptions specific to
    gffutils.Feature.sequence (these are undocumented in the gffutils
    documentation)
    """
    sequence = cds_coord.sequence(fasta)
    cds_len_remainder = len(sequence) % 3
    if cds_len_remainder != 0:
        warnings.warn("{} has length not divisible by 3".format(
            cds_coord.seqid))
        sequence += ("N" * (3 - cds_len_remainder))
    return sequence


def get_cds_codons(cds_coord, fasta):
    """
    Given a GFF feature for a CDS extract the coding
    sequence from the given FASTA file and split into codons.

    See :py:func:`get_cds` and :py:func:`sequence_to_codons`.

    :param cds_coord: GFF feature for a CDS
    :type cds_coord: gffutils.feature.Feature
    :param fasta: FASTA file
    :type fasta: str or unicode
    :return: sequence
    :rtype: str or unicode
    """
    sequence = get_cds(cds_coord, fasta)
    cds_codons = sequence_to_codons(sequence)
    return cds_codons


def get_seqs_cds_codons(fasta, gff):
    """
    Using information within a FASTA file extract information on the
    codons in each coding sequence (as specified via CDS entries in
    the complementary GFF file). A dictionary of the codons for each
    coding sequence, keyed by sequence ID, is returned.

    See :py:func:`get_cds_codons`.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :return: Codons for each coding sequence, keyed by sequence ID
    :rtype: list(dict)
    :raises Exception: Exceptions specific to gffutils.create_db
    (these are undocumented in the gffutils documentation)
    """
    gffdb = gffutils.create_db(gff,
                               dbfn='gff.db',
                               force=True,
                               keep_order=True,
                               merge_strategy='merge',
                               sort_attribute_values=True)
    seqs_cds_codons = {}
    for cds_coord in gffdb.features_of_type('CDS'):
        try:
            cds_codons = get_cds_codons(cds_coord, fasta)
        except Exception as e:
            # Log and continue with other CDSs. A typical exception
            # that can be thrown by get_cds_codons, from
            # gffutils.Feature.sequence is KeyError. This can arise if
            # the GFF file contains information on a sequence that is
            # not in the FASTA file.
            warnings.warn(str(e))
            continue
        if cds_coord.seqid not in seqs_cds_codons:
            seqs_cds_codons[cds_coord.seqid] = []
        seqs_cds_codons[cds_coord.seqid].extend(cds_codons)
    return seqs_cds_codons


def seqs_cds_codons_to_df(seqs_cds_codons):
    """
    Given dictionary of the codons for coding sequences, keyed by
    sequence ID, return a Pandas data frame with the codons, also
    including the position of each codon in its sequence.

    The data frame has columns:

    * :py:const:`GENE`: sequence ID/gene name.
    * :py:const:`CODON`: codon.
    * :py:const:`POS`: codon position in coding sequence (1-indexed).

    :param cds_codons: Codons for each coding sequence, keyed by
    sequence ID
    :type cds_codons: list(dict)
    :return: data frame
    :rtype: pandas.core.frame.DataFrame
    """
    seqs_cds_codons_list = []
    num_seqs_cds_codons = 0
    for seqid, cds_codons in list(seqs_cds_codons.items()):
        num_seqs_cds_codons += len(cds_codons)
        for pos, cds_codon in zip(range(0, len(cds_codons)), cds_codons):
            seqs_cds_codons_list.append({GENE: seqid,
                                         CODON: cds_codon,
                                         POS: pos + 1})
    # Create empty DataFrame so if seqs_cds_codons and
    # seqs_cds_codons_list are empty we still have an empty DataFrame
    # with the column names.
    df = pd.DataFrame(columns=[GENE, CODON, POS])
    df = df.append(pd.DataFrame(seqs_cds_codons_list),
                   ignore_index=True)
    # Validate number of codons
    assert num_seqs_cds_codons == df.shape[0]
    return df


def extract_cds_codons(fasta,
                       gff,
                       cds_codons_file,
                       delimiter="\t"):
    """
    Using information within a FASTA file extract information on the
    codons in each coding sequence (as specified via CDS entries in
    the complementary GFF file) and save a tab-separated values file
    with information on the positions of each codon in each coding
    sequence.

    The tab-separated values file has columns:

    * :py:const:`GENE`: sequence ID/gene name.
    * :py:const:`CODON`: codon.
    * :py:const:`POS`: codon position in coding sequence (1-indexed).

    See :py:func:`get_seqs_cds_codons` and
    :py:func:`seqs_cds_codons_to_df`.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :params cds_codons_file: Coding sequence codons file
    :type cds_codons_file: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    seqs_cds_codons = get_seqs_cds_codons(fasta, gff)
    df = seqs_cds_codons_to_df(seqs_cds_codons)
    provenance.write_provenance_header(__file__, cds_codons_file)
    df[list(df.columns)].to_csv(cds_codons_file,
                                mode='a',
                                sep=delimiter,
                                index=False)
