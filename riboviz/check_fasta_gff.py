"""
FASTA and GFF compatibility functions.
"""
import warnings
from Bio.Seq import Seq
import gffutils


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
        # for all CDS entries in gff,
        # print(cds_coord.seqid)
        # extract CDS
        try:
            cds_seq = cds_coord.sequence(fasta)
        except Exception as e:
            print(str(e))
            continue
        cds_len_remainder = len(cds_seq) % 3
        if cds_len_remainder != 0:
            warnings.warn(
                cds_coord.seqid + " has length that isn't divisible by 3")
            cds_seq += ("N" * (3 - cds_len_remainder))

        cds_trans = Seq(cds_seq).translate()

        if cds_trans[0] != "M":
            print((cds_coord.seqid + " doesn't start with ATG."))
        if cds_trans[-1] != "*":
            print((cds_coord.seqid + " doesn't stop at end."))
        if any([L == "*" for L in cds_trans[:-1]]):
            print((cds_coord.seqid + " has internal STOP."))
