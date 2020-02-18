#! python
"""
Check fasta and gff files are compatible with RiboViz.

Usage:

    python -m riboviz.tools.check_fasta_gff [-h] \
        -f FASTAIN -g GFFIN

Example:

    python -m riboviz.tools.check_fasta_gff \
        -f vignette/input/yeast_YAL_CDS_w_250utrs.fa
        -g vignette/input/yeast_YAL_CDS_w_250utrs.gff3

The script checks that:

* The beginning of every CDS is a start codon (ATG; translates to
  'M').
* The stop of every CDS is a stop codon (TAG, TGA, TAA; translates to
  '*')
* There are no stop codons internal to the CDS.

Some unusual genes (e.g. frameshifts) might not have this.
"""
import argparse
import warnings
from Bio.Seq import Seq
import gffutils
from riboviz import provenance


def check_fasta_gff_files(fasta, gff):
    """
    Check gff and fasta files for compatibility.

    :param fasta: fasta file
    :type fasta: str or unicode
    :param gff: gff file
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


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Check fasta and gff files have start and stop codons as expected")
    parser.add_argument("-f",
                        "--fasta",
                        dest="fasta",
                        required=True,
                        help="fasta file input")
    parser.add_argument("-g",
                        "--gff",
                        dest="gff",
                        required=True,
                        help="gff3 file input")
    options = parser.parse_args()
    return options


def invoke_check_fasta_gff():
    """
    Parse command-line options then invoke "check_fasta_gff_files".
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    fasta = options.fasta
    gff = options.gff
    check_fasta_gff_files(fasta, gff)


if __name__ == "__main__":
    invoke_check_fasta_gff()
