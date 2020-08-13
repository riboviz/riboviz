"""
FASTA and GFF functions.
"""
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
import pandas as pd

GENE = "Gene"
""" Column name (gene name). """
POS = "Pos"
"""
Column name (codon position in coding sequence, 1-indexed by codon).
"""
CODON = "Codon"
""" Column name (codon). """


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
    db = gffutils.create_db(gff, 'gff.db',
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


def print_fasta_cds(cds, fasta):
    """
    Print coding sequences of sequences in a FASTA file.

    :param cds: coding sequence information
    :type cds dict(str or unicode -> list(tuple(int, int)))
    :param fasta: FASTA file
    :type fasta: str or unicode
    """
    with open(fasta, "rt") as f:
        for seq in SeqIO.parse(f, "fasta"):
            # seq methods:
            # 'annotations', 'dbxrefs', 'description', 'features',
            # 'format', 'id', 'letter_annotations', 'lower', 'name',
            # 'reverse_complement', 'seq', 'translate', 'upper'
            if seq.name in cds:
                for (cds_start, cds_end) in cds[seq.name]:
                    print("---- {} ----".format(seq.name))
                    seq_pre_cds = seq.seq[0:cds_start-1]
                    seq_cds = seq.seq[cds_start-1:cds_end]
                    seq_post_cds = seq.seq[cds_end:]
                    print(len(seq_pre_cds))
                    print(seq_pre_cds)
                    print(len(seq_cds))
                    print(seq_cds)
                    print(len(seq_post_cds))
                    print(seq_post_cds)
                    print(len(seq.seq))
                    # Validate sequence was split correctly.
                    assert seq.seq == seq_pre_cds + seq_cds + seq_post_cds


# TODO function comment.
def get_cds_sequence(cds_coord, fasta):
    cds_len = cds_coord.end - cds_coord.start + 1
    print("{} {} {} ({})".format(cds_coord.seqid,
                                 cds_coord.start,
                                 cds_coord.end,
                                 cds_len))
    # cds_coord methods:
    # 'astuple', 'attributes', 'bin', 'calc_bin', 'chrom',
    # 'dialect', 'end', 'extra', 'featuretype', 'file_order',
    # 'frame', 'id', 'keep_order', 'score', 'seqid', 'sequence',
    # 'sort_attribute_values', 'source', 'start', 'stop', 'strand'
    # TODO what exceptions can this throw? Want to catch in caller.
    cds_seq = cds_coord.sequence(fasta)
    # Validate length was calculated correctly.
    assert cds_len == len(cds_seq)  # TODO might this ever not hold?
    cds_len_remainder = cds_len % 3
    # TODO from check_fasta_gff, should this be done?
    if cds_len_remainder != 0:
        warnings.warn("{} has length that isn't divisible by 3".format(
            cds_coord.seqid))
        cds_seq += ("N" * (3 - cds_len_remainder))
    return cds_seq


# TODO function comment.
def split_into_codons(cds_seq):
    codons = [cds_seq[i:i+3] for i in range(0, len(cds_seq), 3)]
    # Validate split was done correctly.
    assert "".join(codons) == cds_seq
    return codons


# TODO function comment.
def codons_list_to_dict(seqid, codons):
    codons_list = []
    for pos, codon in zip(range(0, len(codons)), codons):
        codons_list.append({GENE: seqid,
                            CODON: codon,
                            POS: pos + 1})
    return codons_list

                    
def extract_cds_codons(fasta, gff):
    """
    Using information within a FASTA file extract information on the
    codons in each coding sequence as specified via CDS entries in the
    complementary GFF file. A with columns  :py:const:`GENE` (gene
    name), :py:const:`CODON` (codon),  :py:const:`POS` (codon position
    in coding sequenced, 1-indexed) is returned.

    :param fasta: FASTA file
    :type fasta: str or unicode
    :param gff: GFF file
    :type gff: str or unicode
    :return: S
    :rtype: pandas.core.frame.DataFrame
    """
    gffdb = gffutils.create_db(gff,
                               dbfn='gff.db',
                               force=True,
                               keep_order=True,
                               merge_strategy='merge',
                               sort_attribute_values=True)
    # TODO Extract into new function, return DataFrame.
    all_cds_codons = []
    num_codons = 0
    for cds_coord in gffdb.features_of_type('CDS'):
        try:
            cds_seq = get_cds_sequence(cds_coord, fasta)
        except Exception as e:
            # TODO exit here rather than gulp and continue?
            print(str(e))
            continue
        print(cds_seq)
        codons = split_into_codons(cds_seq)
        all_cds_codons.extend(codons_list_to_dict(cds_coord.seqid, codons))
        num_codons += len(codons)
    df = pd.DataFrame(columns=[GENE, CODON, POS])
    print(df)
    df = pd.DataFrame(all_cds_codons)
    print(df)
    # Validate number of codons
    assert num_codons == df.shape[0]
    return df


if __name__ == "__main__":
    import sys
    cds = extract_cds(sys.argv[2])
    print("---- CDS ----")
    print(cds)
    print_fasta_cds(cds, sys.argv[1])
    cds = extract_cds_codons(sys.argv[1], sys.argv[2])
    file_name = sys.argv[3]
    delimiter = "\t"
    cds[list(cds.columns)].to_csv(
        file_name, mode='a', sep=delimiter, index=False)
