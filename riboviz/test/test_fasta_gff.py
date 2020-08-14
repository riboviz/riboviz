"""
:py:mod:`riboviz.test_fasta_gff` tests.
"""
import pytest
from riboviz import fasta_gff


@pytest.mark.parametrize("seq_codons", [
    ("", []),
    ("A", ["A"]),
    ("GA", ["GA"]),
    ("GAT", ["GAT"]),
    ("GATT", ["GAT", "T"]),
    ("GATTACCA", ["GAT", "TAC", "CA"])])
def test_sequence_to_codons(seq_codons):
    """
    Test :py:func:`riboviz.fasta_gff.sequence_to_codons`
    with values vaid sequences.

    :param seq_codons: sequence and expected codons
    :type seq_codons: tuple(str or unicode, list(str or unicode))
    """
    sequence, codons = seq_codons
    assert fasta_gff.sequence_to_codons(sequence) == codons,\
        "Unexpected codons"


def test_seqs_cds_codons_to_df_empty():
    """
    Test :py:func:`riboviz.fasta_gff.seqs_cds_codons_to_df`
    with no values produces an empty data frame.
    """
    df = fasta_gff.seqs_cds_codons_to_df({})
    assert df.shape == (0, 0), "Unexpected number of rows and columns"


def test_seqs_cds_codons_to_df():
    """
    Test :py:func:`riboviz.fasta_gff.seqs_cds_codons_to_df`
    with values produces a data frame with the expected number
    of values.
    """
    seqs_cds_codons = {
        "G1": fasta_gff.sequence_to_codons("ATGAAATAA"),
        "G2": fasta_gff.sequence_to_codons("ATGGGGCCCTAG")
    }

    df = fasta_gff.seqs_cds_codons_to_df(seqs_cds_codons)

    num_rows, num_columns = df.shape
    assert num_columns == 3, "Unexpected number of columns"
    for column in [fasta_gff.GENE, fasta_gff.POS, fasta_gff.CODON]:
        assert column in list(df.columns),\
            "Missing column {}".format(column)
    assert num_rows == len(seqs_cds_codons["G1"]) \
        + len(seqs_cds_codons["G2"]), \
        "Unexpected number of rows"
    for seqid, cds_codons in seqs_cds_codons.items():
        seq_df = df.loc[df[fasta_gff.GENE] == seqid]
        assert not seq_df.empty, "Missing row for {}".format(seqid)
        num_rows, _ = seq_df.shape
        assert num_rows == len(cds_codons), \
            "Unexpected number of rows for {}".format(seqid)
        for pos, cds_codon in zip(seq_df[fasta_gff.POS],
                                  seq_df[fasta_gff.CODON]):
            # POS is 1-indexed.
            assert cds_codons[pos - 1] == cds_codon
