"""
Identify upstream open reading frames (uORFs).

Usage::

    python -m riboviz.tools.uorf_finde [-h] \
        -f FASTA -g GFF3 -o GFF3_OUTPUT -u FASTA_OUTPUT \
        -m MIN_LENGTH [-s START_CODONS [START_CODONS ...]] \
         [-e STOP_CODONS [STOP_CODONS ...]]

    -h, --help            show this help message and exit
    -f FASTA, --fasta FASTA
                          fasta file input
    -g GFF3, --gff3 GFF3  gff3 file input
    -o GFF3_OUTPUT, --gff3-output GFF3_OUTPUT
                          gff3 file output
    -u FASTA_OUTPUT, --fasta-output FASTA_OUTPUT
                          fasta file output with only UTR5 and CDS
                          regions
    -m MIN_LENGTH, --min-length MIN_LENGTH
                          Minimum number of codons for uORFs
    -s START_CODONS [START_CODONS ...], --start-codons START_CODONS
     [START_CODONS ...] 
                          Start codon (default ATG)
    -e STOP_CODONS [STOP_CODONS ...], --stop-codons STOP_CODONS
    -[STOP_CODONS ...]   
                          Stop codons (default TAA TAG TGA)
"""
import argparse
import os
import gffutils
import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',
                        '--fasta',
                        required=True,
                        help="fasta file input")
    parser.add_argument('-g',
                        '--gff3',
                        required=True,
                        help="gff3 file input")
    parser.add_argument('-o',
                        '--gff3-output',
                        required=True,
                        help="gff3 file output")
    parser.add_argument('-u',
                        '--fasta-output',
                        required=True,
                        help="fasta file output with only UTR5 and CDS regions")
    parser.add_argument('-m',
                        '--min-length',
                        type=int,
                        required=True,
                        help="Minimum number of codons for uORFs")
    parser.add_argument('-s',
                        '--start-codons',
                        default=['ATG'],
                        nargs="+",
                        help="Start codon (default ATG)")
    parser.add_argument('-e',
                        '--stop-codons',
                        default=['TAA', 'TAG', 'TGA'],
                        nargs="+",
                        help="Stop codons (default TAA TAG TGA)")
    options = parser.parse_args()
    return options


def extract_utr5_cds_postion(mygff3):
    db = gffutils.create_db(mygff3,
                            'myGFF.db',
                            merge_strategy="create_unique",
                            keep_order=True,
                            force=True)
    db = gffutils.FeatureDB('myGFF.db')
    position_dict = {}
    for i in db.features_of_type("UTR5"):
        seq_dict = {}
        for j in db.features_of_type("CDS"):
            if i.seqid == j.seqid:
                seq_dict[i.seqid] = [i.start-1, j.end, i.end]
                position_dict.update(seq_dict)
    return position_dict


def extract_fasta_sequence(myfasta, newfasta, position_dict):
    for key in position_dict.keys():
        for record in SeqIO.parse(myfasta, "fasta"):
            if key == record.id:
                sequences = SeqRecord(Seq(str(record.seq[position_dict[key][0]:position_dict[key][1]])), id=record.id)
                with open(newfasta, "a") as output_handle:
                    SeqIO.write(sequences, output_handle, "fasta")


def startstop_codon(newfasta, start_codons, stop_codons):
    new_fasta_iterator = SeqIO.parse(newfasta, "fasta")
    for record in new_fasta_iterator:
        seq = record.seq
        for frame in range(3):
            for i in range(frame, len(seq), 3):
                current_codon1 = seq[i:i+3]
                if current_codon1 in start_codons:
                    start_codon_index = i
                    for j in range(start_codon_index, len(seq), 3):
                        current_codon2 = seq[j:j+3]
                        if current_codon2 in stop_codons:
                            stop_codon_index = j
                            length = stop_codon_index - start_codon_index + 3
                            start_uorf = start_codon_index + 1
                            stop_uorf = stop_codon_index+3
                            yield (record.id, frame+1,
                                   start_uorf, stop_uorf, length)
                            break


def check_length(min_length, initial_list):
    val = [list(ele) for ele in initial_list]
    condition = lambda x: x[4] >= min_length * 3
    filter_length_list = list(filter(condition, val))
    return filter_length_list


def check_uorf_type(filter_length_list, position_dict):
    for key in position_dict.keys():
        for i in filter_length_list:
            if key == i[0]:
                end_utr5 = position_dict[key][2]
                end_cds = position_dict[key][1]
                start_uorf = i[2]
                stop_uorf = i[3]
                if start_uorf < end_utr5:
                    if stop_uorf < end_utr5:
                        i.append("uORF")
                    if stop_uorf == end_cds:
                        i.append("CDS_NTE")
                    if end_utr5 < stop_uorf < end_cds:
                        i.append("overlap_uORF")
    return filter_length_list


def remove_main_orfs(filter_length_type_list):
    condition = lambda x: len(x) == 6
    filter_length_type_orfs_list = list(filter(condition, filter_length_type_list))
    return filter_length_type_orfs_list


def longest_orf(filter_length_type_orfs_list):
    df = pd.DataFrame(filter_length_type_orfs_list)
    df.columns = ['name', 'frame', 'start_uorf', 'stop_uorf', 'length', 'type']
    df1 = df[df.duplicated('stop_uorf', keep=False)].groupby('stop_uorf')['start_uorf'].apply(list).reset_index()
    dic = df1.set_index('stop_uorf').T.to_dict('list')
    longest_pos_dict = {}
    for key, values in dic.items():
        longest_pos_dict[key] = min(dic[key][0])
    g = []
    for x in filter_length_type_orfs_list:
        for key, values in longest_pos_dict.items():
            if x[2] != values and x[3] == key:
                g.append(x)
    final_list = [item for item in filter_length_type_orfs_list if item not in g]
    return final_list


def write_to_gff3(final_list, output_file, newfasta):
    new_fasta_iterator = SeqIO.parse(newfasta, "fasta")
    for record in new_fasta_iterator:
        for i in final_list:
            if record.id == i[0]:
                seq = record.seq
                name = i[0]
                rec = SeqRecord(seq, name)
                qualifiers = {"source": "riboviz",
                              "score": ".",
                              "start_codon": seq[i[2]-1:i[2]+2],
                              "Name": name + "_" + i[5] + "_" + str(i[2]),
                              "frame":i[1]}
                feature = SeqFeature(FeatureLocation(i[2]-1, i[3]),
                                     type=i[5],
                                     strand=1,
                                     qualifiers=qualifiers)
                rec.features = [feature]
                with open(output_file, "a") as out_handle:
                    GFF.write([rec], out_handle)


def invoke_uorfs_finder():
    options = parse_command_line_options()
    fasta = options.fasta
    gff3 = options.gff3
    gff3_output = options.gff3_output
    min_length = options.min_length
    start_codons = options.start_codons
    stop_codons = options.stop_codons
    fasta_output = options.fasta_output
    position_dict = extract_utr5_cds_postion(gff3)
    extract_fasta_sequence(fasta, fasta_output, position_dict)
    initial_list = list(startstop_codon(fasta_output, start_codons, stop_codons))
    filter_length_list = check_length(min_length, initial_list)
    filter_length_type_list = check_uorf_type(filter_length_list,
                                              position_dict)
    filter_length_type_orfs_list = remove_main_orfs(filter_length_type_list)
    final_list = longest_orf(filter_length_type_orfs_list)
    write_to_gff3(final_list, gff3_output, fasta_output)
    # TODO recode. This removes header comments that are present
    # for every record in the GFF file.
    os.system("sed -i '/#/d' {}" .format(gff3_output))


if __name__ == "__main__":
    invoke_uorfs_finder()
