#! python
"""
Creates simple simulated fastq files to test UMI/deduplication,
adaptor trimming, and demultiplexing.

Usage:

    python -m riboviz.tools.create_fastq_examples DIRECTORY

TODO:

- Add barcodes for demultiplexing.
"""

import os
import os.path
from random import choices, seed
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


QUALITY_MEDIUM = list(range(30, 41))
""" List of medium quality scores. """
QUALITY_HIGH = list(range(39, 41))
""" List of high quality scores. """
FASTQ_FORMAT = "fastq"
""" Format string for use with Bio.SeqIO.write. """


def simulate_quality(k, qualities=QUALITY_MEDIUM, weights=None):
    """
    Simulate quality scores. This is a thin wrapper around
    random.choices whose default values represent medium Phred quality
    values.

    See https://docs.python.org/3/library/random.html#random.choices

    :param k: Number of quality scores requested
    :type k: int
    :param qualities: Available quality scores
    :type qualities: list(int)
    :param weights: Optional weightings for qualities
    :type weights: list(ints)
    :return k quality scores, selected at random from qualities
    :rtype: list(int)
    """
    return choices(qualities, k=k, weights=weights)


def make_fastq_record(name, reads, scores=None, qualities=QUALITY_MEDIUM):
    """
    Make a fastq record with sequence readall and name.

    :param name: Name
    :type name: str or unicode
    :param reads: Reads
    :type reads: str or unicode
    :param scores: Quality scores
    :type scores: list(int)
    :param qualities: Available quality scores, used if scores is
    None, in conjunction with simulate_quality to calculate quality
    scores
    :type qualities: list(int)
    :return: fastq record
    :rtype: Bio.SeqRecord.SeqRecord
    """
    if scores is None:
        scores = simulate_quality(len(reads), qualities=qualities)
    record = SeqRecord(Seq(reads, IUPAC.ambiguous_dna),
                       id=name,
                       name=name,
                       description=name)
    record.letter_annotations["phred_quality"] = scores
    return record


def trim_fastq_record_trim_3prime(record, trim, add_trim=False, delimiter="_"):
    """
    Copy fastq record, but trim sequence and quality scores at 3' end
    by given length.

    :param record: fastq record
    :type record: Bio.SeqRecord.SeqRecord
    :param trim: Number of nts to trim by
    :type trim: int
    :param add_trim: Add subsequence that was removed to record ID
    :type add_trim: bool
    :param delimiter: Delimiter to use, if add_trim is True
    :type delimiter: str or unicode
    :return: fastq record
    :rtype: Bio.SeqRecord.SeqRecord
    """
    quality = record.letter_annotations["phred_quality"]
    sequence = str(record.seq)
    record_extension = ""
    if add_trim:
        record_extension = delimiter + sequence[-trim:]
    return make_fastq_record(record.id + record_extension,
                             sequence[0:-trim],
                             quality[0:-trim])


def trim_fastq_record_trim_5prime(record, trim, add_trim=False, delimiter="_"):
    """
    Copy fastq record, but trim sequence and quality scores at 5' end
    by given length.

    :param record: fastq record
    :type record: Bio.SeqRecord.SeqRecord
    :param trim: Number of nts to trim by
    :type trim: int
    :param add_trim: Add subsequence that was removed to record ID
    :type add_trim: bool
    :param delimiter: Delimiter to use, if add_trim is True
    :type delimiter: str or unicode
    :return: fastq record
    :rtype: Bio.SeqRecord.SeqRecord
    """
    quality = record.letter_annotations["phred_quality"]
    sequence = str(record.seq)
    record_extension = ""
    if add_trim:
        record_extension = delimiter + sequence[0:trim]
    return make_fastq_record(record.id + record_extension,
                             sequence[trim:],
                             quality[trim:])


def create_simulated_fastq_files(output_dir):
    """
    Create simulated fastq files.

    :param output_dir: Output directory
    :type output_dir: str or unicode
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    make_fastq_record("SRR", "AAAA")
    seed(42)  # Fix random seed so can repeatedly create same files

    # Components for simulated reads compatible with the vignette
    # files yeast_yAL_CDS_w_250utrs.fa.
    # These are aimed at the Duncan & Mata format with 4nt UMI at
    # each end of read.
    read_a = "ATGGCATCCACCGATTTCTCCAAGATTGAA"  # 30nt starting ORF of YAL003W
    read_ae = "ATGGCATCCACCGATGTCTCCAAGATTGAA"  # 1 error in read A
    read_b = "TCTAGATTAGAAAGATTGACCTCATTAA"  # 28nt immediately following start of ORF of YAL038W

    # UMIs
    umi_x = "CGTA"
    umi_y = "ATAT"
    umi_ye = "ATAA"  # 1 error in umi_y
    umi_z = "CGGC"
    umi_ze = "CTGC"  # 1 error in umi_x
    umi_3_length = 4
    umi_5 = "AAAA"
    umi_5c = "CCCC"
    umi_5_length = 4

    adaptor = "CTGTAGGCACC"  # Adaptor sequence used in vignette data

    post_adaptor_nt = "AC"

    # Simulate raw data with 3' and 5' umi_s, i.e. desired input before
    # adaptor trimming.
    # M.N in the record names note the group the record is expected
    # belong to (M) and its number within that group.
    # After deduplication there should only be 1 member of each group.
    umi5_read_umi3_adaptor_records = [
        make_fastq_record("EWSim-1.1-umi5-reada-umix",
                          umi_5 + read_a + umi_x + adaptor,
                          qualities=QUALITY_HIGH),
        make_fastq_record("EWSim-1.2-umi5-reada-umix",
                          umi_5 + read_a + umi_x + adaptor),
        make_fastq_record("EWSim-1.3-umi5-readae-umix",
                          umi_5 + read_ae + umi_x + adaptor),
        make_fastq_record("EWSim-2.1-umi5-reada-umiy",
                          umi_5 + read_a + umi_y + adaptor)
    ]
    # Extra nt past the adaptor for the shorter read.
    umi5_read_umi3_adaptor_2nt_records = [
        make_fastq_record("EWSim-3.1-umi5-readb-umix",
                          umi_5 + read_b + umi_x + adaptor + post_adaptor_nt),
        make_fastq_record("EWSim-4.1-umi5-readb-umiz",
                          umi_5 + read_b + umi_z + adaptor + post_adaptor_nt,
                          qualities=QUALITY_HIGH),
        make_fastq_record("EWSim-4.2-umi5-readb-umiz",
                          umi_5 + read_b + umi_z + adaptor + post_adaptor_nt),
        make_fastq_record("EWSim-4.3-umi5-readb-umize",
                          umi_5 + read_b + umi_ze + adaptor + post_adaptor_nt),
        make_fastq_record("EWSim-5.1-umi5c-readb-umix",
                          umi_5c + read_b + umi_x + adaptor + post_adaptor_nt)
    ]
    with open(os.path.join(output_dir, "simdata_UMI5and3_4nt_adaptor.fastq"),
              "w") as f:
        SeqIO.write(umi5_read_umi3_adaptor_records, f, FASTQ_FORMAT)
        SeqIO.write(umi5_read_umi3_adaptor_2nt_records, f, FASTQ_FORMAT)

    # Simulate raw data with 3' and 5' umi_s, i.e. desired input after
    # adaptor trimming.
    # Reuse the records above but trim the sequences and quality
    # scores at the 3' end by the length of the adaptor.
    umi5_read_umi3_records = [
        trim_fastq_record_trim_3prime(record, len(adaptor))
        for record in umi5_read_umi3_adaptor_records
    ]
    umi5_read_umi3_2nt_records = [
        trim_fastq_record_trim_3prime(record,
                                      len(adaptor) +
                                      len(post_adaptor_nt))
        for record in umi5_read_umi3_adaptor_2nt_records
    ]
    with open(os.path.join(output_dir, "simdata_UMI5and3_4nt.fastq"),
              "w") as f:
        SeqIO.write(umi5_read_umi3_records, f, FASTQ_FORMAT)
        SeqIO.write(umi5_read_umi3_2nt_records, f, FASTQ_FORMAT)

    # Simulate raw data with only 3' umi_.
    read_umi3_records = [
        make_fastq_record("EWSim-1.1-reada-umix",
                          read_a + umi_x,
                          qualities=QUALITY_HIGH),
        make_fastq_record("EWSim-1.2-reada-umix", read_a + umi_x),
        make_fastq_record("EWSim-1.3-readae-umix", read_ae + umi_x),
        make_fastq_record("EWSim-2.1-reada-umiy", read_a + umi_y),
        make_fastq_record("EWSim-3.1-readb-umix", read_b + umi_x),
        make_fastq_record("EWSim-4.1-readb-umiz",
                          read_b + umi_z,
                          qualities=QUALITY_HIGH),
        make_fastq_record("EWSim-4.2-readb-umiz", read_b + umi_z),
        make_fastq_record("EWSim-4.3-readb-umize", read_b + umi_ze)
    ]
    with open(os.path.join(output_dir, "simdata_UMI3_4nt.fastq"),
              "w") as f:
        SeqIO.write(read_umi3_records, f, FASTQ_FORMAT)

    # Simulate raw data with no UMIs i.e. desired input after
    # adaptor trimming and UMI extraction.
    all_umi5_read_umi3_records = umi5_read_umi3_records \
        + umi5_read_umi3_2nt_records
    # Add UMI to record ID, using "_" delimiter for consistency with
    # UMI-tools.
    extracted_read_umi3_records = [
        trim_fastq_record_trim_5prime(record, umi_5_length, True)
        for record in all_umi5_read_umi3_records
    ]
    # Add UMI to record ID, using "" delimiter for consistency with
    # UMI-tools, which concatenates 5' and 3' UMIs.
    extracted_read_records = [
        trim_fastq_record_trim_3prime(record, umi_3_length, True, "")
        for record in extracted_read_umi3_records
    ]
    with open(os.path.join(output_dir, "simdata_extracted_UMI5and3_4nt.fastq"),
              "w") as f:
        SeqIO.write(extracted_read_records, f, FASTQ_FORMAT)

    extracted_read_records = [
        trim_fastq_record_trim_3prime(record, umi_3_length, True)
        for record in read_umi3_records
    ]
    with open(os.path.join(output_dir, "simdata_extracted_UMI3_4nt.fastq"),
              "w") as f:
        SeqIO.write(extracted_read_records, f, FASTQ_FORMAT)


if __name__ == "__main__":
    create_simulated_fastq_files(sys.argv[1])
