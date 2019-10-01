#! python
"""
simulate_fastq_tests.py

Creates simple simulated fastq files to test UMI/deduplication,
adaptor trimming, and demultiplexing.

To-do:

- Matched functions to make fastq records from a read and quality
  score plus UMI at 3' end:
  - One places the UMI in-line in the read (raw data, i.e. input for
    extract; this is done already)
  - Function to place the UMI in the header (i.e. desired output for
    extract)
  - These should have the same quality score so we can check quality
    is preserved by extract.
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
    # umi_tools keeps the highest-quality read, so we use QUALITY_HIGH
    # for reads we wish to keep.
    umi_5and3_4nt_adaptor_records = [
        make_fastq_record("EWSim1Umi5read_aUmiX.Keep",
                          umi_5 + read_a + umi_x + adaptor,
                          qualities=QUALITY_HIGH),
        make_fastq_record("EWSim2Umi5read_aUmiX.Drop",
                          umi_5 + read_a + umi_x + adaptor),
        make_fastq_record("EWSim3Umi5read_aeUmiX.Drop",
                          umi_5 + read_ae + umi_x + adaptor),
        make_fastq_record("EWSim4Umi5read_aUmiY.Keep",
                          umi_5 + read_a + umi_y + adaptor)
    ]
    # Extra nt past the adaptor for the shorter read.
    umi_5and3_4nt_adaptor_extra_nt_records = [
        make_fastq_record("EWSim5Umi5read_bUmiX.Keep",
                          umi_5 + read_b + umi_x + adaptor + post_adaptor_nt),
        make_fastq_record("EWSim6Umi5read_bUmiZ.Keep",
                          umi_5 + read_b + umi_z + adaptor + post_adaptor_nt,
                          qualities=QUALITY_HIGH),
        make_fastq_record("EWSim7Umi5read_bUmiZ.Drop",
                          umi_5 + read_b + umi_z + adaptor + post_adaptor_nt),
        make_fastq_record("EWSim8Umi5read_bUmiZe.Drop",
                          umi_5 + read_b + umi_ze + adaptor + post_adaptor_nt),
        make_fastq_record("EWSim9Umi5Cread_bUmiX.Keep",
                          umi_5c + read_b + umi_x + adaptor + post_adaptor_nt)
    ]
    with open(os.path.join(output_dir, "simdata_UMI5and3_4nt_adaptor.fastq"),
              "w") as f:
        SeqIO.write(umi_5and3_4nt_adaptor_records, f, FASTQ_FORMAT)
        SeqIO.write(umi_5and3_4nt_adaptor_extra_nt_records, f, FASTQ_FORMAT)

    # Simulate raw data with 3' and 5' umi_s, i.e. desired input after
    # adaptor trimming.
    # Reuse the records above but trim the sequences and quality
    # scores at the 3' end by the length of the adaptor.
    umi_5and3_4nt_records = [
        trim_fastq_record_trim_3prime(record, len(adaptor))
        for record in umi_5and3_4nt_adaptor_records
    ]
    umi_5and3_4nt_extra_nt_records = [
        trim_fastq_record_trim_3prime(record,
                                      len(adaptor) +
                                      len(post_adaptor_nt))
        for record in umi_5and3_4nt_adaptor_extra_nt_records
    ]
    with open(os.path.join(output_dir, "simdata_UMI5and3_4nt.fastq"),
              "w") as f:
        SeqIO.write(umi_5and3_4nt_records, f, FASTQ_FORMAT)
        SeqIO.write(umi_5and3_4nt_extra_nt_records, f, FASTQ_FORMAT)

    # Simulate raw data with only 3' umi_.
    umi_3_4nt_records = [
        make_fastq_record("EWSim1ReadAUmiX.Keep",
                          read_a + umi_x,
                          qualities=QUALITY_HIGH),
        make_fastq_record("EWSim2ReadAUmiX.Drop", read_a + umi_x),
        make_fastq_record("EWSim3ReadAeUmiX.Drop", read_ae + umi_x),
        make_fastq_record("EWSim4ReadAUmiY.Keep", read_a + umi_y),
        make_fastq_record("EWSim5ReadBUmiX.Keep", read_b + umi_x),
        make_fastq_record("EWSim6ReadBUmiZ.Keep",
                          read_b + umi_z,
                          qualities=QUALITY_HIGH),
        make_fastq_record("EWSim7ReadBUmiZ.Drop", read_b + umi_z),
        make_fastq_record("EWSim8ReadBUmiZe.Drop", read_b + umi_ze)
    ]
    with open(os.path.join(output_dir, "simdata_UMI3_4nt.fastq"),
              "w") as f:
        SeqIO.write(umi_3_4nt_records, f, FASTQ_FORMAT)

    # Simulate raw data with no UMIs i.e. desired input after
    # adaptor trimming and UMI extraction.
    all_umi_5and3_4nt_records = umi_5and3_4nt_records \
        + umi_5and3_4nt_extra_nt_records
    # Add UMI to record ID, using "_" delimiter for consistency with
    # UMI-tools.
    extracted_umi_5and3_4nt_records = [
        trim_fastq_record_trim_5prime(record, umi_5_length, True)
        for record in all_umi_5and3_4nt_records
    ]
    # Add UMI to record ID, using "" delimiter for consistency with
    # UMI-tools, which concatenates 5' and 3' UMIs.
    extracted_umi_5and3_4nt_records = [
        trim_fastq_record_trim_3prime(record, umi_3_length, True, "")
        for record in extracted_umi_5and3_4nt_records
    ]
    with open(os.path.join(output_dir, "simdata_extracted_UMI5and3_4nt.fastq"),
              "w") as f:
        SeqIO.write(extracted_umi_5and3_4nt_records, f, FASTQ_FORMAT)

    extracted_umi_3_4nt_records = [
        trim_fastq_record_trim_3prime(record, umi_3_length, True)
        for record in umi_3_4nt_records
    ]
    with open(os.path.join(output_dir, "simdata_extracted_UMI3_4nt.fastq"),
              "w") as f:
        SeqIO.write(extracted_umi_3_4nt_records, f, FASTQ_FORMAT)


if __name__ == "__main__":
    create_simulated_fastq_files(sys.argv[1])
