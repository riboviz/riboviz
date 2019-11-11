#! python
"""
Creates simple simulated FASTQ files to test UMI/deduplication,
adaptor trimming, and demultiplexing.

Usage:

    python -m riboviz.tools.create_fastq_examples DIRECTORY

where `DIRECTORY` is the directory into	which the simulated files are
to be written. The following files are created:

* `simdata_UMI5and3_4nt_adaptor.fastq`: FASTQ file with 9 reads,
  each with a 4nt UMI at the 5' end, a 4nt UMI at the 3' end and a
  11nt adaptor at the 3' end. Reads can be grouped by UMI into 5
  groups.
* `simdata_UMI5and3_4nt.fastq`: FASTQ file identical to the above but
  with the adaptor trimmed.
* `simdata_extracted_UMI5and3_4nt.fastq`: FASTQ file identical to the
  above but with the UMIs extracted and concatenated to the header.
* `simdata_extracted_UMI3_4nt.fastq`: FASTQ file with 8 reads, each
  with a 4nt UMI at the 3' end. Reads can be grouped by UMI into 4
  groups.
* `simdata_UMI3_4nt.fastq`: FASTQ file identical to the above but with
   the UMI extracted and concatenated to the header.

TODO:

- Add barcodes for demultiplexing.
"""

import csv
import os
import os.path
from random import choices, seed
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from riboviz.utils import BARCODE_DELIMITER
from riboviz.utils import UMI_DELIMITER


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


def trim_fastq_record_3prime(record,
                             trim,
                             add_trim=False,
                             delimiter=UMI_DELIMITER):
    """
    Copy fastq record, but trim sequence and quality scores at 3' end
    by given length.

    :param record: fastq record
    :type record: Bio.SeqRecord.SeqRecord
    :param trim: Number of nts to trim by
    :type trim: int
    :param delimiter: Delimiter to use, if add_trim is True
    :type delimiter: str or unicode
    :param add_trim: Add subsequence that was removed to record ID
    :type add_trim: bool
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


def trim_fastq_record_5prime(record,
                             trim,
                             add_trim=False,
                             delimiter=UMI_DELIMITER):
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


def make_fastq_records(tag,
                       read,
                       qualities,
                       umi5="",
                       umi3="",
                       barcode="",
                       adaptor="",
                       post_adaptor_nt=""):
    """
    Create a set of complementary fastq records.

    - A record for sequence: umi5 + read + umi3 + barcode + adaptor +
      post_adaptor_nt.
    - A record as above with the adaptor and post_adaptor_nt trimmed:
      umi5 + read + umi3 + barcode.
      - If adaptor and post_adaptor_nt are both "" then this is
        equivalent to the above record.
    - A record as above with the barcode and UMIs trimmed and added to
      the header with "_" delimiters: read.
      - If both the barcode and both UMIs are "" then this is
        equivalent to the above record.
      - Depending on whether barcode, umi5 and umi3 are "" the header
        will be extended with one of:
        - "<barcode>_<umi5><umi3>"
        - "<barcode>_<umi3>"
        - "_<umi5><umi3>"
        - "<umi3>"

    :param tag: Human-readable tag.
    :type tag: str or unicode
    :param read: Read
    :type read: str or unicode
    :param qualities: Available quality scores
    :type qualities: list(int)
    :param umi5: 5' end UMI
    :type umi5: str or unicode
    :param umi3: 3' end UMI
    :type umi3: str or unicode
    :param barcode: 3' end barcode
    :type barcode: str or unicode
    :param adaptor: 3' end adaptor
    :type adaptor: str or unicode
    :param post_adaptor_nt: 3' end post-adaptor nts
    :type post_adaptor_nt: str or unicode
    :returnL full record, adaptor-trimmed record, barcode- and
    UMI-extracted record
    :rtype: tuple(Bio.SeqRecord.SeqRecord, Bio.SeqRecord.SeqRecord,
    Bio.SeqRecord.SeqRecord)
    """
    sequence = umi5 + read + umi3 + barcode + adaptor + post_adaptor_nt
    record = make_fastq_record(tag, sequence, qualities=qualities)
    # Record after adaptor trimming.
    trim_record = trim_fastq_record_3prime(
        record,
        len(adaptor) + len(post_adaptor_nt))
    if barcode != "":
        # Add barcode to record ID, using "_" delimiter for consistency
        # with UMI-tools.
        barcode_ext_record = trim_fastq_record_3prime(
            trim_record,
            len(barcode),
            True,
            BARCODE_DELIMITER)
    else:
        barcode_ext_record = trim_record
    umi3_delimiter = UMI_DELIMITER
    if umi5 != "":
        # Record after 5' UMI extraction.
        # Add UMI to record ID, using "_" delimiter for consistency
        # with UMI-tools.
        umi5_ext_record = trim_fastq_record_5prime(
            barcode_ext_record,
            len(umi5),
            True,
            UMI_DELIMITER)
        umi3_delimiter = ""
    else:
        umi5_ext_record = barcode_ext_record
    # Record after 3' UMI extraction.
    # Add UMI to record ID, using "_" delimiter for consistency
    # with UMI-tools, unless 5' UMI has been extracted, in which case
    # use "".
    umi3_ext_record = trim_fastq_record_3prime(
        umi5_ext_record,
        len(umi3),
        True,
        umi3_delimiter)
    return (record, trim_record, umi3_ext_record)


def create_fastq_examples(output_dir):
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
    umi_5 = "AAAA"
    umi_5c = "CCCC"

    adaptor = "CTGTAGGCACC"  # Adaptor sequence used in vignette data

    post_adaptor_nt = "AC"

    # Create raw data with 5' and 3' UMIs and an adaptor.
    # M.N in the record names note the group the record is expected
    # belong to (M) and its number within that group.
    # After deduplication there should only be 1 member of each
    # group.
    config_5_3_adaptor = [
        ["EWSim-1.1-umi5-reada-umix", umi_5, read_a, umi_x, QUALITY_HIGH],
        ["EWSim-1.2-umi5-reada-umix", umi_5, read_a, umi_x, QUALITY_MEDIUM],
        ["EWSim-1.3-umi5-readae-umix", umi_5, read_ae, umi_x, QUALITY_MEDIUM],
        ["EWSim-2.1-umi5-reada-umiy", umi_5, read_a, umi_y, QUALITY_MEDIUM]
    ]
    # Create raw data with 5' and 3' UMIs and an adaptor plus an
    # extra nt past the "", adaptor for the shorter read.
    config_5_3_post_adaptor_nt = [
        ["EWSim-3.1-umi5-readb-umix", umi_5, read_b, umi_x, QUALITY_MEDIUM],
        ["EWSim-4.1-umi5-readb-umiz", umi_5, read_b, umi_z, QUALITY_HIGH],
        ["EWSim-4.2-umi5-readb-umiz", umi_5, read_b, umi_z, QUALITY_MEDIUM],
        ["EWSim-4.3-umi5-readb-umize", umi_5, read_b, umi_ze, QUALITY_MEDIUM],
        ["EWSim-5.1-umi5c-readb-umix", umi_5c, read_b, umi_x, QUALITY_MEDIUM]
    ]
    records = [
        make_fastq_records(tag, read, qualities, umi5, umi3, "", adaptor)
        for [tag, umi5, read, umi3, qualities] in config_5_3_adaptor]
    records_post_adaptor_nt = [
        make_fastq_records(tag, read, qualities, umi5, umi3, "",
                           adaptor, post_adaptor_nt)
        for [tag, umi5, read, umi3, qualities] in config_5_3_post_adaptor_nt]
    records.extend(records_post_adaptor_nt)
    file_names = ["simdata_UMI5and3_4nt_adaptor.fastq",
                  "simdata_UMI5and3_4nt.fastq",
                  "simdata_extracted_UMI5and3_4nt.fastq"]
    for file_name, fastq_records in zip(file_names, zip(*records)):
        with open(os.path.join(output_dir, file_name), "w") as f:
            SeqIO.write(fastq_records, f, FASTQ_FORMAT)

    # Simulate raw data with only 3' umi.
    config_3 = [
        ["EWSim-1.1-reada-umix", read_a, umi_x, QUALITY_HIGH],
        ["EWSim-1.2-reada-umix", read_a, umi_x, QUALITY_MEDIUM],
        ["EWSim-1.3-readae-umix", read_ae, umi_x, QUALITY_MEDIUM],
        ["EWSim-2.1-reada-umiy", read_a, umi_y, QUALITY_MEDIUM],
        ["EWSim-3.1-readb-umix", read_b, umi_x, QUALITY_MEDIUM],
        ["EWSim-4.1-readb-umiz", read_b, umi_z, QUALITY_HIGH],
        ["EWSim-4.2-readb-umiz", read_b, umi_z, QUALITY_MEDIUM],
        ["EWSim-4.3-readb-umize", read_b, umi_ze, QUALITY_MEDIUM],
    ]
    records = [
        make_fastq_records(tag, read, qualities, "", umi3, "", adaptor)
        for [tag, read, umi3, qualities] in config_3]
    file_names = ["simdata_UMI3_4nt_adaptor.fastq",
                  "simdata_UMI3_4nt.fastq",
                  "simdata_extracted_UMI3_4nt.fastq"]
    for file_name, fastq_records in zip(file_names, zip(*records)):
        with open(os.path.join(output_dir, file_name), "w") as f:
            SeqIO.write(fastq_records, f, FASTQ_FORMAT)

    # Create multiplexed data.
    # Use same data as 5' and 3' UMIs and an adaptor but with
    # barcodes.
    # Barcodes (keys) each with list of barcodes with 1-nt and 2-nt
    # mismatches.
    barcodes = {'ACG': ['ACT', 'TAG'],
                'GAC': ['GTC', 'TGA'],
                'CGA': ['TGA', 'CTT']}
    barcode_names = list(barcodes.keys())
    # Barcode that will be unassigned during demultiplexing.
    barcodes['TTT'] = []

    barcode_format = "-bar{:01d}.{:01d}"

    for file_name in ["simdata_multiplex.fastq",
                      "simdata_multiplex_trim.fastq",
                      "simdata_multiplex_trim_ext.fastq",
                      "simdata_multiplex_barcodes.tsv"]:
        file_path = os.path.join(output_dir, file_name)
        if os.path.exists(file_path):
            os.remove(file_path)
    with open(os.path.join(output_dir,
                           "simdata_multiplex_barcodes.tsv"), "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["SampleID", "TagRead"])
        for index, barcode in enumerate(barcode_names):
            writer.writerow(["Tag{:01d}".format(index), barcode])
    for barcode_index, barcode in enumerate(barcodes):
        matches = [barcode] + barcodes[barcode]
        for index, match in enumerate(matches):
            records = [
                make_fastq_records(tag +
                                   barcode_format.format(barcode_index, index),
                                   read, qualities,
                                   umi5, umi3, match,
                                   adaptor, "")
                for [tag, umi5, read, umi3, qualities] in config_5_3_adaptor]
            records_post_adaptor_nt = [
                make_fastq_records(tag +
                                   barcode_format.format(barcode_index, index),
                                   read, qualities,
                                   umi5, umi3, match,
                                   adaptor, post_adaptor_nt)
                for [tag, umi5, read, umi3, qualities] in config_5_3_post_adaptor_nt]
            records.extend(records_post_adaptor_nt)
            file_names = ["simdata_multiplex.fastq",
                          "simdata_multiplex_trim.fastq",
                          "simdata_multiplex_trim_ext.fastq"]
            for file_name, fastq_records in zip(file_names, zip(*records)):
                with open(os.path.join(output_dir, file_name), "a") as f:
                    SeqIO.write(fastq_records, f, FASTQ_FORMAT)

    # TODO GZIP


if __name__ == "__main__":
    create_fastq_examples(sys.argv[1])
