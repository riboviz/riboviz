#! python
"""
Demultiplex fastq files using UMI-tools-compliant barcodes present
within the fastq headers. These headers are assumed to be of form:

    @..._<BARCODE>_...

where the barcode is the first section delimited by underscores. If
another delimiter was used then that can be specified.

Inputs:

* `-s|--sample-sheet`: Sample sheet filename, tab-delimited text
    format with SampleID and TagRead (barcode) columns
* `-1|--read1`: Read 1 filename, fastq[.gz] format
* `-2|--read2`: Read 2 pair filename, fastq[.gz] format
  (must be consistent with Read 1 filename) (optional)
  If provided then the read files should have read pairs in
  corresponding positions.
* `-m|--mismatches`: Number of mismatches permitted in barcode
    (optional, default 1)
* `-o|--outdir`: Output directory (optional, default output)
* `-d|--delimiter`: Barcode delimiter (optional, default "_")

Outputs:

* If `1` only was provided:
  - A file SampleID.fastq[.gz] with assigned reads.
  - A file, Unassigned.fastq[.gz], with information on unassigned reads.
* If `1` and `2` were provided:
  - Files SampleID_R1.fastq[.gz] and SampleID_R2.fastq[.gz] with assigned
    reads.
  - Files, Unassigned_R1.fastq[.gz] and Unassigned_R2.fastq[.gz] with
    information on unassigned reads.
* If the input file(s) had were of type fastq.gz then the output files
  will be of type fastq.gz.
* A file, num_reads.tsv, with SampleID, TagRead and NumReads columns,
  specifying the number of reads for each SampleID and TagRead in the
  original sample sheet.

Usage:

    python -m riboviz.tools.demultiplex_fastq -h

Known issue:

If the number of mismatches (`-m|--mismatches`) is less than the
Hamming distance between the barcodes (`TagReads` within the sample
sheet) then a read will be assigned to the first barcode that matches
even if this is not the closest barcode in terms of Hamming distance.

For example, imagine we had a barcode in a read, AGA, and our barcodes
in our samplesheet are AAA, CCC, GGG, TTT. The Hamming distances
between the barcode and the sample barcodes are as follows:

* d(AGA, AAA) = 1
* d(AGA, GGG) = 2
* d(AGA, TTT) = 3
* d(AGA, CCC) = 3

If mismatches is 2 then AGA could be assigned to AAA or GGG, depending
on the ordering of barcodes in the sample sheet, even though AAA is
closest in terms of Hamming distance.

If mismatches is 3 then AGA could be assigned to AAA, GGG, TTT
depending on the ordering of barcodes in the sample sheet.

Caution should be taken if the Hamming distance of the barcodes in the
sample sheet is less than the number of mismatches times 2. In the
above two examples, the Hamming distance between each of the sample
barcodes is 3 is less than the number of mismatches times 2, which is
4 and 6 respectively.
"""

import argparse
import gzip
import os
from itertools import islice
from riboviz import utils
from riboviz.utils import BARCODE_DELIMITER
from riboviz.utils import SAMPLE_ID
from riboviz.utils import TAG_READ
from riboviz.utils import NUM_READS


NUM_READS_FILE = "num_reads.tsv"
""" Number of reads summary file name """


def assign_sample(fastq_record1,
                  fastq_record2,
                  barcode,
                  read1_split_fh,
                  read2_split_fh,
                  is_paired_end,
                  mismatches,
                  delimiter):
    """
    Check if fastq record matches barcode for a sample and, if so, add
    the record to the fastq output file(s) for that sample.

    :param fastq_record1: fastq record
    :type fastq_record1: list(str or unicode)
    :param fastq_record2: fastq record for paired read (or None if
    none)
    :type fastq_record2: list(str or unicode)
    :param barcode: Barcode to match fastq record against
    :type barcode: str or unicode
    :param read1_split_fh: Read 1 output file handle
    :type read1_split_fh: io.IOBase
    :param read2_split_fh: Read 2 output file handle (None if none)
    :type read2_split_fh: io.IOBase
    :param is_paired_end: Are paired reads being used? (if True then
    fastq_record2 is assumed to have a fastq record and read2_split_fh
    is assumed to have an output file handle)
    :type is_paired_end: bool
    :param mismatches: Number of mismatches permitted in barcode
    :type mismatches: int
    :param delimiter: Barcode delimiter
    :type delimiter: str or unicode
    :returns: True if fastq record matches barcode
    :rtype: Bool
    """
    is_assigned = False
    if utils.barcode_matches(fastq_record1[0], barcode, mismatches, delimiter):
        is_assigned = True
        read1_split_fh.writelines(fastq_record1)
        if is_paired_end:
            read2_split_fh.writelines(fastq_record2)
    return is_assigned


def assign_samples(fastq_record1,
                   fastq_record2,
                   barcodes,
                   read1_split_fhs,
                   read2_split_fhs,
                   is_paired_end,
                   num_samples,
                   num_reads,
                   mismatches,
                   delimiter):
    """
    Iterate through each sample and check if fastq record matches
    barcode for a sample and, if so, add the record to the fastq
    output file(s) for that barcode, and update the count in num_reads
    for the sample.

    :param fastq_record1: fastq record
    :type fastq_record1: list(str or unicode)
    :param fastq_record2: fastq record for paired read (or None)
    :type fastq_record2: list(str or unicode)
    :param barcodes: Barcodes to match fastq record against
    :type barcodes: list(str or unicode)
    :param read1_split_fhs: Read 1 output file handles (assumed to be
    of the same length as barcodes)
    :type read1_split_fhs: list(io.IOBase)
    :param read2_split_fhs: Read 2 output file handles (None if none)
    :type read2_split_fhs: list(io.IOBase)
    :param is_paired_end: Are paired reads being used? (if True then
    fastq_record2 is assumed to have a fastq record and
    read2_split_fhs is assumed to have output file handles and be of
    the same length as read1_split_fhs)
    :type is_paired_end: bool
    :param num_samples: Number of samples
    :type num_samples: int
    :param num_reads: Number of reads matching barcodes for each
    sample
    :type num_reads: list(int)
    :param mismatches: Number of mismatches permitted in barcode
    :type mismatches: int
    :param delimiter: Barcode delimiter
    :type delimiter: str or unicode
    :returns: True if fastq record matches barcode for a sample
    :rtype: Bool
    """
    is_assigned = False
    for sample in range(num_samples):
        read2_split_fh = None
        if is_paired_end:
            read2_split_fh = read2_split_fhs[sample]
        is_assigned = assign_sample(fastq_record1,
                                    fastq_record2,
                                    barcodes[sample],
                                    read1_split_fhs[sample],
                                    read2_split_fh,
                                    is_paired_end,
                                    mismatches,
                                    delimiter)
        if is_assigned:
            num_reads[sample] += 1
            break
    return is_assigned


def demultiplex(sample_sheet_file,
                read1_file,
                read2_file=None,
                mismatches=1,
                out_dir="output",
                delimiter=BARCODE_DELIMITER):
    """
    Demultiplex reads from fastq[.gz] by inline barcodes.

    :param sample_sheet_file: Sample sheet filename, tab-delimited
    text format with SampleID and TagRead columns
    :type sample_sheet_file: str or unicode
    :param read1_file: Read 1 filename (fastq[.gz] format)
    :type read1_file: str or unicode
    :param read2_file: Read 2 pair filename (fastq[.gz] format,
    must be the same as read1_file format)
    :type read2_file: str or unicode
    :param mismatches: Number of mismatches permitted in barcode
    :type mismatches: int
    :param out_dir: Output directory
    :type out_dir: str or unicode
    :param delimiter: Barcode delimiter
    :type delimiter: str or unicode
    """
    print(("Demultiplexing reads for file: " + read1_file))
    print(("Using sample sheet: " + sample_sheet_file))

    sample_sheet = utils.load_sample_sheet(sample_sheet_file)
    num_samples = sample_sheet.shape[0]
    num_reads = [0] * num_samples
    num_unassigned_reads = 0
    total_reads = 0
    sample_ids = list(sample_sheet[SAMPLE_ID])
    barcodes = list(sample_sheet[TAG_READ])
    length_tag = len(barcodes[1])
    print(("Number of samples: {}".format(num_samples)))
    print(("Allowed mismatches: {}".format(mismatches)))
    print(("Tag length: {}".format(length_tag)))
    print(("Barcode delimiter: {}".format(delimiter)))

    if not os.path.isfile(read1_file):
        raise FileNotFoundError(
            "Error: read 1 file {} does not exist".format(read1_file))

    file_type = os.path.splitext(read1_file)[1]
    is_gz = file_type.lower().endswith(".gz")
    if is_gz:
        file_type = ".fastq.gz"
        open_file = gzip.open
    else:
        file_type = ".fastq"
        open_file = open
    sample_format = "{}" + file_type
    unassigned_format = "Unassigned{}" + file_type

    read1_fh = open_file(read1_file, 'rt')
    is_paired_end = read2_file is not None
    if is_paired_end:
        if not os.path.isfile(read2_file):
            raise FileNotFoundError(
                "Error: read 2 file {} does not exist".format(
                    read2_file))
        read2_fh = open_file(read2_file, 'rt')

    if not os.path.exists(out_dir):
        try:
            os.mkdir(out_dir)
        except Exception:
            raise IOError(
                "Error: output directory {} cannot be created".format(out_dir))

    num_reads_file = os.path.join(out_dir, NUM_READS_FILE)
    if is_paired_end:
        read1_split_fhs = [
            open_file(os.path.join(out_dir,
                                   sample_format.format(sample_id + "_R1")),
                      "wt")
            for sample_id in sample_ids]
        read1_unassigned_fh = open_file(
            os.path.join(out_dir, unassigned_format.format("_R1")), "wt")
        read2_split_fhs = [
            open_file(os.path.join(out_dir,
                                   sample_format.format(sample_id + "_R2")),
                      "wt")
            for sample_id in sample_ids]
        read2_unassigned_fh = open_file(
            os.path.join(out_dir, unassigned_format.format("_R2")), "wt")
    else:
        read1_split_fhs = [
            open_file(os.path.join(out_dir,
                                   sample_format.format(sample_id)),
                      "wt")
            for sample_id in sample_ids]
        read1_unassigned_fh = open_file(
            os.path.join(out_dir, unassigned_format.format("")), "wt")
        read2_split_fhs = None

    while True:
        # Get fastq record/read (4 lines)
        fastq_record1 = list(islice(read1_fh, 4))
        if not fastq_record1:
            break
        if is_paired_end:
            fastq_record2 = list(islice(read2_fh, 4))
        else:
            fastq_record2 = None
        # Count number of processed reads, output every millionth.
        total_reads += 1
        if (total_reads % 1000000) == 0:
            print(("{} reads processed".format(total_reads)))
        # Assign read to a SampleID,
        # TagRead is 1st read with less than threshold mismatches.
        # Beware: this could cause problems if many mismatches.
        is_assigned = assign_samples(fastq_record1,
                                     fastq_record2,
                                     barcodes,
                                     read1_split_fhs,
                                     read2_split_fhs,
                                     is_paired_end,
                                     num_samples,
                                     num_reads,
                                     mismatches,
                                     delimiter)
        if not is_assigned:
            # Write unassigned read to file.
            # Note: unassigned reads are not trimmed.
            read1_unassigned_fh.writelines(fastq_record1)
            if is_paired_end:
                read2_unassigned_fh.writelines(fastq_record2)
            num_unassigned_reads += 1

    # Close output handles and fastq file.
    for fh in read1_split_fhs:
        fh.close()
    read1_unassigned_fh.close()
    read1_fh.close()
    if is_paired_end:
        for fh in read2_split_fhs:
            fh.close()
        read2_unassigned_fh.close()
        read2_fh.close()

    print(("All {} reads processed".format(total_reads)))

    # Output number of reads by sample to file.
    sample_sheet[NUM_READS] = num_reads
    utils.save_deplexed_sample_sheet(sample_sheet,
                                     num_unassigned_reads,
                                     num_reads_file)
    print(("Done"))


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Demultiplex reads from fastq[.gz] by inline barcodes")
    parser.add_argument("-s",
                        "--sample-sheet",
                        dest="sample_sheet_file",
                        nargs='?',
                        help="Sample sheet filename, tab-delimited text format with SampleID and TagRead columns")
    parser.add_argument("-1",
                        "--read1",
                        dest="read1_file",
                        nargs='?',
                        help="Read 1 filename, fastq[.gz] format")
    parser.add_argument("-2",
                        "--read2",
                        dest="read2_file",
                        default=None,
                        nargs='?',
                        help="Read 2 pair filename, fastq[.gz] format")
    parser.add_argument("-m",
                        "--mismatches",
                        dest="mismatches",
                        default=1,
                        type=int,
                        help="Number of mismatches permitted in barcode")
    parser.add_argument("-o",
                        "--outdir",
                        dest="out_dir",
                        nargs='?',
                        default="output",
                        help="Output directory")
    parser.add_argument("-d",
                        "--delimiter",
                        dest="delimiter",
                        nargs='?',
                        default=BARCODE_DELIMITER,
                        help="Barcode delimiter")
    options = parser.parse_args()
    return options


def main():
    """
    Parse command-line options then invoke "demultiplex".
    """
    options = parse_command_line_options()
    sample_sheet_file = options.sample_sheet_file
    read1_file = options.read1_file
    read2_file = options.read2_file
    mismatches = options.mismatches
    out_dir = options.out_dir
    delimiter = options.delimiter
    demultiplex(sample_sheet_file,
                read1_file,
                read2_file,
                mismatches,
                out_dir,
                delimiter)


if __name__ == "__main__":
    main()
