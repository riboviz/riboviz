#! python
"""
Demultiplex fastq files using UMI-tools-compliant barcodes present
within the fastq headers. These headers are assumed to be of form:

    @..._<BARCODE>_...

where the barcode is the first section delimited by underscores. If
another delimiter was used then that can be specified.

Inputs:

* `-ss|--samplesheet`: Sample sheet filename, tab-delimited text
    format with SampleID and TagRead (barcode) columns
* `-r1|--read1`: Read 1 filename, fastq.gz format
* `-r2|--read2`: Read 2 pair filename, fastq.gz format (optional)
  If provided then the read files should have read pairs in
  corresponding positions.
* `-m|--mismatches`: Number of mismatches permitted in barcode
    (optional, default 1)
* `-o|--outdir`: Output directory (optional, default output)
* `-d|--delimiter`: Barcode delimiter (optional, default "_")

Outputs:

* If `r1` only was provided:
  - A file SampleID.fastq.gz with assigned reads.
  - A file, Unassigned.fastq.gz, with information on unassigned reads.
* If `r1` and `r2` were provided:
  - Files SampleID_R1.fastq.gz and SampleID_R2.fastq.gz with assigned
    reads.
  - Files, Unassigned_R1.fastq.gz and Unassigned_R2.fastq.gz with
    information on unassigned reads.
* A file, nreads.txt, with SampleID, TagRead and nreads columns,
  specifying the number of reads for each SampleID and TagRead in the
  original sample sheet.

Usage:

    python demultiplex_fastq.py -h
"""

import argparse
import gzip
import os
from itertools import islice
import pandas as pd


DELIMITER = "_"
""" Default barcode delmiter in fastq header """


def barcode_mismatch(record, barcode, mismatches=0, delimiter=DELIMITER):
    """
    Returns True if fastq record header includes barcode, up to a given
    number of mismatches. The header is assumed to be of form:

        @...<DELIMITER><BARCODE><DELIMITER>...

    :param record: Record
    :type record: str or unicode
    :param barcode: Barcode
    :type barcode: str or unicode
    :param mismatches: Number of mismatches
    :type mismatches: int
    :param delimiter: Barcode delimiter
    :type delimiter: str or unicode
    :returns: True or False
    :rtype: bool
    """
    chunks = record.split(delimiter)
    if len(chunks) == 1:
        return False
    candidate = chunks[1]
    mismatch = sum([candidate[i] != barcode[i] for i in range(len(barcode))])
    return mismatch <= mismatches


def demultiplex(sample_sheet_file,
                read1_file,
                read2_file=None,
                mismatches=1,
                out_dir="output",
                delimiter=DELIMITER):
    """
    Demultiplex reads from fastq.gz by inline barcodes.

    :param sample_sheet_file: Sample sheet filename, tab-delimited
    text format with SampleID and TagRead columns
    :type sample_sheet_file: str or unicode
    :param read1_file: Read 1 filename (fastq.gz format)
    :type read1_file: str or unicode
    :param read2_file: Read 2 pair filename (fastq.gz format)
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

    if not os.path.isfile(sample_sheet_file):
        raise FileNotFoundError(
            "Error: sample sheet file {} does not exist".format(
                sample_sheet_file))

    sample_sheet = pd.read_csv(sample_sheet_file,
                               comment="#",
                               delimiter="\t")
    num_samples = sample_sheet.shape[0]
    num_reads = [0] * num_samples
    num_unassigned_reads = 0
    total_reads = 0
    sample_ids = list(sample_sheet.SampleID)
    tag_reads = list(sample_sheet.TagRead)
    length_tag = len(tag_reads[1])
    print(("Number of samples: {}".format(num_samples)))
    print(("Allowed mismatches: {}".format(mismatches)))
    print(("Tag length: {}".format(length_tag)))
    print(("Barcode delimiter: {}".format(delimiter)))

    if not os.path.isfile(read1_file):
        raise FileNotFoundError(
            "Error: read 1 file {} does not exist".format(read1_file))
    read1_fh = gzip.open(read1_file, 'rt')

    is_paired_end = read2_file is not None
    if is_paired_end:
        if not os.path.isfile(read2_file):
            raise FileNotFoundError(
                "Error: read 2 file {} does not exist".format(
                    read2_file))
        read2_fh = gzip.open(read2_file, 'rt')

    try:
        os.mkdir(out_dir)
    except Exception:
        raise IOError(
            "Error: output directory {} cannot be created".format(out_dir))

    num_reads_file = os.path.join(out_dir, "nreads.txt")
    if is_paired_end:
        read1_split_fhs = [
            gzip.open(os.path.join(out_dir, sample_id + "_R1.fastq.gz"), "wt")
            for sample_id in sample_ids]
        read1_unassigned_fh = gzip.open(
            os.path.join(out_dir, "Unassigned_R1.fastq.gz"), "wt")
        read2_split_fhs = [
            gzip.open(os.path.join(out_dir, sample_id + "_R2.fastq.gz"), "wt")
            for sample_id in sample_ids]
        read2_unassigned_fh = gzip.open(
            os.path.join(out_dir, "Unassigned_R2.fastq.gz"), "wt")
    else:
        read1_split_fhs = [
            gzip.open(os.path.join(out_dir, sample_id + ".fastq.gz"), "wt")
            for sample_id in sample_ids]
        read1_unassigned_fh = gzip.open(
            os.path.join(out_dir, "Unassigned.fastq.gz"), "wt")

    while True:
        # Get fastq record/read (4 lines)
        fastq_record1 = list(islice(read1_fh, 4))
        if not fastq_record1:
            break
        if is_paired_end:
            fastq_record2 = list(islice(read2_fh, 4))
        # Count number of processed reads, output every millionth.
        total_reads += 1
        if (total_reads % 1000000) == 0:
            print(("{} reads processed".format(total_reads)))
        # Assign read to a SampleID,
        # TagRead is 1st read with less than threshold mismatches.
        # Beware: this could cause problems if many mismatches.
        is_assigned = False

        for sample in range(num_samples):
            # Check if initial segment of read matches sample barcode.
            if barcode_mismatch(fastq_record1[0],
                                tag_reads[sample],
                                mismatches,
                                delimiter):
                is_assigned = True
                read1_split_fhs[sample].writelines(fastq_record1)
                if is_paired_end:
                    read2_split_fhs[sample].writelines(fastq_record2)
                num_reads[sample] += 1
                # Stop testing against other barcodes
                break
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
    sample_sheet["nreads"] = num_reads
    sample_sheet[['SampleID', 'TagRead', 'nreads']].to_csv(
        num_reads_file, sep="\t", index=False)
    # Append unassigned reads to file.
    with open(num_reads_file, "a") as nrf:
        nrf.write("Unassd\tNNNNNNNNN\t{}\n".format(num_unassigned_reads))
        nrf.write("TOTAL\t\t{}".format(total_reads))

    print(("Done"))


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Demultiplex reads from fastq.gz by inline barcodes")
    parser.add_argument("-ss",
                        "--samplesheet",
                        dest="sample_sheet_file",
                        nargs='?',
                        help="Sample sheet filename, tab-delimited text format with SampleID and TagRead columns")
    parser.add_argument("-r1",
                        "--read1",
                        dest="read1_file",
                        nargs='?',
                        help="Read 1 filename, fastq.gz format")
    parser.add_argument("-r2",
                        "--read2",
                        dest="read2_file",
                        default=None,
                        nargs='?',
                        help="Read 2 pair filename, fastq.gz format")
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
                        default=DELIMITER,
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
