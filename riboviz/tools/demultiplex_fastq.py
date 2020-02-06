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
from riboviz import barcodes_umis
from riboviz import demultiplex_fastq
from riboviz import provenance


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
                        required=True,
                        help="Sample sheet filename, tab-delimited text format with SampleID and TagRead columns")
    parser.add_argument("-1",
                        "--read1",
                        dest="read1_file",
                        required=True,
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
                        default=barcodes_umis.BARCODE_DELIMITER,
                        help="Barcode delimiter")
    options = parser.parse_args()
    return options


def invoke_demultiplex_fastq():
    """
    Parse command-line options then invoke "demultiplex".
    """
    print(provenance.get_provenance_str(__file__))
    options = parse_command_line_options()
    sample_sheet_file = options.sample_sheet_file
    read1_file = options.read1_file
    read2_file = options.read2_file
    mismatches = options.mismatches
    out_dir = options.out_dir
    delimiter = options.delimiter
    demultiplex_fastq.demultiplex(sample_sheet_file,
                                  read1_file,
                                  read2_file,
                                  mismatches,
                                  out_dir,
                                  delimiter)


if __name__ == "__main__":
    invoke_demultiplex_fastq()
