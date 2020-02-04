"""
Scan input, temporary and output directories and count the number of
reads (sequences) processed by specific stages of a RiboViz
workflow. The scan is based on the directory structure and file
patterns used by RiboViz.

Usage:

        count_reads.py [-h] -i INPUT_DIR -t TMP_DIR -o OUTPUT_DIR \
                       -r READS_FILE

Arguments:

* '-h', '--help':  show this help message and exit
* '-i INPUT_DIR', '--input-dir INPUT_DIR': Input files directory
* '-t TMP_DIR', '--tmp-dir TMP_DIR': Temporary files directory
* '-o OUTPUT_DIR', '--output OUTPUT_DIR': Output files directory
* '-r READS_FILE', '--reads-file READS_FILE': Reads file (output)

The following information is included:

* Input files: number of reads in the FASTQ files used as inputs.
* 'cutadapt': number of reads in the FASTQ file output.
* 'riboviz.tools.demultiplex_fastq': number of reads in the FASTQ
  files output, as recorded in the 'num_reads.tsv' TSV file output.
* 'hisat2': number of reads in the SAM file and FASTQ file output.
* 'riboviz.tools.trim_5p_mismatch': number of reads in the SAM file
  output as recorded in the 'trim_5p_mismatch.tsv' summary file
  output, or the SAM file itself, if the TSV file cannot be found.
* 'umi_tools dedup': number of reads in the BAM file output.

The output file is a TSV file with columns:

* 'SampleName': Name of the sample to which this file belongs. This is
  an empty value if the step was not sample-specific
  (e.g. demultiplexing a multiplexed FASTQ file).
* 'Program': Program that wrote the file. The special token
  'input' denotes input files.
* 'File': Path to file.
* 'NumReads': Number of reads in the file.
* 'Description': Human-readable description of the file contents.
"""
import argparse
from riboviz import count_reads
from riboviz import provenance


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Scan RiboViz input, temporary and output directories and count the number of reads (sequences) processed at specific stages of a RiboViz workflow")
    parser.add_argument("-i",
                        "--input-dir",
                        dest="input_dir",
                        required=True,
                        help="Input files directory")
    parser.add_argument("-t",
                        "--tmp-dir",
                        dest="tmp_dir",
                        required=True,
                        help="Temporary files directory")
    parser.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="Output files directory")
    parser.add_argument("-r",
                        "--reads-file",
                        dest="reads_file",
                        required=True,
                        help="Reads file (output)")
    options = parser.parse_args()
    return options


def invoke_count_reads():
    """
    Parse command-line options then invoke "count_reads".
    """
    print(provenance.get_provenance_str(__file__))
    options = parse_command_line_options()
    input_dir = options.input_dir
    tmp_dir = options.tmp_dir
    output_dir = options.output_dir
    reads_file = options.reads_file
    count_reads.count_reads(input_dir, tmp_dir, output_dir, reads_file)


if __name__ == "__main__":
    invoke_count_reads()
