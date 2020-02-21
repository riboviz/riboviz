"""
Process a workflow files log file and count the number of reads
(sequences) processed by specific stages of a RiboViz workflow.

Usage:

    count_reads.py [-h] -i WORKFLOW_FILES_LOG_FILE -o READS_FILE

Arguments:

* '-h', '--help': show this help message and exit
* '-i WORKFLOW_FILES_LOG_FILE', '--input WORKFLOW_FILES_LOG_FILE':
  workflow files log file (input)
* '-o READS_FILE', '--output READS_FILE': reads file (output)

The input file is assumed to be a TSV file with
riboviz.workflow_files_logger-consistent columns:

* 'SampleName': Name of the sample to which this file belongs. This is
  an empty value if the step was not sample-specific (e.g. creating
  index files or demultiplexing a multiplexed FASTQ file).
* 'Program': Program that read/wrote the file. The special token
  'input' denotes input files.
* 'File': Path to file read/written.
* 'Read/Write': 'read' if the file was read, 'write' if the file was
  written.

The files logged in the workflow files log file must exist.

The following information is included:

* Input files: number of reads in the FASTQ files used as inputs.
* 'cutadapt': number of reads in the FASTQ file output.
* 'riboviz.tools.demultiplex_fastq': number of reads in the FASTQ
  files output, as recorded in the 'num_reads.tsv' file output.
* 'hisat2': number of reads in the SAM file and FASTQ file output.
* 'riboviz.tools.trim_5p_mismatch': number of reads in the SAM fil
  output as recorded in the TSV summary file output.
* 'umi_tools dedup': number of reads in the BAM file output.

The output file is a TSV file with columns:

* 'SampleName'
* 'Program'
* 'File'
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
        description="Process a workflow files log file and count the number of reads (sequences) processed at specific stages of a RiboViz workflow")
    parser.add_argument("-i",
                        "--input",
                        dest="workflow_file",
                        required=True,
                        help="Workflow files log file (input)")
    parser.add_argument("-o",
                        "--output",
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
    workflow_file = options.workflow_file
    reads_file = options.reads_file
    count_reads.count_reads(workflow_file, reads_file)


if __name__ == "__main__":
    invoke_count_reads()
