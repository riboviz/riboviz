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
* 'Description': Human-readable description of the step at which this
  file was read or written.

The files logged in the workflow files log file must exist.

Read counts for files produced at the following stages are calculated:

* Input files:
  - Number of reads in sample FASTQ files used as inputs (if
    non-multiplexed samples are used).
  - Number of reads in multiplexed FASTQ file (if multiplexed
    samples are used).
* 'cutadapt':
  - Number of reads in FASTQ file output as recorded in the FASTQ
    file output.
* 'riboviz.tools.demultiplex_fastq'
  - Number of reads in demultiplexed FASTQ files and unassigned reads,
    as recorded in the 'num_reads.tsv' file output. This file is used
    to save having to traverse each of the output FASTQ files.
* 'hisat2':
  - Number of reads in SAM file and FASTQ files output.
* 'riboviz.tools.trim_5p_mismatch':
  - Number of reads in SAM file output as recorded in the TSV summary
    file output.This file is used to save having to traverse each of
    the output SAM files.
* 'umi_tools dedup':
  - Number of reads in output.

The output file is a TSV file with columns:

* SampleName
* Program
* File
* NumReads: Number of reads in the file.
* Description
"""
import argparse
from riboviz import count_reads


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
    options = parse_command_line_options()
    workflow_file = options.workflow_file
    reads_file = options.reads_file
    count_reads.count_reads(workflow_file, reads_file)


if __name__ == "__main__":
    invoke_count_reads()
