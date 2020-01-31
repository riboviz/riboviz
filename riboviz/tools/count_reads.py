"""
Process a workflow files log file and count the number of reads (sequences)
processed at specific stages of the 'riboviz.tools.prep_riboviz' workflow

Usage:

    count_reads.py [-h] -i WORKFLOW_FILES_LOG_FILE -o READS_FILE

Arguments:

* '-h', '--help': show this help message and exit
* '-i WORKFLOW_FILES_LOG_FILE', '--input WORKFLOW_FILES_LOG_FILE':
  workflow files log file (input)
* '-o READS_FILE', '--output READS_FILE': reads file (output)

The following information is extracted:

* Input files
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
"""
import argparse
import os
import sys
import pandas as pd
from riboviz import fastq
from riboviz import provenance
from riboviz import sam_bam
from riboviz import workflow_files_logger


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Process a workflow files log file and count the number of reads (sequences) processed at specific stages of the 'riboviz.tools.prep_riboviz' workflow")
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


def get_file_ext(file_name):
    """
    Given a file name return full file extension, everything after the
    first "." in the file name. For example, for 'a.fastq.gz' return
    'fastq.gz', for 'a.fastq' return 'fastq', for 'a' return ''. The
    extension is returned in lower-case.

    :param file_name: File name
    :type file_name: str or unicode
    :return: extension
    :rtype: str or unicode
    """
    file_type = ".".join(os.path.basename(file_name).split(".")[1:])
    return file_type.lower()


FASTQ_EXTS = ["fastq", "fq", "fastq.gz", "fq.gz", "fastq.gzip", "fq.gzip"]


def filter_file_type(df):
    print("XXX " + df)
    print(type(df))
    print(df[workflow_files_logger.FILE])
    file_name = df[workflow_files_logger.FILE]
    file_type = ".".join(os.path.basename(file_name).split(".")[1:])
    print(file_type)
    return file_type in FASTQ_EXTS


def get_program_files(df, program, read_write):
    filtered = df[(df[workflow_files_logger.PROGRAM] == program)
                  & (df[workflow_files_logger.READ_WRITE] == read_write)]
    projected = filtered[[workflow_files_logger.SAMPLE_NAME,
                          workflow_files_logger.FILE]]
    return projected


def count_reads_sample(df):
    program = "cutadapt"
    print("-------- {} --------".format(program))
    cutadapt = get_program_files(df, program,
                                 workflow_files_logger.WRITE)
    # print(cutadapt)
    if len(cutadapt) > 0:
        # cutadapt outputs one file, the FASTQ file.
        file_name = cutadapt.iloc[0][workflow_files_logger.FILE]
        counts = 0
        # TODO uncomment
        # counts = fastq.count_sequences(file_name)
        print("{}: {}".format(file_name, counts))
    else:
        print("No entries")

    program = "umi_tools dedup"
    print("-------- {} --------".format(program))
    umi_tools_dedup = get_program_files(df, program,
                                        workflow_files_logger.WRITE)
    # print(umi_tools_dedup)
    bams = umi_tools_dedup[umi_tools_dedup[
        workflow_files_logger.FILE].str.lower().str.endswith("bam")]
    for _, row in bams.iterrows():
        file_name = row[workflow_files_logger.FILE]
        counts = 0
        # TODO uncomment
        # counts = sam_bam.count_sequences(file_name)
        # TODO handle both values - (sequences, primary sequences)
        print("{}: {}".format(file_name, counts))

    program = "hisat2"
    print("-------- {} --------".format(program))
    hisat2 = get_program_files(df, program, workflow_files_logger.WRITE)
    # print(hisat2)
    sams = hisat2[hisat2[
        workflow_files_logger.FILE].str.lower().str.endswith("sam")]
    for _, row in sams.iterrows():
        file_name = row[workflow_files_logger.FILE]
        counts = 0
        # TODO uncomment
        # counts = sam_bam.count_sequences(file_name)
        # TODO handle both values - (sequences, primary sequences)
        print("{}: {}".format(file_name, counts))
    fqs = hisat2[hisat2[
        workflow_files_logger.FILE].str.lower().str.endswith(tuple(FASTQ_EXTS))]
    for _, row in fqs.iterrows():
        file_name = row[workflow_files_logger.FILE]
        counts = 0
        # TODO uncomment
        # counts = fastq.count_sequences(file_name)
        print("{}: {}".format(file_name, counts))

    # TODO
    # * Number of reads in demultiplexed FASTQ files and unassigned reads,
    #   as recorded in the 'num_reads.tsv' file output. This file is used
    #   to save having to traverse each of the output FASTQ files.
    # * Compare to explicit count of reads in .fq files.
    program = "riboviz.tools.demultiplex_fastq"
    print("-------- {} --------".format(program))
    demultiplex_fastq = get_program_files(df, program, workflow_files_logger.WRITE)
    print(demultiplex_fastq)

    # TODO
    # * Number of reads in SAM file output as recorded in the TSV summary
    #   file output..
    # * Compare to explicit count of reads in .sam file.
    program = "riboviz.tools.trim_5p_mismatch"
    print("-------- {} --------".format(program))
    trim_5p_mismatch = get_program_files(df, program, workflow_files_logger.WRITE)
    print(trim_5p_mismatch)

    # TODO
    # * Number of reads in sample FASTQ files used as inputs (if
    #   non-multiplexed samples are used).
    # * Number of reads in multiplexed FASTQ file (if multiplexed
    #   samples are used).
    program = "inputs"
    print("-------- {} --------".format(program))
    inputs = df[df[workflow_files_logger.PROGRAM] ==
                workflow_files_logger.INPUT][[workflow_files_logger.SAMPLE_NAME,
                                              workflow_files_logger.FILE]]
    print(inputs)


def count_reads(workflow_file, reads_file):
    df = pd.read_csv(workflow_file, sep="\t", comment="#")
    # Pandas treats missing cells as nan, so convert to ""
    for column in workflow_files_logger.HEADER:
        df.fillna(value={column: ""}, inplace=True)
    print("----- samples")
    samples = df[workflow_files_logger.SAMPLE_NAME].unique()
    print(samples)
    for sample in samples:
        print("========{}========".format(sample))
        print("========{}========".format(sample))
        print("========{}========".format(sample))
        sample_df = df[df[workflow_files_logger.SAMPLE_NAME] == sample]
        print(sample_df)
        count_reads_sample(sample_df)
    # https://cmdlinetips.com/2018/02/how-to-subset-pandas-dataframe-based-on-values-of-a-column/
    # .fastq, .fastq.gz, .fq, .fq.gz, .fq.gzip, .fastq.gzip
    # .sam, .bam, .tsv

    # TODO - Output TSV file:
    #
    # SampleName | Human summary | Program | File | Count |
    #
    # Include a provenance header comment.


def invoke_count_reads():
    """
    Parse command-line options then invoke "count_reads".
    """
    print(provenance.get_provenance_str(__file__))
    options = parse_command_line_options()
    workflow_file = options.workflow_file
    reads_file = options.reads_file
    count_reads(workflow_file, reads_file)


if __name__ == "__main__":
    invoke_count_reads()
