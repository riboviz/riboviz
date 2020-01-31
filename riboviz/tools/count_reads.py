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


# TODO Move to riboviz.fastq.
FASTQ_EXTS = ["fastq", "fq", "fastq.gz", "fq.gz", "fastq.gzip", "fq.gzip"]


# TODO Move to riboviz.utils and add tests.
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


def get_program_files(df, program, read_write):
    """
    Filter data to extract rows whose "Program" and "Read/Write"
    column values match those given.

    :param df: Data to filter
    :type df: pandas.core.frame.DataFrame
    :param program: Value of "Program" column to filter on
    :type program: str or unicode
    :param read_write: Value of "Read/Write" column to filter on
    :type read_write: str or unicode
    :rtype: pandas.core.frame.DataFrame
    """
    filtered = df[
        (df[workflow_files_logger.PROGRAM] == program)
        & (df[workflow_files_logger.READ_WRITE] == read_write)
    ]
    return filtered


# TODO rename
def log_entry(row, counts):
    """
    Log an entry.

    :param row: Data
    :type row: pandas.core.frame.Series
    :param counts: Counts
    :type counts: int
    """
    sample = row[workflow_files_logger.SAMPLE_NAME]
    program = row[workflow_files_logger.PROGRAM]
    file_name = row[workflow_files_logger.FILE]
    print("{} {} {} {}".format(sample, program, file_name, counts))


def count_reads_inputs(df):
    program = "inputs"
    entries = df[df[workflow_files_logger.PROGRAM] ==
                 workflow_files_logger.INPUT]
    if len(entries) == 0:
        return
    fqs = entries[entries[
        workflow_files_logger.FILE].str.lower().str.endswith(tuple(FASTQ_EXTS))]
    if len(fqs) == 0:
        return
    print("-------- {} --------".format(program))
    for _, row in fqs.iterrows():
        file_name = row[workflow_files_logger.FILE]
        counts = fastq.count_sequences(file_name)
        log_entry(row, counts)


def count_reads_cutadapt(df):
    program = "cutadapt"
    entries = get_program_files(df, program,
                                workflow_files_logger.WRITE)
    if len(entries) == 0:
        return
    print("-------- {} --------".format(program))
    # cutadapt outputs one file, the FASTQ file.
    row = entries.iloc[0]
    file_name = row[workflow_files_logger.FILE]
    counts = fastq.count_sequences(file_name)
    log_entry(row, counts)


def count_reads_demultiplex_fastq(df):
    # TODO
    # Number of reads in demultiplexed FASTQ files and unassigned reads,
    # as recorded in the 'num_reads.tsv' file output. This file is used
    # to save having to traverse each of the output FASTQ files.
    program = "riboviz.tools.demultiplex_fastq"
    entries = get_program_files(df, program,
                                workflow_files_logger.WRITE)
    if len(entries) == 0:
        return
    print("-------- {} --------".format(program))
    fqs = entries[entries[
        workflow_files_logger.FILE].str.lower().str.endswith(tuple(FASTQ_EXTS))]
    for _, row in fqs.iterrows():
        file_name = row[workflow_files_logger.FILE]
        counts = fastq.count_sequences(file_name)
        log_entry(row, counts)


def count_reads_hisat2(df):
    program = "hisat2"
    entries = get_program_files(df, program,
                                workflow_files_logger.WRITE)
    if len(entries) == 0:
        return
    print("-------- {} --------".format(program))
    sams = entries[entries[
        workflow_files_logger.FILE].str.lower().str.endswith("sam")]
    for _, row in sams.iterrows():
        file_name = row[workflow_files_logger.FILE]
        # TODO handle both values - (sequences, primary sequences)
        counts = sam_bam.count_sequences(file_name)
        log_entry(row, counts)
    fqs = entries[entries[
        workflow_files_logger.FILE].str.lower().str.endswith(tuple(FASTQ_EXTS))]
    for _, row in fqs.iterrows():
        file_name = row[workflow_files_logger.FILE]
        counts = fastq.count_sequences(file_name)
        log_entry(row, counts)


def count_reads_trim_5p_mismatch(df):
    # TODO
    # Number of reads in SAM file output as recorded in the TSV summary
    # file output.
    program = "riboviz.tools.trim_5p_mismatch"
    entries = get_program_files(df, program,
                                workflow_files_logger.WRITE)
    if len(entries) == 0:
        return
    print("-------- {} --------".format(program))
    sams = entries[entries[
        workflow_files_logger.FILE].str.lower().str.endswith("sam")]

    # trim_5p_mismatch outputs one file, the SAM file.
    row = sams.iloc[0]
    file_name = row[workflow_files_logger.FILE]
    # TODO handle both values - (sequences, primary sequences)
    counts = sam_bam.count_sequences(file_name)
    log_entry(row, counts)


def count_reads_umi_tools_dedup(df):
    program = "umi_tools dedup"
    entries = get_program_files(df, program,
                                workflow_files_logger.WRITE)
    if len(entries) == 0:
        return
    print("-------- {} --------".format(program))
    bams = entries[entries[
        workflow_files_logger.FILE].str.lower().str.endswith("bam")]
    for _, row in bams.iterrows():
        file_name = row[workflow_files_logger.FILE]
        counts = sam_bam.count_sequences(file_name)
        # TODO handle both values - (sequences, primary sequences)
        log_entry(row, counts)


def count_reads(workflow_file, reads_file):
    df = pd.read_csv(workflow_file, sep="\t", comment="#")
    # Pandas treats missing cells as nan, so convert to ""
    for column in workflow_files_logger.HEADER:
        df.fillna(value={column: ""}, inplace=True)
    # Note: steps not corresponding to a sample will have SAMPLE_NAME
    # "".
    samples = df[workflow_files_logger.SAMPLE_NAME].unique()
    print("Samples:" + str(samples))
    for sample in samples:
        print("======== {} ========".format(sample))
        print("======== {} ========".format(sample))
        print("======== {} ========".format(sample))
        sample_df = df[df[workflow_files_logger.SAMPLE_NAME] == sample]
        count_reads_inputs(sample_df)
        count_reads_cutadapt(sample_df)
        count_reads_demultiplex_fastq(sample_df)
        count_reads_hisat2(sample_df)
        count_reads_trim_5p_mismatch(sample_df)
        count_reads_umi_tools_dedup(sample_df)

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
