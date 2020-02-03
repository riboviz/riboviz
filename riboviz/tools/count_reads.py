"""
Process a workflow files log file and count the number of reads
(sequences) processed at specific stages of the
'riboviz.tools.prep_riboviz' workflow.

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

The output file is a TSV file with columns:

* SampleName
* Program
* File
* NumReads: Number of reads in the file.
* Description
"""
import argparse
import pandas as pd
from riboviz import fastq
from riboviz import provenance
from riboviz import sam_bam
from riboviz import workflow_files_logger


NUM_READS = "NumReads"
""" NumReads column name """
HEADER = [workflow_files_logger.SAMPLE_NAME,
          workflow_files_logger.PROGRAM,
          workflow_files_logger.FILE,
          NUM_READS,
          workflow_files_logger.DESCRIPTION]
""" Number of reads file header """


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


def count_reads_inputs(df):
    """
    Extract "input" rows from DataFrame, count number of reads in
    each FASTQ file used as an input, and return list of Series, one
    per row, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File, NumReads,
    Description columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    associated FASTQ files
    :rtype: list(pandas.core.frame.Series)
    """
    input_df = df[
        df[workflow_files_logger.PROGRAM] == workflow_files_logger.INPUT
    ]
    if input_df.empty:
        return []
    files_df = input_df[input_df[
        workflow_files_logger.FILE].str.lower().str.endswith(tuple(fastq.FASTQ_EXTS))]
    if files_df.empty:
        return []
    count_rows = []
    for _, row in files_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        row[NUM_READS] = fastq.count_sequences(file_name)
        count_rows.append(row)
    return count_rows


def count_reads_cutadapt(df):
    """
    Extract "cutadapt" rows from DataFrame, count number of reads
    in each FASTQ file output, and return list of Series, one per row,
    with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File, NumReads,
    Description columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] == "cutadapt"]
    if program_df.empty:
        return []
    # cutadapt outputs one file, the FASTQ file.
    row = list(program_df.iterrows())[0][1]
#    row = program_df.iloc[0]
    file_name = row[workflow_files_logger.FILE]
    row[NUM_READS] = fastq.count_sequences(file_name)
    return [row]


def count_reads_demultiplex_fastq(df):
    """
    Extract "riboviz.tools.demultiplex_fastq" rows from DataFrame,
    count number of reads in each demultiplexed FASTQ file output, and
    return list of Series, one per row, with NumReads set to the
    number of reads in the files.

    :param df: DataFrame with SampleName, Program, File, NumReads,
    Description columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    # TODO
    # Number of reads in demultiplexed FASTQ files and unassigned reads,
    # as recorded in the 'num_reads.tsv' file output. This file is used
    # to save having to traverse each of the output FASTQ files.
    program_df = df[df[workflow_files_logger.PROGRAM] ==
                    "riboviz.tools.demultiplex_fastq"]
    if program_df.empty:
        return []
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(tuple(fastq.FASTQ_EXTS))]
    count_rows = []
    for _, row in files_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        row[NUM_READS] = fastq.count_sequences(file_name)
        count_rows.append(row)
    return count_rows


def count_reads_hisat2(df):
    """
    Extract "hisat2" rows from DataFrame, count number of reads in
    each FASTQ and SAM file output, and return list of Series, one per
    row, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File, NumReads,
    Description columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] == "hisat2"]
    if program_df.empty:
        return []
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith("sam")]
    count_rows = []
    for _, row in files_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        # TODO handle both values - (sequences, primary sequences)
        counts = sam_bam.count_sequences(file_name)
        row[NUM_READS] = counts
        count_rows.append(row)
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(tuple(fastq.FASTQ_EXTS))]
    for _, row in files_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        row[NUM_READS] = fastq.count_sequences(file_name)
        count_rows.append(row)
    return count_rows


def count_reads_trim_5p_mismatch(df):
    """
    Extract "riboviz.tools.trim_5p_mismatch" rows from DataFrame,
    count number of reads in the SAM file output, and return list of
    Series, one per row, with NumReads set to the number of reads in
    the files.

    :param df: DataFrame with SampleName, Program, File, NumReads,
    Description columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    # TODO
    # Number of reads in SAM file output as recorded in the TSV summary
    # file output.
    # This file is used to save having to traverse each of
    # the output SAM files.
    program_df = df[df[workflow_files_logger.PROGRAM] ==
                    "riboviz.tools.trim_5p_mismatch"]
    if program_df.empty:
        return []
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith("sam")]
    # trim_5p_mismatch outputs one SAM file.
    row = list(files_df.iterrows())[0][1]
    file_name = row[workflow_files_logger.FILE]
    # TODO handle both values - (sequences, primary sequences)
    counts = sam_bam.count_sequences(file_name)
    row[NUM_READS] = counts
    return [row]


def count_reads_umi_tools_dedup(df):
    """
    Extract "umi_tools dedup" program_df from DataFrame, count number of
    reads in the BAM file output, and return list of Series, one per
    row, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File, NumReads,
    Description columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] ==
                    "umi_tools dedup"]
    if program_df.empty:
        return []
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith("bam")]
    # umi_tools dedup outputs one BAM file.
    row = list(files_df.iterrows())[0][1]
    file_name = row[workflow_files_logger.FILE]
    # TODO handle both values - (sequences, primary sequences)
    counts = sam_bam.count_sequences(file_name)
    row[NUM_READS] = counts
    return [row]


def count_input_write_reads(df):
    """

    df is assumed to be a Pandas DataFrame with
    riboviz.workflow_files_logger-consistent columns:

    * 'SampleName': Name of the sample to which this file
      belongs. This is an empty value if the step was not
      sample-specific (e.g. creating index files or demultiplexing a
      multiplexed FASTQ file).
    * 'Program': Program that read/wrote the file. The special token
      'input' denotes input files.
    * 'File': Path to file read/written.
    * 'Read/Write': 'read' if the file was read, 'write' if the file
      was written.
    * 'Description': Human-readable description of the step at which
      this file was read or written.

    A DataFrame with the following columns is returned:

    * SampleName
    * Program
    * File
    * NumReads: Number of reads in the file.
    * Description

    :param df: DataFrame with riboviz.workflow_files_logger-consistent
    data
    :type df: pandas.core.frame.DataFrame
    :return: DataFrame with number of reads in each input file and
    output file
    :rtype: pandas.core.frame.DataFrame
    """
    # Filter to get files recorded as "input" or "write" only.
    logs_df = df[
        (df[workflow_files_logger.PROGRAM] == workflow_files_logger.INPUT)
        |
        (df[workflow_files_logger.READ_WRITE] == workflow_files_logger.WRITE)
    ]
    reads_df = pd.DataFrame(columns=HEADER)
    # Replace Read/Write column with NumReads/
    index = logs_df.columns.get_loc(workflow_files_logger.READ_WRITE)
    logs_df = logs_df.drop(columns=workflow_files_logger.READ_WRITE)
    logs_df.insert(index, NUM_READS, 0)
    # Note: steps not corresponding to a sample will have SAMPLE_NAME
    # == "".
    samples = logs_df[workflow_files_logger.SAMPLE_NAME].unique()
    for sample in samples:
        sample_df = logs_df[logs_df[workflow_files_logger.SAMPLE_NAME] == sample]
        rows = count_reads_inputs(sample_df)
        if rows:
            reads_df = reads_df.append(rows, ignore_index=True)
        rows = count_reads_cutadapt(sample_df)
        if rows:
            reads_df = reads_df.append(rows, ignore_index=True)
        rows = count_reads_demultiplex_fastq(sample_df)
        if rows:
            reads_df = reads_df.append(rows, ignore_index=True)
        rows = count_reads_hisat2(sample_df)
        if rows:
            reads_df = reads_df.append(rows, ignore_index=True)
        rows = count_reads_trim_5p_mismatch(sample_df)
        if rows:
            reads_df = reads_df.append(rows, ignore_index=True)
        rows = count_reads_umi_tools_dedup(sample_df)
        if rows:
            reads_df = reads_df.append(rows, ignore_index=True)
    return reads_df


def count_reads(workflow_file, reads_file, delimiter="\t", comment="#"):
    """
    Process workflow files log file and count reads in FASTQ, SAM and
    BAM files that were input or written and save results.

    workflow_file is assumed to be a TSV file with
    riboviz.workflow_files_logger-consistent columns:

    * 'SampleName': Name of the sample to which this file
      belongs. This is an empty value if the step was not
      sample-specific (e.g. creating index files or demultiplexing a
      multiplexed FASTQ file).
    * 'Program': Program that read/wrote the file. The special token
      'input' denotes input files.
    * 'File': Path to file read/written.
    * 'Read/Write': 'read' if the file was read, 'write' if the file
      was written.
    * 'Description': Human-readable description of the step at which
      this file was read or written.

    reads_file is a TSV file with columns:

    * SampleName
    * Program
    * File
    * NumReads: Number of reads in the file.
    * Description

    :param workflow_file: Workflow files log file input
    :type workflow_file: str or unicode
    :param reads_file: Reads file output
    :type reads_file: str or unicode
    :param comment: Comment prefix
    :type comment: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    """
    df = pd.read_csv(workflow_file, sep=delimiter, comment=comment)
    # Pandas treats missing cells as nan, so convert to "".
    for column in workflow_files_logger.HEADER:
        df.fillna(value={column: ""}, inplace=True)
    reads_df = count_input_write_reads(df)
    provenance.write_provenance_header(__file__, reads_file)
    reads_df[list(reads_df.columns)].to_csv(
        reads_file, mode='a', sep=delimiter, index=False)


def invoke_count_reads():
    """
    Parse command-line options then invoke "count_reads".
    """
    options = parse_command_line_options()
    workflow_file = options.workflow_file
    reads_file = options.reads_file
    count_reads(workflow_file, reads_file)


if __name__ == "__main__":
    invoke_count_reads()
