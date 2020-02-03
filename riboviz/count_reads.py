"""
Process a workflow files log file and count the number of reads
(sequences) processed by specific stages of a RiboViz workflow.

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
    file output. This file is used to save having to traverse each of
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
import os
import pandas as pd
from riboviz import demultiplex_fastq
from riboviz import fastq
from riboviz import provenance
from riboviz import sam_bam
from riboviz import sample_sheets
from riboviz import trim_5p_mismatch
from riboviz import workflow_files_logger


NUM_READS = "NumReads"
""" NumReads column name """
HEADER = [workflow_files_logger.SAMPLE_NAME,
          workflow_files_logger.PROGRAM,
          workflow_files_logger.FILE,
          NUM_READS,
          workflow_files_logger.DESCRIPTION]
""" Number of reads file header """


def count_reads_inputs(df):
    """
    Extract "input" rows from DataFrame, count number of reads in
    each FASTQ file used as an input, and return list of Series, one
    per row, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File, Description
    columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    input_df = df[
        df[workflow_files_logger.PROGRAM] == workflow_files_logger.INPUT
    ]
    if input_df.empty:
        return []
    files_df = input_df[input_df[
        workflow_files_logger.FILE].str.lower().str.endswith(
            tuple(fastq.FASTQ_EXTS))]
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

    :param df: DataFrame with SampleName, Program, File, Description
    columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] == "cutadapt"]
    if program_df.empty:
        return []
    # cutadapt outputs one file, the FASTQ file.
    row = program_df.iloc[0]
    file_name = row[workflow_files_logger.FILE]
    row[NUM_READS] = fastq.count_sequences(file_name)
    return [row]


def count_reads_demultiplex_fastq(df):
    """
    Extract "riboviz.tools.demultiplex_fastq" rows from DataFrame,
    count number of reads in each demultiplexed FASTQ file output, as
    recorded in the TSV file output by "demultiplex_fastq" (to save
    having to traverse the output FASTQ files) and return list of
    Series, one per row, with NumReads set to the number of reads in
    the files.

    :param df: DataFrame with SampleName, Program, File, Description
    columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] ==
                    "riboviz.tools.demultiplex_fastq"]
    if program_df.empty:
        return []
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(
            tuple(fastq.FASTQ_EXTS))]
    count_rows = []
    # TODO edit and cut down when happy with Alternate approach.
    for _, row in files_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        row[NUM_READS] = fastq.count_sequences(file_name)
        count_rows.append(row)
    # Alternate approach: get counts from demultiplex_fastq summary
    # file.
    deplex_files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(
            demultiplex_fastq.NUM_READS_FILE)]
    # demultiplex_fastq outputs one TSV file.
    deplex_file_row = deplex_files_df.iloc[0]
    deplex_file_name = deplex_file_row[workflow_files_logger.FILE]
    deplex_df = pd.read_csv(deplex_file_name, delimiter="\t", comment="#")
    for _, row in files_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        # Extract sample ID from file name.
        tag = os.path.basename(file_name).split(".")[0]
        # Look up sample ID in summary file data.
        tag_row = deplex_df[
            deplex_df[sample_sheets.SAMPLE_ID] == tag].iloc[0]
        row[NUM_READS] = tag_row[sample_sheets.NUM_READS]
        count_rows.append(row)
    return count_rows


def count_reads_hisat2(df):
    """
    Extract "hisat2" rows from DataFrame, count number of reads in
    each FASTQ and SAM file output, and return list of Series, one per
    row, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File, Description
    columns
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
        # TODO decide whether to record sequences or primary_sequences
        sequences, primary_sequences = sam_bam.count_sequences(file_name)
        row[NUM_READS] = sequences
        count_rows.append(row)
        primary_row = row.copy()
        primary_row[NUM_READS] = primary_sequences
        count_rows.append(primary_row)
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(
            tuple(fastq.FASTQ_EXTS))]
    for _, row in files_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        row[NUM_READS] = fastq.count_sequences(file_name)
        count_rows.append(row)
    return count_rows


def count_reads_trim_5p_mismatch(df):
    """
    Extract "riboviz.tools.trim_5p_mismatch" rows from DataFrame,
    count number of reads in the SAM file output, as recorded in the
    TSV file output by "trim_5p_mismatch" (to save having to traverse
    the output SAM file) and return list of Series, one per row, with
    NumReads set to the number of reads in the output.

    :param df: DataFrame with SampleName, Program, File, Description
    columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] ==
                    "riboviz.tools.trim_5p_mismatch"]
    if program_df.empty:
        return []
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith("sam")]
    # TODO edit and cut down when happy with Alternate approach.
    # trim_5p_mismatch outputs one SAM file.
    row = program_df.iloc[0]
    file_name = row[workflow_files_logger.FILE]
    sequences, primary_sequences = sam_bam.count_sequences(file_name)
    row[NUM_READS] = sequences
    primary_row = row.copy()
    primary_row[NUM_READS] = primary_sequences
    # Alternate approach: get count from trim_5p_mismatch summary
    # file.
    files_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(
            trim_5p_mismatch.TRIM_5P_MISMATCH_FILE)]
    # trim_5p_mismatch outputs one TSV file.
    trim_file_row = files_df.iloc[0]
    trim_file_name = trim_file_row[workflow_files_logger.FILE]
    trim_df = pd.read_csv(trim_file_name, delimiter="\t", comment="#")
    trim_row = trim_df.iloc[0]
    summary_row = row.copy()
    summary_row[NUM_READS] = trim_row[trim_5p_mismatch.NUM_WRITTEN]
    return [row, primary_row, summary_row]


def count_reads_umi_tools_dedup(df):
    """
    Extract "umi_tools dedup" program_df from DataFrame, count number of
    reads in the BAM file output, and return list of Series, one per
    row, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File, Description
    columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] ==
                    "umi_tools dedup"]
    if program_df.empty:
        return []
    # umi_tools dedup outputs one BAM file.
    row = program_df.iloc[0]
    file_name = row[workflow_files_logger.FILE]
    sequences, primary_sequences = sam_bam.count_sequences(file_name)
    row[NUM_READS] = sequences
    primary_row = row.copy()
    primary_row[NUM_READS] = primary_sequences
    return [row, primary_row]


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
    # Remove Read/Write column with NumReads so rows (Series) from
    # logs_df can be reused when producing the processed DataFrame.
    logs_df = logs_df.drop(columns=workflow_files_logger.READ_WRITE)
    # Stages not corresponding to a specific sample have SAMPLE_NAME
    # == "".
    input_df = logs_df[logs_df[workflow_files_logger.SAMPLE_NAME]
                       == '']
    # Process inputs and outputs that are non-sample-specific.
    total_rows = []
    for count_fn in [count_reads_inputs, count_reads_cutadapt]:
        rows = count_fn(input_df)
        total_rows.extend(rows)
    # Process outputs from demultiplex_fastq - a special case since it
    # takes in non-sample-specific inputs and produces outputs that
    # are both sample-specific and non-sample specific.
    rows = count_reads_demultiplex_fastq(logs_df)
    total_rows.extend(rows)
    # Process sample-specific inputs and outputs.
    samples = list(logs_df[workflow_files_logger.SAMPLE_NAME].unique())
    samples.remove("")
    for sample in samples:
        sample_df = logs_df[logs_df[workflow_files_logger.SAMPLE_NAME]
                            == sample]
        for count_fn in [count_reads_inputs,
                         count_reads_cutadapt,
                         count_reads_hisat2,
                         count_reads_trim_5p_mismatch,
                         count_reads_umi_tools_dedup]:
            rows = count_fn(sample_df)
            total_rows.extend(rows)
    if total_rows:
        reads_df = reads_df.append(total_rows, ignore_index=True)
    return reads_df


SAMPLE_COUNT_FNS = [count_reads_cutadapt,
                    count_reads_hisat2,
                    count_reads_trim_5p_mismatch,
                    count_reads_umi_tools_dedup]


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
