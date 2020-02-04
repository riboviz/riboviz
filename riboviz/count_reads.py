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

The files logged in the workflow files log file must exist.

The following information is included:

* Input files: number of reads in the FASTQ files used as inputs.
* 'cutadapt': number of reads in the FASTQ file output.
* 'riboviz.tools.demultiplex_fastq': number of reads in the FASTQ
  files output, as recorded in the 'num_reads.tsv' file output.
* 'hisat2': number of reads in the SAM file and FASTQ file output.
* 'riboviz.tools.trim_5p_mismatch': number of reads in the SAM file
  output as recorded in the TSV summary file output.
* 'umi_tools dedup': number of reads in the BAM file output.

The output file is a TSV file with columns:

* 'SampleName'
* 'Program'
* 'File'
* 'NumReads': Number of reads in the file.
* 'Description': Human-readable description of the file contents.
"""
import os
import pandas as pd
from riboviz import demultiplex_fastq
from riboviz import fastq
from riboviz import provenance
from riboviz import sam_bam
from riboviz import sample_sheets
from riboviz import trim_5p_mismatch
from riboviz import workflow
from riboviz import workflow_files_logger


NUM_READS = "NumReads"
""" NumReads column name """
DESCRIPTION = "Description"
""" Description column name """
HISAT2_DESCRIPTIONS = {
    workflow.NON_RRNA_FQ: "rRNA or other contaminating reads removed by alignment to rRNA index files",
    workflow.RRNA_MAP_SAM: "Reads with rRNA and other contaminating reads removed by alignment to rRNA index files",
    workflow.ORF_MAP_SAM: "Reads aligned to ORFs index files",
    workflow.UNALIGNED_FQ: "Unaligned reads removed by alignment of remaining reads to ORFs index files"
}
""" Mapping of HISAT2 file names to escriptions """
HEADER = [workflow_files_logger.SAMPLE_NAME,
          workflow_files_logger.PROGRAM,
          workflow_files_logger.FILE,
          NUM_READS,
          DESCRIPTION]
""" Number of reads file header """


def count_reads_inputs(df):
    """
    Extract "input" rows from DataFrame, count number of reads in
    each FASTQ file used as an input, and return list of Series, one
    per file, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads and Description set
    :rtype: list(pandas.core.frame.Series)
    """
    fq_df = df[
        (df[workflow_files_logger.PROGRAM] == workflow_files_logger.INPUT)
        &
        (df[workflow_files_logger.FILE].str.lower().str.endswith(
            tuple(fastq.FASTQ_EXTS)))]
    count_rows = []
    for _, row in fq_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        row[NUM_READS] = fastq.count_sequences(file_name)
        row[DESCRIPTION] = "Original reads"
        count_rows.append(row)
    return count_rows


def count_reads_cutadapt(df):
    """
    Extract "cutadapt" rows from DataFrame, count number of reads
    in the FASTQ file output, and return list of single Series, for
    this file, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads and Description set
    :rtype: list(pandas.core.frame.Series)
    """
    fq_df = df[
        (df[workflow_files_logger.PROGRAM] == "cutadapt")
        &
        (df[workflow_files_logger.FILE].str.lower().str.endswith(
            tuple(fastq.FASTQ_EXTS)))]
    if fq_df.empty:
        return []
    row = fq_df.iloc[0]
    file_name = row[workflow_files_logger.FILE]
    row[NUM_READS] = fastq.count_sequences(file_name)
    row[DESCRIPTION] = "Reads after removal of sequencing library adapters"
    return [row]


def count_reads_demultiplex_fastq(df):
    """
    Extract "riboviz.tools.demultiplex_fastq" rows from DataFrame,
    count number of reads in each demultiplexed FASTQ file output, as
    recorded in the TSV file output by "demultiplex_fastq" (to save
    having to traverse the output FASTQ files) and return list of
    Series, one per file, with NumReads set to the number of reads in
    the files.

    :param df: DataFrame with SampleName, Program, File columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads and Description set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] ==
                    "riboviz.tools.demultiplex_fastq"]
    fq_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(
            tuple(fastq.FASTQ_EXTS))]
    if fq_df.empty:
        return []
    count_rows = []
    # Get demultiplex_fastq TSV summary file.
    tsv_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(
            demultiplex_fastq.NUM_READS_FILE)]
    tsv_file_name = tsv_df.iloc[0][workflow_files_logger.FILE]
    deplex_df = pd.read_csv(tsv_file_name, delimiter="\t", comment="#")
    # Iterate through FASTQ entries and lookup counts in TSV summary
    # file data.
    for _, row in fq_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        if file_name.lower().endswith(tuple(fastq.FASTQ_EXTS)):
            # Extract sample ID from file name.
            tag = os.path.basename(file_name).split(".")[0]
            # Look up sample ID in summary file data.
            tag_df = deplex_df[
                deplex_df[sample_sheets.SAMPLE_ID] == tag]
            row[NUM_READS] = tag_df.iloc[0][sample_sheets.NUM_READS]
            row[DESCRIPTION] = "Demultiplexed reads"
            count_rows.append(row)
    return count_rows


def count_reads_hisat2(df):
    """
    Extract "hisat2" rows from DataFrame, count number of reads in
    each FASTQ and SAM file output, and return list of Series, one per
    file, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads and Description set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] == "hisat2"]
    count_rows = []
    # Iterate through rows and check extensions in turn so Series
    # output list is consistent in order with rows in input
    # DataFrame.
    for _, row in program_df.iterrows():
        file_name = row[workflow_files_logger.FILE]
        local_file_name = os.path.basename(file_name)
        if local_file_name in HISAT2_DESCRIPTIONS:
            description = HISAT2_DESCRIPTIONS[local_file_name]
        else:
            description = None
        if file_name.lower().endswith("sam"):
            sequences, _ = sam_bam.count_sequences(file_name)
            row[NUM_READS] = sequences
            if description is None:
                descriptin = "Aligned reads"
            row[DESCRIPTION] = description
            count_rows.append(row)
        elif file_name.lower().endswith(tuple(fastq.FASTQ_EXTS)):
            row[NUM_READS] = fastq.count_sequences(file_name)
            if description is None:
                descriptin = "Unaligned reads"
            row[DESCRIPTION] = description
            count_rows.append(row)
    return count_rows


def count_reads_trim_5p_mismatch(df):
    """
    Extract "riboviz.tools.trim_5p_mismatch" rows from DataFrame,
    count number of reads in the SAM file output, as recorded in the
    TSV file output by "trim_5p_mismatch" (to save having to traverse
    the output SAM file) and return list of single Series, for this
    file, with NumReads set to the number of reads in the output.

    :param df: DataFrame with SampleName, Program, File columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads and Description set
    :rtype: list(pandas.core.frame.Series)
    """
    program_df = df[df[workflow_files_logger.PROGRAM] ==
                    "riboviz.tools.trim_5p_mismatch"]
    sam_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith("sam")]
    if sam_df.empty:
        return []
    # trim_5p_mismatch outputs one SAM file.
    row = sam_df.iloc[0]
    tsv_df = program_df[program_df[
        workflow_files_logger.FILE].str.lower().str.endswith(
            trim_5p_mismatch.TRIM_5P_MISMATCH_FILE)]
    # trim_5p_mismatch outputs one TSV file.
    trim_file_name = tsv_df.iloc[0][workflow_files_logger.FILE]
    trim_data_df = pd.read_csv(trim_file_name, delimiter="\t", comment="#")
    trim_row = trim_data_df.iloc[0]
    row[NUM_READS] = trim_row[trim_5p_mismatch.NUM_WRITTEN]
    row[DESCRIPTION] = "Reads after trimming of 5' mismatches and removal of those with more than 2 mismatches"
    return [row]


def count_reads_umi_tools_dedup(df):
    """
    Extract "umi_tools dedup" program_df from DataFrame, count number
    of reads in the BAM file output, and return single Series, for
    this file, with NumReads set to the number of reads in the files.

    :param df: DataFrame with SampleName, Program, File columns
    :type df: pandas.core.frame.DataFrame
    :return: Filtered rows with NumReads and Description set
    :rtype: list(pandas.core.frame.Series)
    """
    bam_df = df[
        (df[workflow_files_logger.PROGRAM] == "umi_tools dedup")
        &
        (df[workflow_files_logger.FILE].str.lower().str.endswith("bam"))]
    if bam_df.empty:
        return []
    # umi_tools dedup outputs one BAM file.
    row = bam_df.iloc[0]
    file_name = row[workflow_files_logger.FILE]
    sequences, _ = sam_bam.count_sequences(file_name)
    row[NUM_READS] = sequences
    row[DESCRIPTION] = "Deduplicated reads"
    return [row]


def count_input_write_reads(df):
    """
    Process workflow files log file data frame and count reads in
    FASTQ, SAM and BAM files that were input or written.

    df is assumed to be a DataFrame with
    riboviz.workflow_files_logger-consistent columns: SampleName,
    Program, File, Read/Write.

    The DataFrame returned has columns: SampleName, Program, File,
    NumReads, Description.

    :param df: DataFrame with riboviz.workflow_files_logger-consistent
    data
    :type df: pandas.core.frame.DataFrame
    :return: DataFrame
    :rtype: pandas.core.frame.DataFrame
    """
    # Filter to get files recorded as "input" or "write" only.
    inputs_writes_df = df[
        (df[workflow_files_logger.PROGRAM] == workflow_files_logger.INPUT)
        |
        (df[workflow_files_logger.READ_WRITE] == workflow_files_logger.WRITE)
    ]
    reads_df = pd.DataFrame(columns=HEADER)
    # Remove Read/Write column with NumReads so rows (Series) from
    # inputs_writes_df can be reused when producing the processed
    # DataFrame.
    inputs_writes_df = inputs_writes_df.drop(
        columns=workflow_files_logger.READ_WRITE)
    # Process entries not corresponding to a specific sample.
    generics_df = inputs_writes_df[
        inputs_writes_df[workflow_files_logger.SAMPLE_NAME] == '']
    total_rows = []
    for count_fn in [count_reads_inputs, count_reads_cutadapt]:
        rows = count_fn(generics_df)
        total_rows.extend(rows)
    # Process outputs from demultiplex_fastq - a special case since it
    # takes in non-sample-specific inputs and produces outputs that
    # are both sample-specific and non-sample specific.
    rows = count_reads_demultiplex_fastq(inputs_writes_df)
    total_rows.extend(rows)
    # Process entries corresponding to specific samples.
    samples = list(inputs_writes_df[
        workflow_files_logger.SAMPLE_NAME].unique())
    samples.remove("")
    for sample in samples:
        sample_df = inputs_writes_df[
            inputs_writes_df[workflow_files_logger.SAMPLE_NAME] == sample]
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


def count_reads(workflow_file, reads_file, delimiter="\t", comment="#"):
    """
    Process workflow files log file and count reads in FASTQ, SAM and
    BAM files that were input or written and save results.

    workflow_file is assumed to be a TSV file with
    riboviz.workflow_files_logger-consistent columns: SampleName,
    Program, File, Read/Write.

    reads_file is a TSV file with columns: SampleName, Program, File,
    NumReads, Description.

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
