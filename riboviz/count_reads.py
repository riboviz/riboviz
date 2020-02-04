"""
Scan input, temporary and output directories and count the number of
reads (sequences) processed by specific stages of a RiboViz
workflow. The scan is based on the directory structure and file
patterns used by RiboViz.

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
import glob
import os
import pandas as pd
from riboviz import demultiplex_fastq
from riboviz import fastq
from riboviz import file_names
from riboviz import provenance
from riboviz import sam_bam
from riboviz import sample_sheets
from riboviz import trim_5p_mismatch
from riboviz.tools import demultiplex_fastq as demultiplex_fastq_tools_module
from riboviz.tools import trim_5p_mismatch as trim_5p_mismatch_tools_module

SAMPLE_NAME = "SampleName"
""" Sample name column name """
PROGRAM = "Program"
""" Program column name """
FILE = "File"
""" File column name """
NUM_READS = "NumReads"
""" NumReads column name """
DESCRIPTION = "Description"
""" Description column name """
HEADER = [SAMPLE_NAME, PROGRAM, FILE, NUM_READS, DESCRIPTION]
""" Number of reads file header """
INPUT = "input"
"""
Special value for PROGRAM field to denote input files that do not
originate from any step in a workflow
"""


def input_fq(input_dir):
    """
    Extract reads FASTQ input files.

    :param input_dir: Directory
    :type input_dir: str or unicode
    :return: List of Series with fields SampleName, Program, File,
    NumReads, Description or []
    :rtype: list(pandas.core.frame.Series)
    """
    # FASTQ files may differ in their extensions, find all.
    fq_files = [glob.glob(os.path.join(input_dir, "*" + ext))
                for ext in fastq.FASTQ_EXTS]
    # Flatten list of lists
    fq_files = [f for files in fq_files for f in files]
    if not fq_files:
        return []
    fq_files.sort()
    rows = []
    for fq_file in fq_files:
        num_reads = fastq.count_sequences(fq_file)
        row = pd.DataFrame([["", INPUT, fq_file, num_reads, INPUT]],
                           columns=HEADER)
        rows.append(row)
    return rows


def cutadapt_fq(tmp_dir, sample=""):
    """
    Extract reads from FASTQ file output by "cutadapt".

    :param tmp_dir: Directory
    :type tmp_dir: str or unicode
    :param sample: Sample name / subdirectory of tmp_dir, or "" if a
    multiplexed file was processed
    :type sample: str or unicode
    :param fq_file_name: FASTQ file name
    :type fq_file_name: str or unicode
    :param description: Description of this step
    :type description: str or unicode
    :return: Series with fields SampleName, Program, File, NumReads,
    Description or None
    :rtype: pandas.core.frame.Series
    """
    fq_files = glob.glob(os.path.join(
        tmp_dir, sample, "*" + file_names.ADAPTER_TRIM_FQ))
    if sample == "":
        # If using with a multiplexed FASTQ file then there may be a
        # file with extension "_extract_trim.fq" which also will be
        # caught by the glob above, so remove this file name.
        umi_files = glob.glob(os.path.join(
            tmp_dir, sample, "*" + file_names.UMI_EXTRACT_FQ))
        fq_files = [file_name for file_name in fq_files
                    if file_name not in umi_files]
    if not fq_files:
        return None
    fq_file = fq_files[0]  # Only 1 match expected.
    num_reads = fastq.count_sequences(fq_file)
    description = "Reads after removal of sequencing library adapters"
    row = pd.DataFrame([[sample, "cutadapt", fq_file, num_reads,
                         description]], columns=HEADER)
    return row


def umi_tools_deplex_fq(tmp_dir):
    """
    Extract reads from FASTQ file output by "cutadapt".

    :param tmp_dir: Directory
    :type tmp_dir: str or unicode
    :param sample: Sample name / subdirectory of tmp_dir, or "" if a
    multiplexed file was processed
    :type sample: str or unicode
    :param fq_file_name: FASTQ file name
    :type fq_file_name: str or unicode
    :param description: Description of this step
    :type description: str or unicode
    :return: List of Series with fields SampleName, Program, File,
    NumReads, Description or []
    :rtype: list(pandas.core.frame.Series)
    """
    deplex_dirs = glob.glob(os.path.join(
        tmp_dir, file_names.DEPLEX_DIR_FORMAT.format("*")))
    if not deplex_dirs:
        return []
    rows = []
    for deplex_dir in deplex_dirs:
        tsv_files = glob.glob(
            os.path.join(deplex_dir,
                         demultiplex_fastq.NUM_READS_FILE))
        if not tsv_files:
            pass  # TODO: Check if none and handle
        num_reads_file = tsv_files[0]
        deplex_df = pd.read_csv(num_reads_file, delimiter="\t", comment="#")
        fq_files = [glob.glob(os.path.join(deplex_dir, "*" + ext))
                    for ext in fastq.FASTQ_EXTS]
        # Flatten
        fq_files = [f for files in fq_files for f in files]
        if not fq_files:
            continue
        fq_files.sort()
        print(fq_files)
        for fq_file in fq_files:
            tag = os.path.basename(fq_file).split(".")[0]
            tag_df = deplex_df[
                deplex_df[sample_sheets.SAMPLE_ID] == tag]
            num_reads = tag_df.iloc[0][sample_sheets.NUM_READS]
            description = "Demultiplexed reads"
            row = pd.DataFrame(
                [[tag,
                  demultiplex_fastq_tools_module.__name__,
                  fq_file, num_reads, description]],
                columns=HEADER)
            rows.append(row)
    return rows


def hisat2_fq(tmp_dir, sample, fq_file_name, description):
    """
    Extract reads from FASTQ file output by "hisat2".

    :param tmp_dir: Directory
    :type tmp_dir: str or unicode
    :param sample: Sample name / subdirectory of tmp_dir
    :type sample: str or unicode
    :param fq_file_name: FASTQ file name
    :type fq_file_name: str or unicode
    :param description: Description of this step
    :type description: str or unicode
    :return: Series with fields SampleName, Program, File, NumReads,
    Description or None
    :rtype: pandas.core.frame.Series
    """
    fq_files = glob.glob(os.path.join(tmp_dir, sample, fq_file_name))
    if not fq_files:
        return None
    fq_file = fq_files[0]  # Only 1 match expected.
    num_reads = fastq.count_sequences(fq_file)
    row = pd.DataFrame([[sample, "hisat2", fq_file, num_reads,
                         description]], columns=HEADER)
    return row


def hisat2_sam(tmp_dir, sample, sam_file_name, description):
    """
    Extract reads from SAM file output by "hisat2".

    :param tmp_dir: Directory
    :type tmp_dir: str or unicode
    :param sample: Sample name / subdirectory of tmp_dir
    :type sample: str or unicode
    :param sam_file_name: SAM file name
    :type sam_file_name: str or unicode
    :param description: Description of this step
    :type description: str or unicode
    :return: Series with fields SampleName, Program, File, NumReads,
    Description or None
    :rtype: pandas.core.frame.Series
    """
    sam_files = glob.glob(os.path.join(tmp_dir, sample, sam_file_name))
    if not sam_files:
        return None
    sam_file = sam_files[0]  # Only 1 match expected.
    sequences, _ = sam_bam.count_sequences(sam_file)
    row = pd.DataFrame([[sample, "hisat2", sam_file, sequences,
                         description]], columns=HEADER)
    return row


def trim_5p_mismatch_sam(tmp_dir, sample):
    """
    Extract reads from SAM file output by "trim_5p_mismatch", using
    the information in the associated 'trim_5p_mismatch.tsv' summary
    file, or, if this can't be found, the SAM file itself.

    :param tmp_dir: Directory
    :type tmp_dir: str or unicode
    :param sample: Sample name / subdirectory of tmp_dir
    :type sample: str or unicode
    :return: Series with fields SampleName, Program, File, NumReads,
    Description or None
    :rtype: pandas.core.frame.Series
    """
    # Look for the SAM file.
    sam_files = glob.glob(os.path.join(
        tmp_dir, sample, file_names.ORF_MAP_CLEAN_SAM))
    if not sam_files:
        return None
    sam_file = sam_files[0]  # Only 1 match expected.
    # Look for trim_5p_mismatch.tsv.
    tsv_files = glob.glob(os.path.join(
        tmp_dir, sample, trim_5p_mismatch.TRIM_5P_MISMATCH_FILE))
    if tsv_files:
        tsv_file = tsv_files[0]
        trim_data = pd.read_csv(tsv_file, delimiter="\t", comment="#")
        trim_row = trim_data.iloc[0]
        sequences = trim_row[trim_5p_mismatch.NUM_WRITTEN]
    else:
        # Traverse SAM file.
        sequences, _ = sam_bam.count_sequences(sam_file)
    description = "Reads after trimming of 5' mismatches and removal of those with more than 2 mismatches"
    row = pd.DataFrame([[sample,
                         trim_5p_mismatch_tools_module.__name__,
                         sam_file, sequences, description]],
                       columns=HEADER)
    return row


def umi_tools_dedup_bam(tmp_dir, output_dir, sample):
    """
    Extract reads from BAM file output by "umi_tools dedup". A check
    is made for a file, '<tmp_dir>/<sample>/pre_dedup.bam', and, if
    this exists, then reads are counted from
    '<out_dir>/<sample>/<sample>.bam'.

    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    :param output_dir: Output directory
    :type output_dir: str or unicode
    :param sample: Sample name / subdirectory of tmp_dir
    :type sample: str or unicode
    :return: Series with fields SampleName, Program, File, NumReads,
    Description or None
    :rtype: pandas.core.frame.Series
    """
    # Look for pre_dedup.bam
    files = glob.glob(os.path.join(
        tmp_dir, sample, file_names.PRE_DEDUP_BAM))
    if not files:
        # Deduplication was not done.
        return None
    # Look for the BAM file output.
    files = glob.glob(os.path.join(
        output_dir, sample, file_names.BAM_FORMAT.format(sample)))
    if not files:
        return None
    file_name = files[0]  # Only 1 match expected.
    sequences, _ = sam_bam.count_sequences(file_name)
    description = "Deduplicated reads"
    row = pd.DataFrame(
        [[sample, "umi_tools dedup", file_name, sequences, description]],
        columns=HEADER)
    return row


def count_reads_df(input_dir, tmp_dir, output_dir):
    """
    Scan input, temporary and output directories and count the number
    of reads (sequences) processed by specific stages of a RiboViz
    workflow. The scan is based on the directory structure and file
    patterns used by RiboViz.

    The DataFrame returned has columns: SampleName, Program, File,
    NumReads, Description.

    :param input_dir: Input files directory
    :type input_dir: str or unicode
    :param tmp_dir: Temporary files directory
    :type tmp_dir: str or unicode
    :param output_dir: Output files directory
    :type output_dir: str or unicode
    :return: DataFrame
    :rtype: pandas.core.frame.DataFrame
    """
    df = pd.DataFrame(columns=HEADER)
    rows = []
    rows.extend(input_fq(input_dir))
    rows.append(cutadapt_fq(tmp_dir))
    rows.extend(umi_tools_deplex_fq(tmp_dir))
    tmp_samples = [f.name for f in os.scandir(tmp_dir) if f.is_dir()]
    tmp_samples.sort()
    for sample in tmp_samples:
        rows.append(cutadapt_fq(tmp_dir, sample))
        rows.append(hisat2_fq(tmp_dir, sample, file_names.NON_RRNA_FQ,
                              "rRNA or other contaminating reads removed by alignment to rRNA index files"))
        rows.append(hisat2_sam(tmp_dir, sample, file_names.RRNA_MAP_SAM,
                               "Reads with rRNA and other contaminating reads removed by alignment to rRNA index files"))
        rows.append(hisat2_fq(tmp_dir, sample, file_names.UNALIGNED_FQ,
                              "Reads aligned to ORFs index files"))
        rows.append(hisat2_sam(tmp_dir, sample, file_names.ORF_MAP_SAM,
                               "Unaligned reads removed by alignment of remaining reads to ORFs index files"))
        rows.append(trim_5p_mismatch_sam(tmp_dir, sample))
        rows.append(umi_tools_dedup_bam(tmp_dir, output_dir, sample))
    rows = [row for row in rows if row is not None]
    df = df.append(rows)
    return df


def count_reads(input_dir, tmp_dir, output_dir, reads_file):
    """
    Scan input, temporary and output directories and count the number
    of reads (sequences) processed by specific stages of a RiboViz
    workflow. The scan is based on the directory structure and file
    patterns used by RiboViz.

    reads_file is a TSV file with columns: SampleName, Program, File,
    NumReads, Description.

    :param input_dir: Input files directory
    :type input_dir: str or unicode
    :param tmp_dir: Temporary files directory
    :type tmp_dir: str or unicode
    :param output_dir: Output files directory
    :type output_dir: str or unicode
    :param reads_file: Reads file output
    :type reads_file: str or unicode
    """
    reads_df = count_reads_df(input_dir, tmp_dir, output_dir)
    provenance.write_provenance_header(__file__, reads_file)
    reads_df[list(reads_df.columns)].to_csv(
        reads_file, mode='a', sep="\t", index=False)
