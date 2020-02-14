"""
Scan input, temporary and output directories and count the number of
reads (sequences) processed by specific stages of a RiboViz
workflow. The scan is based on the configuration, directory structure
and file patterns used by RiboViz.

The following information is included:

* Input files: number of reads in the FASTQ files used as inputs.
* 'cutadapt': number of reads in the FASTQ file output.
* 'riboviz.tools.demultiplex_fastq': FASTQ files output by
  "demultiplex_fastq", using the information in the associated
  'num_reads.tsv' summary files, or, if these can't be found, the
  FASTQ files themselves.
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
import yaml
import pandas as pd
from riboviz import demultiplex_fastq
from riboviz import fastq
from riboviz import params
from riboviz import provenance
from riboviz import sam_bam
from riboviz import sample_sheets
from riboviz import trim_5p_mismatch
from riboviz import utils
from riboviz import workflow_files
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


def input_fq(config_file, input_dir):
    """
    Extract reads FASTQ input files.

    The configuration file is checked to see if it has an 'fq_files'
    key whose value is a dict from sample names to sample files. If
    present this mapping is used to determine which input files in
    input_dir are have their reads counted and the sample names for
    these files.

    If there is no 'fq_files' key but there is a 'multiplex_fq_files'
    key then the value of this key is assumed to be a list of
    multiplexed input files and their reads are counted.

    If both keys exist then both sets of input files are traversed.

    If neither key exists then no input files are traversed.

    :param config_file: Configuration file
    :type config_file: str or unicode
    :param input_dir: Directory
    :type input_dir: str or unicode
    :return: List of Series with fields SampleName, Program, File, \
    NumReads, Description or []
    :rtype: list(pandas.core.frame.Series)
    """
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    rows = []
    if utils.value_in_dict(params.FQ_FILES, config):
        sample_files = [(sample_name, os.path.join(input_dir, file_name))
                        for sample_name, file_name in
                        list(config[params.FQ_FILES].items())]
    else:
        sample_files = []
    if utils.value_in_dict(params.MULTIPLEX_FQ_FILES, config):
        multiplex_files = [("", os.path.join(input_dir, file_name))
                           for file_name in config[params.MULTIPLEX_FQ_FILES]]
    else:
        multiplex_files = []
    files = sample_files + multiplex_files
    for (sample_name, file_name) in files:
        print(file_name)
        try:
            num_reads = fastq.count_sequences(file_name)
            row = pd.DataFrame(
                [[sample_name, INPUT, file_name, num_reads, INPUT]],
                columns=HEADER)
            rows.append(row)
        except Exception as e:
            print(e)
            continue
    return rows


def cutadapt_fq(tmp_dir, sample=""):
    """
    Extract reads from FASTQ file output by "cutadapt".

    :param tmp_dir: Directory
    :type tmp_dir: str or unicode
    :param sample: Sample name / subdirectory of tmp_dir, or "" if \
    a multiplexed file was processed
    :type sample: str or unicode
    :param fq_file_name: FASTQ file name
    :type fq_file_name: str or unicode
    :param description: Description of this step
    :type description: str or unicode
    :return: Series with fields SampleName, Program, File, NumReads, \
    Description or None
    :rtype: pandas.core.frame.Series
    """
    fq_files = glob.glob(os.path.join(
        tmp_dir, sample, "*" + workflow_files.ADAPTER_TRIM_FQ))
    # If using with FASTQ files then there may be a
    # file with extension "_extract_trim.fq" which also will be
    # caught by the glob above, so remove this file name.
    umi_files = glob.glob(os.path.join(
        tmp_dir, sample, "*" + workflow_files.UMI_EXTRACT_FQ))
    fq_files = [file_name for file_name in fq_files
                if file_name not in umi_files]
    if not fq_files:
        return None
    fq_file = fq_files[0]  # Only 1 match expected.
    print(fq_file)
    try:
        num_reads = fastq.count_sequences(fq_file)
    except Exception as e:
        print(e)
        return None
    description = "Reads after removal of sequencing library adapters"
    row = pd.DataFrame([[sample, "cutadapt", fq_file, num_reads,
                         description]], columns=HEADER)
    return row


def umi_tools_deplex_fq(tmp_dir):
    """
    Extract reads from FASTQ files output by "demultiplex_fastq",
    using the information in the associated 'num_reads.tsv' summary
    files, or, if these can't be found, the FAST files themselves

    :param tmp_dir: Directory
    :type tmp_dir: str or unicode
    :param sample: Sample name / subdirectory of tmp_dir, or "" if \
    a multiplexed file was processed
    :type sample: str or unicode
    :param fq_file_name: FASTQ file name
    :type fq_file_name: str or unicode
    :param description: Description of this step
    :type description: str or unicode
    :return: List of Series with fields SampleName, Program, File, \
    NumReads, Description or []
    :rtype: list(pandas.core.frame.Series)
    """
    deplex_dirs = glob.glob(os.path.join(
        tmp_dir, workflow_files.DEPLEX_DIR_FORMAT.format("*")))
    if not deplex_dirs:
        return []
    description = "Demultiplexed reads"
    rows = []
    for deplex_dir in deplex_dirs:
        fq_files = [glob.glob(os.path.join(deplex_dir, "*" + ext))
                    for ext in fastq.FASTQ_EXTS]
        # Flatten
        fq_files = [f for files in fq_files for f in files]
        if not fq_files:
            continue
        fq_files.sort()
        tsv_files = glob.glob(
            os.path.join(deplex_dir,
                         demultiplex_fastq.NUM_READS_FILE))
        is_tsv_problem = False
        if tsv_files:
            num_reads_file = tsv_files[0]
            print(num_reads_file)
            try:
                deplex_df = pd.read_csv(num_reads_file,
                                        delimiter="\t",
                                        comment="#")
                for fq_file in fq_files:
                    tag = os.path.basename(fq_file).split(".")[0]
                    tag_df = deplex_df[
                        deplex_df[sample_sheets.SAMPLE_ID] == tag]
                    num_reads = tag_df.iloc[0][sample_sheets.NUM_READS]
                    row = pd.DataFrame(
                        [[tag,
                          demultiplex_fastq_tools_module.__name__,
                          fq_file, num_reads, description]],
                        columns=HEADER)
                    rows.append(row)
            except Exception as e:
                print(e)
                is_tsv_problem = True
        if is_tsv_problem or not tsv_files:
            # Traverse FASTQ files directly.
            for fq_file in fq_files:
                print(fq_file)
                tag = os.path.basename(fq_file).split(".")[0]
                try:
                    num_reads = fastq.count_sequences(fq_file)
                except Exception as e:
                    print(e)
                    continue
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
    :return: Series with fields SampleName, Program, File, NumReads, \
    Description or None
    :rtype: pandas.core.frame.Series
    """
    fq_files = glob.glob(os.path.join(tmp_dir, sample, fq_file_name))
    if not fq_files:
        return None
    fq_file = fq_files[0]  # Only 1 match expected
    print(fq_file)
    try:
        num_reads = fastq.count_sequences(fq_file)
    except Exception as e:
        print(e)
        return None
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
    :return: Series with fields SampleName, Program, File, NumReads, \
    Description or None
    :rtype: pandas.core.frame.Series
    """
    sam_files = glob.glob(os.path.join(tmp_dir, sample, sam_file_name))
    if not sam_files:
        return None
    sam_file = sam_files[0]  # Only 1 match expected.
    print(sam_file)
    try:
        sequences, _ = sam_bam.count_sequences(sam_file)
    except Exception as e:
        print(e)
        return None
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
    :return: Series with fields SampleName, Program, File, NumReads, \
    Description or None
    :rtype: pandas.core.frame.Series
    """
    # Look for the SAM file.
    sam_files = glob.glob(os.path.join(
        tmp_dir, sample, workflow_files.ORF_MAP_CLEAN_SAM))
    if not sam_files:
        return None
    sam_file = sam_files[0]  # Only 1 match expected.
    # Look for trim_5p_mismatch.tsv.
    tsv_files = glob.glob(os.path.join(
        tmp_dir, sample, trim_5p_mismatch.TRIM_5P_MISMATCH_FILE))
    is_tsv_problem = False
    if tsv_files:
        tsv_file = tsv_files[0]
        print(tsv_file)
        try:
            trim_data = pd.read_csv(tsv_file, delimiter="\t", comment="#")
            trim_row = trim_data.iloc[0]
            sequences = trim_row[trim_5p_mismatch.NUM_WRITTEN]
        except Exception as e:
            print(e)
            is_tsv_problem = True
    if is_tsv_problem or not tsv_files:
        # Traverse SAM file directly.
        print(sam_file)
        try:
            sequences, _ = sam_bam.count_sequences(sam_file)
        except Exception as e:
            print(e)
            return None
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
    :return: Series with fields SampleName, Program, File, NumReads, \
    Description or None
    :rtype: pandas.core.frame.Series
    """
    # Look for pre_dedup.bam
    files = glob.glob(os.path.join(
        tmp_dir, sample, workflow_files.PRE_DEDUP_BAM))
    if not files:
        # Deduplication was not done.
        return None
    # Look for the BAM file output.
    files = glob.glob(os.path.join(
        output_dir, sample, sam_bam.BAM_FORMAT.format(sample)))
    if not files:
        return None
    file_name = files[0]  # Only 1 match expected.
    print(file_name)
    try:
        sequences, _ = sam_bam.count_sequences(file_name)
    except Exception as e:
        print(e)
        return None
    description = "Deduplicated reads"
    row = pd.DataFrame(
        [[sample, "umi_tools dedup", file_name, sequences, description]],
        columns=HEADER)
    return row


def count_reads_df(config_file, input_dir, tmp_dir, output_dir):
    """
    Scan input, temporary and output directories and count the number
    of reads (sequences) processed by specific stages of a RiboViz
    workflow. The scan is based on the directory structure and file
    patterns used by RiboViz.

    The DataFrame returned has columns: SampleName, Program, File,
    NumReads, Description.

    :param config_file: Configuration file
    :type config_file: str or unicode
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
    rows.extend(input_fq(config_file, input_dir))
    rows.append(cutadapt_fq(tmp_dir))
    rows.extend(umi_tools_deplex_fq(tmp_dir))
    tmp_samples = [f.name for f in os.scandir(tmp_dir) if f.is_dir()]
    tmp_samples.sort()
    for sample in tmp_samples:
        rows.append(cutadapt_fq(tmp_dir, sample))
        rows.append(hisat2_fq(tmp_dir, sample, workflow_files.NON_RRNA_FQ,
                              "rRNA or other contaminating reads removed by alignment to rRNA index files"))
        rows.append(hisat2_sam(tmp_dir, sample, workflow_files.RRNA_MAP_SAM,
                               "Reads with rRNA and other contaminating reads removed by alignment to rRNA index files"))
        rows.append(hisat2_fq(tmp_dir, sample, workflow_files.UNALIGNED_FQ,
                              "Unaligned reads removed by alignment of remaining reads to ORFs index files"))
        rows.append(hisat2_sam(tmp_dir, sample, workflow_files.ORF_MAP_SAM,
                               "Reads aligned to ORFs index files"))
        rows.append(trim_5p_mismatch_sam(tmp_dir, sample))
        rows.append(umi_tools_dedup_bam(tmp_dir, output_dir, sample))
    rows = [row for row in rows if row is not None]
    df = df.append(rows)
    return df


def count_reads(config_file, input_dir, tmp_dir, output_dir, reads_file):
    """
    Scan input, temporary and output directories and count the number
    of reads (sequences) processed by specific stages of a RiboViz
    workflow. The scan is based on the configuration, directory
    structure and file patterns used by RiboViz.

    reads_file is a TSV file with columns: SampleName, Program, File,
    NumReads, Description.

    :param config_file: Configuration file
    :type config_file: str or unicode
    :param input_dir: Input files directory
    :type input_dir: str or unicode
    :param tmp_dir: Temporary files directory
    :type tmp_dir: str or unicode
    :param output_dir: Output files directory
    :type output_dir: str or unicode
    :param reads_file: Reads file output
    :type reads_file: str or unicode
    """
    reads_df = count_reads_df(config_file, input_dir, tmp_dir, output_dir)
    provenance.write_provenance_header(__file__, reads_file)
    reads_df[list(reads_df.columns)].to_csv(
        reads_file, mode='a', sep="\t", index=False)
