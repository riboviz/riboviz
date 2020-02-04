#!/usr/bin/env python
"""
RiboViz workflow-related constants, types and functions.
"""

import collections
import glob
import logging
import os
import os.path
from riboviz import params
from riboviz import process_utils
from riboviz import logging_utils
from riboviz import sam_bam
from riboviz import sample_sheets
from riboviz import workflow_r
from riboviz import workflow_files_logger
from riboviz import demultiplex_fastq as demultiplex_fastq_module
from riboviz.tools import demultiplex_fastq as demultiplex_fastq_tools_module
from riboviz.tools import trim_5p_mismatch as trim_5p_mismatch_tools_module
from riboviz.utils import value_in_dict


RunConfigTuple = collections.namedtuple(
    "RunConfigTuple", ["r_scripts",
                       "cmd_file",
                       "workflow_files_log_file",
                       "is_dry_run",
                       "logs_dir",
                       "nprocesses"])
""" Run-related configuration """

DEFAULT_CMD_FILE = "run_riboviz_vignette.sh"
""" Default command file """

ADAPTER_TRIM_FQ = "trim.fq"
""" Adapter trimmed reads """
NON_RRNA_FQ = "nonrRNA.fq"
""" Non-rRNA reads """
RRNA_MAP_SAM = "rRNA_map.sam"
""" rRNA-mapped reads """
ORF_MAP_SAM = "orf_map.sam"
""" ORF-mapped reads """
ORF_MAP_CLEAN_SAM = "orf_map_clean.sam"
""" ORF-mapped reads with mismatched nts trimmed """
UNALIGNED_FQ = "unaligned.fq"
""" Unaligned reads """
UMI_EXTRACT_FQ = "extract_trim.fq"
""" Adapter trimmed reads with UMIs extracted """
PRE_DEDUP_BAM = "pre_dedup.bam"
""" BAM file prior to deduplication """
PRE_DEDUP_GROUPS_TSV = "pre_dedup_groups.tsv"
""" UMI groups before deduplication """
POST_DEDUP_GROUPS_TSV = "post_dedup_groups.tsv"
""" UMI groups after deduplication """
DEDUP_STATS_PREFIX = "dedup_stats"
""" UMI deduplication statistics file prefix """
DEPLEX_DIR_FORMAT = "{}_deplex"
""" Demultiplexed data directory """
ADAPTER_TRIM_FQ_FORMAT = "{}_trim.fq"
""" Adapter trimmed multiplexed reads"""
UMI_EXTRACT_FQ_FORMAT = "{}_extract_trim.fq"
""" Adapter trimmed multiplexed reads with UMIs extracted """
BAM_FORMAT = "{}.bam"
""" Reads mapped to transcripts """
BAM_BAI_FORMAT = "{}.bam.bai"
""" Reads mapped to transcripts (index) """
H5_FORMAT = "{}.h5"
""" Length-sensitive alignments in compressed h5 format """
MINUS_BEDGRAPH = "minus.bedgraph"
""" bedgraph of reads from minus strand """
PLUS_BEDGRAPH = "plus.bedgraph"
""" bedgraph of reads from plus strand """
WORKFLOW_FILES_LOG_FILE = "workflow_files.tsv"
""" Workflow files log file """


logging_utils.configure_logging()
LOGGER = logging.getLogger(__name__)
""" Logger """


def create_directory(directory, cmd_file, is_dry_run=False):
    """
    Add bash command to create `directory` to `cmd_file` and, if
    `is_dry_run` is `True` create the directory.
     a directory and, optionally,

    :param directory: Directory
    :type directory: str or unicode
    :param cmd_file: Commands file
    :type cmd_file: str or unicode
    :param is_dry_run: Don't execute workflow commands (useful for
    seeing what commands would be executed)
    :type is_dry_run: bool
    """
    with open(cmd_file, "a") as f:
        f.write("mkdir -p %s\n" % directory)
    if not is_dry_run:
        if not os.path.exists(directory):
            os.makedirs(directory)


def build_indices(fasta, index_dir, ht_prefix, log_file, run_config):
    """
    Build indices for alignment using hisat2-build.
    Index files have name <ht_prefix>.<N>.ht2.

    :param fasta: FASTA file to be indexed
    :type fasta: str or unicode
    :param index_dir: Index directory
    :type index_dir: str or unicode
    :param ht_prefix: Prefix of HT2 index files
    :type ht_prefix: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if hisat2-build cannot be found
    :raise AssertionError: if hisat2-build returns non-zero exit
    code
    """
    LOGGER.info("Build indices for alignment (%s). Log: %s", fasta, log_file)
    # Explicitly invoke --version as hisat2-build does not log its own
    # version during a processing run.
    cmd = ["hisat2-build", "--version"]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    index_file_path = os.path.join(index_dir, ht_prefix)
    cmd = ["hisat2-build", fasta, index_file_path]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        index_files = glob.glob(index_file_path + "*")
        index_files.sort()
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "hisat2-build",
            [fasta],
            index_files)


def cut_adapters(sample_id, adapter, original_fq, trimmed_fq,
                 log_file, run_config):
    """
    Cut out sequencing library adapters using cutadapt.

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param adapter: Adapter to trim
    :type adapter: str or unicode
    :param original_fq: FASTQ file to have adapters trimmed
    :type original_fq: str or unicode
    :param trimmed_fq: FASTQ file with trimmed adapters
    :type trimmed_fq: str or unicode
    :param log_file Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if cutadapt cannot be found
    :raise AssertionError: if cutadapt returns non-zero exit
    code
    """
    LOGGER.info("Cut out sequencing library adapters. Log: %s", log_file)
    cmd = ["cutadapt", "--trim-n", "-O", "1", "-m", "5",
           "-a", adapter, "-o", trimmed_fq, original_fq]
    cmd += ["-j", str(0)]  # Request all available processors
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "cutadapt",
            [original_fq],
            [trimmed_fq],
            sample_id)


def extract_barcodes_umis(sample_id, original_fq, extract_fq, regexp,
                          log_file, run_config):
    """
    Extract barcodes and UMIs using "umi_tools extract".

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param original_fq: FASTQ file to have UMIs and barcodes extracted
    :type original_fq: str or unicode
    :param extract_fq: FASTQ file with UMIs and barcodes extracted
    :type extract_fq: str or unicode
    :param regexp: umi_tools-compliant regular expression to extract
    barcodes and UMIs
    :type regexp: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if umi_tools cannot be found
    :raise AssertionError: if umi_tooles returns non-zero exit
    code
    """
    LOGGER.info("Extract barcodes and UMIs. Log: %s", log_file)
    # Command to be run from within Python differs subtly from the
    # command logged for the following reason:
    # Assume config["umi_regexp"] is, for example
    #     ^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$
    # If cmd is configured to include:
    #     "--bc-pattern=" + "^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$"
    # then "umi_tools extract" is invoked successfuly from within
    # Python. However the command logged (for rerunning as a bash
    # script) includes:
    #     --bc-pattern=^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$
    # which does not run in bash. It fails with, for
    # example:
    #     syntax error near unexpected token ('
    # If cmd is configured to include:
    #     "--bc-pattern=" + "\"^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$\""
    # then "umi_tools extract" fails as the escaped quotes are
    # treated as part of the regular expression. However the
    # command logged (for rerunning as a bash script) includes:
    #     --bc-pattern="^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$"
    # which does run in bash.
    # Hence, the command run from within Python differs
    # from that logged by the abscence of the quotes round the
    # regular expression.
    pattern_parameter = "--bc-pattern=" + regexp
    cmd = ["umi_tools", "extract", "-I", original_fq, pattern_parameter,
           "--extract-method=regex", "-S", extract_fq]
    cmd_to_log = ["--bc-pattern=" + "\"" + regexp + "\""
                  if c == pattern_parameter else c for c in cmd]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run,
                                     cmd_to_log)
    if not run_config.is_dry_run:
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "umitools extract",
            [original_fq],
            [extract_fq],
            sample_id)


def map_to_r_rna(sample_id, fastq, index_dir, ht_prefix, mapped_sam,
                 unmapped_fastq, log_file, run_config):
    """
    Remove rRNA or other contaminating reads by alignment to rRNA
    index files using hisat2.

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param fastq: FASTQ file
    :type fastq: str or unicode
    :param index_dir: Index directory
    :type index_dir: str or unicode
    :param ht_prefix: Prefix of HT2 rRNA index files for alignment
    :type ht_prefix: str or unicode
    :param mapped_sam: SAM file for mapped reads
    :type mapped_sam: str or unicode
    :param unmapped_fastq: FASTQ file for unmapped reads
    :type unmapped_fastq: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if hisat2-build cannot be found
    :raise AssertionError: if hisat2-build returns non-zero exit
    code
    """
    LOGGER.info(
        "Remove rRNA or other contaminating reads by alignment to rRNA index files. Log: %s",
        log_file)
    # Explicitly invoke --version as hisat2 does not log its own
    # version during a processing run.
    cmd = ["hisat2", "--version"]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    index_file_path = os.path.join(index_dir, ht_prefix)
    cmd = ["hisat2", "-p", str(run_config.nprocesses), "-N", "1",
           "--un", unmapped_fastq, "-x", index_file_path,
           "-S", mapped_sam, "-U", fastq]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        index_files = glob.glob(index_file_path + "*")
        index_files.sort()
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "hisat2",
            [fastq] + index_files,
            [unmapped_fastq, mapped_sam],
            sample_id)


def map_to_orf(sample_id, fastq, index_dir, ht_prefix, mapped_sam,
               unmapped_fastq, log_file, run_config):
    """
    Align remaining reads to ORFs index files using hisat2.

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param fastq: FASTQ file
    :type fastq: str or unicode
    :param index_dir: Index directory
    :type index_dir: str or unicode
    :param ht_prefix: Prefix of HT2 ORF index files for alignment
    :type ht_prefix: str or unicode
    :param mapped_sam: SAM file for mapped reads
    :type mapped_sam: str or unicode
    :param unmapped_fastq: FASTQ file for unmapped reads
    :type unmapped_fastq: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if hisat2-build cannot be found
    :raise AssertionError: if hisat2-build returns non-zero exit
    code
    """
    LOGGER.info(
        "Align remaining reads to ORFs index files using hisat2. Log: %s",
        log_file)
    # Explicitly invoke --version as hisat2 does not log its own
    # version during a processing run.
    cmd = ["hisat2", "--version"]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    index_file_path = os.path.join(index_dir, ht_prefix)
    cmd = ["hisat2", "-p", str(run_config.nprocesses), "-k", "2",
           "--no-spliced-alignment", "--rna-strandness",
           "F", "--no-unal", "--un", unmapped_fastq,
           "-x", index_file_path, "-S", mapped_sam,
           "-U", fastq]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        index_files = glob.glob(index_file_path + "*")
        index_files.sort()
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "hisat2",
            [fastq] + index_files,
            [unmapped_fastq, mapped_sam],
            sample_id)


def trim_5p_mismatches(sample_id, orf_map_sam, orf_map_sam_clean,
                       summary_file, log_file, run_config):
    """
    Trim 5' mismatches from reads and remove reads with more than 2
    mismatches using trim_5p_mismatches.py

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param orf_map_sam: ORF-mapped reads
    :type orf_map_sam: str or unicode
    :param orf_map_sam_clean: Cleaned ORF-mapped reads
    :type orf_map_sam_clean: str or unicode
    :param summary_file: TSV summary file with "num_processed",
    "num_discarded", "num_trimmed" and "num_written" columns
    :type summary_file: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if python cannot be found
    :raise AssertionError: if python returns non-zero exit code
    """
    LOGGER.info(
        "Trim 5' mismatches from reads and remove reads with more than 2 mismatches. Log: %s",
        log_file)
    cmd = ["python", "-m", trim_5p_mismatch_tools_module.__name__,
           "-m", "2", "-i", orf_map_sam, "-o", orf_map_sam_clean,
           "-s", summary_file]
    process_utils.run_logged_command(
        cmd, log_file, run_config.cmd_file, run_config.is_dry_run)
    if not run_config.is_dry_run:
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            trim_5p_mismatch_tools_module.__name__,
            [orf_map_sam],
            [orf_map_sam_clean, summary_file],
            sample_id)


def sort_bam(sample_id, sam_file, bam_file, log_file, run_config):
    """
    Convert SAM to BAM and sort on genome using "samtools view" and
    "samtools sort".

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param sam_file: SAM file
    :type sam_file: str or unicode
    :param bam_file: BAM file
    :type bam_file: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if samtools cannot be found
    :raise AssertionError: if samtools returns non-zero exit code
    """
    LOGGER.info(
        "Convert SAM to BAM and sort on genome. Log: %s", log_file)
    # Explicitly invoke --version as samtools does not log its own
    # version during a processing run.
    cmd = ["samtools", "--version"]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    cmd_view = ["samtools", "view", "-b", sam_file]
    cmd_sort = ["samtools", "sort", "-@", str(run_config.nprocesses),
                "-O", "bam", "-o", bam_file, "-"]
    process_utils.run_logged_pipe_command(cmd_view, cmd_sort,
                                          log_file,
                                          run_config.cmd_file,
                                          run_config.is_dry_run)
    if not run_config.is_dry_run:
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "samtools view | samtools sort",
            [sam_file],
            [bam_file],
            sample_id)


def index_bam(sample_id, bam_file, log_file, run_config):
    """
    Index BAM file using "samtools index".

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param bam_file: BAM file
    :type bam_file: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if samtools cannot be found
    :raise AssertionError: if samtools returns non-zero exit code
    """
    LOGGER.info("Index BAM file. Log: %s", log_file)
    # Explicitly invoke --version as samtools does not log its own
    # version during a processing run.
    cmd = ["samtools", "--version"]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    cmd = ["samtools", "index", bam_file]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "samtools index",
            [bam_file],
            [sam_bam.BAM_BAI_FORMAT.format(bam_file)],
            sample_id)


def group_umis(sample_id, bam_file, groups_file, log_file, run_config):
    """
    Idenfity UMI groups using "umi_tools group".

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param bam_file: BAM file, input
    :type bam_file: str or unicode
    :param groups_file: Groups file, output
    :type groups_file: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise OSError: if a third-party tool cannot be found
    :raise FileNotFoundError: if umi_tools cannot be found
    :raise AssertionError: if umi_tools returns non-zero exit code
    """
    LOGGER.info("Idenfity UMI groups. Log: %s", log_file)
    cmd = ["umi_tools", "group", "-I", bam_file,
           "--group-out", groups_file]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "umi_tools group",
            [bam_file, sam_bam.BAM_BAI_FORMAT.format(bam_file)],
            [groups_file],
            sample_id)


def deduplicate_umis(sample_id, bam_file, dedup_bam_file,
                     stats_prefix, log_file, run_config):
    """
    Deduplicate UMIs using "umi_tools dedup".

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param bam_file: BAM file, input
    :type bam_file: str or unicode
    :param dedup_bam_file: Deduplicated BAM file, output
    :type dedup_bam_file: str or unicode
    :param stats_prefix: File prefix for deduplication statistics
    :type stats_prefix: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise OSError: if a third-party tool cannot be found
    :raise FileNotFoundError: if umi_tools cannot be found
    :raise AssertionError: if umi_tools returns non-zero exit code
    """
    LOGGER.info("Deduplicate UMIs. Log: %s", log_file)
    cmd = ["umi_tools", "dedup", "-I", bam_file, "-S", dedup_bam_file,
           "--output-stats=" + stats_prefix]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        stats_files = glob.glob(stats_prefix + "*")
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "umi_tools dedup",
            [bam_file, sam_bam.BAM_BAI_FORMAT.format(bam_file)],
            [dedup_bam_file] + stats_files,
            sample_id)


def make_bedgraph(sample_id, bam_file, bedgraph_file, is_plus,
                  log_file, run_config):
    """
    Calculate transcriptome coverage and save as a bedgraph using
    "bedtools genomecov".

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param bam_file: BAM file, input
    :type bam_file: str or unicode
    :param bedgraph_file: Bedgraph file, output
    :type bedgraph_file: str or unicode
    :param is_plus: Is bedgraph to be created for plus strand (True)
    or minus strand (False)
    :type is_plus: bool
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise OSError: if a third-party tool cannot be found
    :raise FileNotFoundError: if bedtools cannot be found
    :raise AssertionError: if bedtools returns non-zero exit code
    """
    if is_plus:
        strand = "+"
    else:
        strand = "-"
    LOGGER.info(
        "Calculate transcriptome coverage for %s strand and save as a bedgraph. Log: %s",
        strand, log_file)
    # Explicitly invoke --version as bedtools does not log its own
    # version during a processing run.
    cmd = ["bedtools", "--version"]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    cmd = ["bedtools", "genomecov", "-ibam", bam_file,
           "-trackline", "-bga", "-5", "-strand", strand]
    process_utils.run_logged_redirect_command(cmd, bedgraph_file,
                                              log_file,
                                              run_config.cmd_file,
                                              run_config.is_dry_run)
    if not run_config.is_dry_run:
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            "bedtools",
            [bam_file, sam_bam.BAM_BAI_FORMAT.format(bam_file)],
            [bedgraph_file],
            sample_id)


def bam_to_h5(sample_id, bam_file, h5_file, orf_gff_file, config,
              log_file, run_config):
    """
    Make length-sensitive alignments in H5 format using bam_to_h5.R.

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param bam_file: BAM file, input
    :type bam_file: str or unicode
    :param h5_file: H5 file, output
    :type h5_file: str or unicode
    :param orf_gff_file: GFF2/GFF3 file for ORFs
    :type orf_gff_file: str or unicode
    :param config: RiboViz configuration
    :type config: dict
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise KeyError: if a configuration parameter is mssing
    :raise FileNotFoundError: if Rscript cannot be found
    :raise AssertionError: if Rscript returns non-zero exit code
    """
    LOGGER.info("Make length-sensitive alignments in H5 format. Log: %s",
                log_file)
    secondary_id = config[params.SECONDARY_ID]
    if secondary_id is None:
        secondary_id = "NULL"
    cmd = ["Rscript", "--vanilla",
           os.path.join(run_config.r_scripts,
                        workflow_r.BAM_TO_H5_R),
           "--num-processes=" + str(run_config.nprocesses),
           "--min-read-length=" + str(config[params.MIN_READ_LENGTH]),
           "--max-read-length=" + str(config[params.MAX_READ_LENGTH]),
           "--buffer=" + str(config[params.BUFFER]),
           "--primary-id=" + config[params.PRIMARY_ID],
           "--secondary-id=" + secondary_id,
           "--dataset=" + config[params.DATASET],
           "--bam-file=" + bam_file,
           "--hd-file=" + h5_file,
           "--orf-gff-file=" + orf_gff_file,
           "--is-riboviz-gff=" + str(config[params.IS_RIBOVIZ_GFF]),
           "--stop-in-cds=" + str(config[params.STOP_IN_CDS])]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            workflow_r.BAM_TO_H5_R,
            [bam_file, sam_bam.BAM_BAI_FORMAT.format(bam_file), orf_gff_file],
            [h5_file],
            sample_id)


def generate_stats_figs(sample_id, h5_file, out_dir, config, log_file,
                        run_config):
    """
    Create summary statistics, and analyses and QC plots for both RPF
    and mRNA datasets using generate_stats_figs.R.

    :param sample_id: Sample ID
    :type sample_id: str or unicode
    :param h5_file: H5 file
    :type h5_file: str or unicode
    :param out_dir: Directory for output files
    :type out_dir: str or unicode
    :param config: RiboViz configuration
    :type config: dict
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise KeyError: if a configuration parameter is mssing
    :raise FileNotFoundError: if Rscript cannot be found
    :raise AssertionError: if Rscript returns non-zero exit code
    """
    LOGGER.info(
        "Create summary statistics, and analyses and QC plots for both RPF and mRNA datasets. Log: %s",
        log_file)
    cmd = ["Rscript", "--vanilla",
           os.path.join(run_config.r_scripts,
                        workflow_r.GENERATE_STATS_FIGS_R),
           "--num-processes=" + str(run_config.nprocesses),
           "--min-read-length=" + str(config[params.MIN_READ_LENGTH]),
           "--max-read-length=" + str(config[params.MAX_READ_LENGTH]),
           "--buffer=" + str(config[params.BUFFER]),
           "--primary-id=" + config[params.PRIMARY_ID],
           "--dataset=" + config[params.DATASET],
           "--hd-file=" + h5_file,
           "--orf-fasta-file=" + config[params.ORF_FASTA_FILE],
           "--rpf=" + str(config[params.RPF]),
           "--output-dir=" + out_dir,
           "--do-pos-sp-nt-freq=" + str(config[params.DO_POS_SP_NT_FREQ])]
    # Add optional flags and values.
    flags = zip([params.T_RNA_FILE, params.CODON_POSITIONS_FILE,
                 params.FEATURES_FILE, params.ORF_GFF_FILE,
                 params.ASITE_DISP_LENGTH_FILE],
                ["t-rna-file", "codon-positions-file",
                 "features-file", "orf-gff-file",
                 "asite-disp-length-file"])
    flag_files = []
    for (flag, parameter) in flags:
        if value_in_dict(flag, config):
            flag_file = config[flag]
            cmd.append("--" + parameter + "=" + flag_file)
            flag_files.append(flag_file)
    if value_in_dict(params.COUNT_THRESHOLD, config):
        cmd.append("--count-threshold=" +
                   str(config[params.COUNT_THRESHOLD]))
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        tsv_files = glob.glob(os.path.join(out_dir, "*.tsv"))
        pdf_files = glob.glob(os.path.join(out_dir, "*.pdf"))
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            workflow_r.GENERATE_STATS_FIGS_R,
            [h5_file, config[params.ORF_FASTA_FILE]] + flag_files,
            tsv_files + pdf_files,
            sample_id)


def collate_tpms(out_dir, samples, log_file, run_config, tpms_file=None):
    """
    Collate TPMs across sample results using collate_tpms.R.

    :param out_dir Output files directory
    :type out_dir: str or unicode
    :param samples: Sample names
    :type samples: list(str or unicode)
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :param tpms_file: TPMS file relative to out_dir, if omitted then
    default, chosen by collate_tpms.R, is used
    :type tpms_file: str or unicode
    :raise FileNotFoundError: if Rscript cannot be found
    :raise AssertionError: if Rscript returns non-zero exit code
    """
    LOGGER.info("Collate TPMs across sample results. Log: %s", log_file)
    cmd = ["Rscript", "--vanilla",
           os.path.join(run_config.r_scripts,
                        workflow_r.COLLATE_TPMS_R),
           "--sample-subdirs=" + str(True),
           "--output-dir=" + out_dir]
    if tpms_file is not None:
        cmd.append("--tpms-file=" + tpms_file)
    cmd += samples
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        if tpms_file is None:
            tpms_file = workflow_r.TPMS_COLLATED_TSV
        tpms_file = os.path.join(out_dir, tpms_file)
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            workflow_r.COLLATE_TPMS_R,
            [os.path.join(out_dir, sample, workflow_r.TPMS_TSV) for sample in samples],
            [tpms_file])


def demultiplex_fastq(fastq, barcodes_file, deplex_dir, log_file,
                      run_config):
    """
    Demultiplex reads using demultiplex_fastq.py.

    :param fastq: FASTQ file to demultiplex
    :type fastq: str or unicode
    :param barcodes_file: Sample sheet filename, tab-delimited
    text format with SampleID and TagRead columns, where TagReads are
    the barcodes to use to demultiplex fastq
    :type barcodes_file: str or unicode
    :param deplex_dir: Directory to write demultiplexed files
    :type deplex_dir: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if python cannot be found
    :raise AssertionError: if python returns non-zero exit code
    """
    LOGGER.info("Demultiplex reads. Log: %s", log_file)
    cmd = ["python", "-m", demultiplex_fastq_tools_module.__name__,
           "-1", fastq, "-s", barcodes_file, "-o", deplex_dir,
           "-m", "2"]
    process_utils.run_logged_command(cmd, log_file,
                                     run_config.cmd_file,
                                     run_config.is_dry_run)
    if not run_config.is_dry_run:
        deplex_files = glob.glob(os.path.join(deplex_dir, "*.*"))
        num_reads_file = glob.glob(
            os.path.join(deplex_dir,
                         demultiplex_fastq_module.NUM_READS_FILE))[0]
        unassigned_fq_file = glob.glob(
            os.path.join(deplex_dir,
                         sample_sheets.UNASSIGNED_TAG + ".*"))[0]
        # Remove specific files to leave only fastq files.
        deplex_files.remove(num_reads_file)
        deplex_files.remove(unassigned_fq_file)
        workflow_files_logger.log_files(
            run_config.workflow_files_log_file,
            demultiplex_fastq_tools_module.__name__,
            [fastq, barcodes_file],
            [num_reads_file, unassigned_fq_file])
        deplex_files.sort()
        for file_name in deplex_files:
            # Extract sample ID from file name.
            tag = os.path.basename(file_name).split(".")[0]
            workflow_files_logger.log_files(
                run_config.workflow_files_log_file,
                demultiplex_fastq_tools_module.__name__,
                [],
                [file_name],
                tag)
