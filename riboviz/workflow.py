#!/usr/bin/env python
"""
RiboViz workflow-related constants, types and functions.
"""

import collections
from datetime import datetime
import errno
import logging
import os
import os.path
from riboviz import process_utils
from riboviz import logging_utils
from riboviz.utils import value_in_dict


CmdConfigTuple = collections.namedtuple(
    "CmdConfigTuple", ["cmd_file", "is_dry_run"])
""" Command-related configuration """
OutputConfigTuple = collections.namedtuple(
    "OutputConfigTuple",
    ["index_dir", "tmp_dir", "out_dir", "logs_dir"])
""" Output directories """


logging_utils.configure_logging()
LOGGER = logging.getLogger(__name__)
""" Logger """


def setup_output_directories(config, cmd_config):
    """
    Parse configuration for dir_index, dir_tmp, dir_out, dir_logs,
    create paths, and, if is_dry_run is True, make all output
    directories if they do not exist.

    :param config: RiboViz configuration
    :type config: dict
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :return: Output directories
    :rtype: OutputConfigTuple
    :raise KeyError: if a configuration parameter is mssing
    :raise AssertionError: if there is a problem configuring a
    directory
    """
    index_dir = config["dir_index"]
    tmp_dir = config["dir_tmp"]
    out_dir = config["dir_out"]
    base_logs_dir = config["dir_logs"]
    logs_dir = os.path.join(
        base_logs_dir, datetime.now().strftime('%Y%m%d-%H%M%S'))
    dirs = [index_dir, tmp_dir, out_dir, base_logs_dir, logs_dir]
    for d in dirs:
        with open(cmd_config.cmd_file, "a") as f:
            f.write("mkdir -p %s\n" % d)
    if not cmd_config.is_dry_run:
        for d in dirs:
            if not os.path.exists(d):
                os.makedirs(d)
    return OutputConfigTuple(index_dir, tmp_dir, out_dir, logs_dir)


def build_indices(fasta, ht_prefix, log_file, cmd_config):
    """
    Build indices for alignment via invocation of hisat2-build.
    Index files have name <ht_prefix>.<N>.ht2.

    :param fasta: FASTA file to be indexed
    :type fasta: str or unicode
    :param ht_prefix: Prefix of HT2 index files
    :type ht_prefix: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if fasta or hisat2-build cannot be found
    :raise AssertionError: if hisat2-build returns non-zero exit
    code
    """
    LOGGER.info("Create indices (%s). Log: %s", fasta, log_file)
    if not os.path.exists(fasta):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), fasta)
    cmd = ["hisat2-build", fasta, ht_prefix]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def cut_adapters(adapter, original_fq, trimmed_fq, log_file, cmd_config):
    """
    Trim adapters using cutadapt.

    :param adapter: Adapter to trim
    :type adapter: str or unicode
    :param original_fq: FASTQ file to have adapters trimmed
    :type original_fq: str or unicode
    :param trimmed_fq: FASTQ file with trimmed adapters
    :type trimmed_fq: str or unicode
    :param log_file Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if cutadapt cannot be found
    :raise AssertionError: if cutadapt returns non-zero exit
    code
    """
    LOGGER.info("Cut Illumina adapters. Log: %s", log_file)
    cmd = ["cutadapt", "--trim-n", "-O", "1", "-m", "5",
           "-a", adapter, "-o", trimmed_fq, original_fq]
    # cutadapt allows all available processors to be requested.
    cmd += ["-j", str(0)]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def extract_barcodes_umis(
        original_fq, extract_fq, regexp, log_file, cmd_config):
    """
    Extract barcodes and UMIs using "umi_tools extract".

    :param original_fq: FASTQ file to have UMIs and barcodes extracted
    :type original_fq: str or unicode
    :param extract_fq: FASTQ file with UMIs and barcodes extracted
    :type extract_fq: str or unicode
    :param regexp: umi_tools-compliant regular expression to extract
    barcodes and UMIs
    :type regexp: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
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
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run,
        cmd_to_log)


def map_to_r_rna(fastq, index, mapped_sam, unmapped_fastq, nprocesses,
                 log_file, cmd_config):
    """
    Map reads to rRNA.

    :param fastq: FASTQ file
    :type fastq: str or unicode
    :param index: rRNA index file for alignment
    :type index: str or unicode
    :param mapped_sam: SAM file for mapped reads
    :type mapped_sam: str or unicode
    :param unmapped_fastq: FASTQ file for unmapped reads
    :type unmapped_fastq: str or unicode
    :param nprocesses: Number of processes
    :type nprocesses: int
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if fasta or hisat2-build cannot be found
    :raise AssertionError: if hisat2-build returns non-zero exit
    code
    """
    LOGGER.info("Map reads to rRNA. Log: %s", log_file)
    cmd = ["hisat2", "-p", str(nprocesses), "-N", "1",
           "--un", unmapped_fastq, "-x", index,
           "-S", mapped_sam, "-U", fastq]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def map_to_orf(fastq, index, mapped_sam, unmapped_fastq, nprocesses,
               log_file, cmd_config):
    """
    Map reads to ORF.

    :param fastq: FASTQ file
    :type fastq: str or unicode
    :param index: ORF index file for alignment
    :type index: str or unicode
    :param mapped_sam: SAM file for mapped reads
    :type mapped_sam: str or unicode
    :param unmapped_fastq: FASTQ file for unmapped reads
    :type unmapped_fastq: str or unicode
    :param nprocesses: Number of processes
    :type nprocesses: int
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if fasta or hisat2-build cannot be found
    :raise AssertionError: if hisat2-build returns non-zero exit
    code
    """
    LOGGER.info(
        "Map to ORFs with up to 2 alignments. Log: %s", log_file)
    cmd = ["hisat2", "-p", str(nprocesses), "-k", "2",
           "--no-spliced-alignment", "--rna-strandness",
           "F", "--no-unal", "--un", unmapped_fastq,
           "-x", index, "-S", mapped_sam,
           "-U", fastq]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def trim_5p_mismatches(orf_map_sam, orf_map_sam_clean, py_scripts,
                       log_file, cmd_config):
    """
    Trim 5' mismatches.

    :param orf_map_sam: ORF-mapped reads
    :type orf_map_sam: str or unicode
    :param orf_map_sam_clean: Cleaned ORF-mapped reads
    :type orf_map_sam_clean: str or unicode
    :param py_scripts: Python scripts directory
    :type py_scripts: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if python cannot be found
    :raise AssertionError: if python returns non-zero exit code
    """
    LOGGER.info(
        "Trim 5' mismatched nt and remove reads with >1 mismatch. Log: %s",
        log_file)
    cmd = ["python", os.path.join(py_scripts, "trim_5p_mismatch.py"),
           "-mm", "2", "-in", orf_map_sam, "-out", orf_map_sam_clean]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def sort_bam(sam_file, bam_file, nprocesses, log_file, cmd_config):
    """
    Convert SAM to BAM and sort on genome.

    :param sam_file: SAM file
    :type sam_file: str or unicode
    :param bam_file: BAM file
    :type bam_file: str or unicode
    :param nprocesses: Number of processes
    :type nprocesses: int
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if samtools cannot be found
    :raise AssertionError: if samtools returns non-zero exit code
    """
    LOGGER.info(
        "Convert SAM to BAM and sort on genome. Log: %s", log_file)
    cmd_view = ["samtools", "view", "-b", sam_file]
    cmd_sort = ["samtools", "sort", "-@", str(nprocesses),
                "-O", "bam", "-o", bam_file, "-"]
    process_utils.run_logged_pipe_command(
        cmd_view, cmd_sort, log_file, cmd_config.cmd_file,
        cmd_config.is_dry_run)


def index_bam(bam_file, log_file, cmd_config):
    """
    Index BAM file.

    :param bam_file: BAM file
    :type bam_file: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if samtools cannot be found
    :raise AssertionError: if samtools returns non-zero exit code
    """
    LOGGER.info("Index BAM file. Log: %s", log_file)
    cmd = ["samtools", "index", bam_file]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def group_umis(bam_file, groups_file, log_file, cmd_config):
    """
    Run "umi_tools group" on a BAM file.

    :param bam_file: BAM file, input
    :type bam_file: str or unicode
    :param groups_file: Groups file, output
    :type groups_file: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise OSError: if a third-party tool cannot be found
    :raise FileNotFoundError: if umi_tools cannot be found
    :raise AssertionError: if umi_tools returns non-zero exit code
    """
    LOGGER.info("Identify UMI groups. Log: %s", log_file)
    cmd = ["umi_tools", "group", "-I", bam_file,
           "--group-out", groups_file]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def deduplicate_umis(
        bam_file, dedup_bam_file, stats_prefix, log_file, cmd_config):
    """
    Run "umi_tools dedup" on a BAM file.

    :param bam_file: BAM file, input
    :type bam_file: str or unicode
    :param dedup_bam_file: Deduplicated BAM file, output
    :type dedup_bam_file: str or unicode
    :param stats_prefix: File prefix for deduplication statistics
    :type stats_prefix: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise OSError: if a third-party tool cannot be found
    :raise FileNotFoundError: if umi_tools cannot be found
    :raise AssertionError: if umi_tools returns non-zero exit code
    """
    LOGGER.info("Deduplicate UMIs. Log: %s", log_file)
    cmd = ["umi_tools", "dedup", "-I", bam_file, "-S", dedup_bam_file,
           "--output-stats=" + stats_prefix]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def make_bedgraph(bam_file, bedgraph_file, is_plus, log_file, cmd_config):
    """
    Calculate transcriptome coverage and save as a bedgraph.

    :param bam_file: BAM file, input
    :type bam_file: str or unicode
    :param bedgraph_file: Bedgraph file, output
    :type bedgraph_file: str or unicode
    :param is_plus: Is bedgraph to be created for plus strand (True)
    or minus strand (False)
    :type is_plus: bool
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
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
    cmd = ["bedtools", "genomecov", "-ibam", bam_file,
           "-trackline", "-bga", "-5", "-strand", strand]
    process_utils.run_logged_redirect_command(
        cmd, bedgraph_file, log_file, cmd_config.cmd_file,
        cmd_config.is_dry_run)


def bam_to_h5(bam_file, h5_file, orf_gff_file, config, nprocesses,
              r_scripts, log_file, cmd_config):
    """
    Make length-sensitive alignments in H5 format

    :param bam_file: BAM file, input
    :type bam_file: str or unicode
    :param h5_file: H5 file, output
    :type h5_file: str or unicode
    :param orf_gff_file: GFF2/GFF3 file for ORFs
    :type orf_gff_file: str or unicode
    :param config: RiboViz configuration
    :type config: dict
    :param nprocesses: Number of processes
    :type nprocesses: int
    :param r_scripts: R scripts directory
    :type r_scripts: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise KeyError: if a configuration parameter is mssing
    :raise FileNotFoundError: if Rscript cannot be found
    :raise AssertionError: if Rscript returns non-zero exit code
    """
    LOGGER.info("Make length-sensitive alignments in H5 format. Log: %s",
                log_file)
    second_id = config["SecondID"]
    if second_id is None:
        second_id = "NULL"
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "bam_to_h5.R"),
           "--Ncores=" + str(nprocesses),
           "--MinReadLen=" + str(config["MinReadLen"]),
           "--MaxReadLen=" + str(config["MaxReadLen"]),
           "--Buffer=" + str(config["Buffer"]),
           "--PrimaryID=" + config["PrimaryID"],
           "--SecondID=" + second_id,
           "--dataset=" + config["dataset"],
           "--bamFile=" + bam_file,
           "--hdFile=" + h5_file,
           "--orf_gff_file=" + orf_gff_file,
           "--ribovizGFF=" + str(config["ribovizGFF"]),
           "--StopInCDS=" + str(config["StopInCDS"])]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def generate_stats_figs(h5_file, out_dir, prefix, config, nprocesses,
                        r_scripts, log_file, cmd_config):
    """
    Create summary statistics and analyses plots.

    :param h5_file: H5 file
    :type h5_file: str or unicode
    :param out_dir: Directory for output files
    :type out_dir: str or unicode
    :param prefix: Output file name prefix
    :type prefix: str or unicode
    :param config: RiboViz configuration
    :type config: dict
    :param nprocesses: Number of processes
    :type nprocesses: int
    :param r_scripts: R scripts directory
    :type r_scripts: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise KeyError: if a configuration parameter is mssing
    :raise FileNotFoundError: if Rscript cannot be found or any files
    specified in t_rna, codon_pos, features_file configuration
    parameters cannot be found
    :raise AssertionError: if Rscript returns non-zero exit code
    """
    LOGGER.info("Create summary statistics and analyses plots. Log: %s",
                log_file)
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "generate_stats_figs.R"),
           "--Ncores=" + str(nprocesses),
           "--MinReadLen=" + str(config["MinReadLen"]),
           "--MaxReadLen=" + str(config["MaxReadLen"]),
           "--Buffer=" + str(config["Buffer"]),
           "--PrimaryID=" + config["PrimaryID"],
           "--dataset=" + config["dataset"],
           "--hdFile=" + h5_file,
           "--out_prefix=" + prefix,
           "--orf_fasta=" + config["orf_fasta"],
           "--rpf=" + str(config["rpf"]),
           "--dir_out=" + out_dir,
           "--do_pos_sp_nt_freq=" + str(config["do_pos_sp_nt_freq"])]
    for flag in ["t_rna", "codon_pos", "features_file"]:
        if value_in_dict(flag, config):
            flag_file = config[flag]
            if not os.path.exists(flag_file):
                raise FileNotFoundError(errno.ENOENT,
                                        os.strerror(errno.ENOENT),
                                        flag_file)
            cmd.append("--" + flag + "=" + flag_file)
    for flag in ["orf_gff_file", "count_threshold"]:
        if value_in_dict(flag, config):
            cmd.append("--" + flag + "=" + str(config[flag]))
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def collate_tpms(out_dir, samples, r_scripts, log_file, cmd_config):
    """
    Collate TPMs across sample results.

    :param out_dir Output files directory
    :type out_dir: str or unicode
    :param samples: Sample names
    :type samples: list(str or unicode)
    :param r_scripts: R scripts directory
    :type r_scripts: str or unicode
    :param log_file: Log file
    :type log_file: str or unicode
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if Rscript cannot be found
    :raise AssertionError: if Rscript returns non-zero exit code
    """
    LOGGER.info("Collate TPMs across all processed samples. Log: %s",
                log_file)
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "collate_tpms.R"),
           "--dir_out=" + out_dir]
    cmd += samples
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def demultiplex_fastq(fastq, barcodes_file, deplex_dir, log_file,
                      cmd_config):
    """
    Demultiplex reads.

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
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :raise FileNotFoundError: if python cannot be found
    :raise AssertionError: if python returns non-zero exit code
    """
    LOGGER.info("Demultiplex reads. Log: %s", log_file)
    cmd = ["python", "-m", "riboviz.tools.demultiplex_fastq",
           "-r1", fastq, "-ss", barcodes_file, "-o", deplex_dir,
           "-m", "2"]
    process_utils.run_logged_command(
        cmd, log_file, cmd_config.cmd_file, cmd_config.is_dry_run)


def get_sample_log_file(logs_dir, sample, step, index):
    """
    Get name of log file for a specific processing step applied to a
    specific sample.

    :param logs_dir Log files directory
    :type logs_dir: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param step: Name of processing step
    :type step: str or unicode
    :param index: Index of processing step, 1..N
    :type index: int
    :return: file name
    :rtype: str or unicode
    """
    return os.path.join(logs_dir, "%s_%02d_%s.log" % (sample, index, step))
