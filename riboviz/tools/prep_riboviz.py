#!/usr/bin/env python
"""
RiboViz workflow.

Usage:

    PYTHONPATH=. python riboviz/tools/prep_riboviz.py \
        [--dry-run] \
        <R_SCRIPTS_DIRECTORY>\
        <YAML_CONFIG_FILE>

Example:

    PYTHONPATH=. python riboviz/tools/prepRiboviz.py \
        rscripts/ \
        vignette/vignette_config.yaml

Prepare ribosome profiling data for RiboViz or other analysis:

* Reads configuration information from YAML configuration file.
* Builds hisat2 indices if requested (config["build_indices"] ==
  True) using "hisat2 build" and saves these into an index directory
  (config["dir_index"]).
* Processes all fastq.gz files (config["dir_in"]). For each fastq.gz
  file:
  - Cuts out sequencing library adapters (config["adapters"],
    default "CTGTAGGCACC") using "cutadapt".
  - Extracts UMIs using "umi_tools extract", if requested
    (config["extract_umis"] == True), using a UMI-tools-compliant
    regular expression pattern (config["umi_regexp"]).
  - Removes rRNA or other contaminating reads by alignment to
    rRNA index file (config["rRNA_index"]) using "hisat2".
  - Aligns remaining reads to ORFs index file
    (config["orf_index"]). using "hisat2".
  - Trims 5' mismatches from reads and remove reads with more than 2
    mismatches using trim_5p_mismatch.py.
  - Outputs UMI groups pre-deduplication using "umi_tools group" if
    requested (config["dedup_umis"] == True and
    config["group_umis"] == True)
  - Deduplicates UMIs using "umi_tools dedup", if requested
    (config["dedup_umis"] == True)
  - Outputs UMI groups post-deduplication using "umi_tools group" if
    requested (config["dedup_umis"] == True and
    config["group_umis"] == True)
  - Exports bedgraph files for plus and minus strands, if requested
    (config["make_bedgraph"] == True) using "bedtools genomecov".
  - Makes length-sensitive alignments in compressed h5 format using
    "bam_to_h5.R".
  - Generates summary statistics, and analyses and QC plots for both
    RPF and mRNA datasets using "generate_stats_figs.R". This
    includes estimated read counts, reads per base, and transcripts
    per million for each ORF in each sample.
* Collates TPMs across all processed fastq.gz files, using
  "collate_tpms.R".
* The workflow can parallelize partos of its operation over many
  processes (config["nprocesses"]):
  - This value is used to configure "hisat2", "samtools sort",
    "bam_to_h5.R" and "generate_stats_figs.R".
  - For "cutadapt", the number of available processors
    on the host will be used.
* Writes all intermediate files into a temporary directory
  (config["dir_tmp"]).
* Writes all output files into an output directory
  (config["dir_out"]).

Exit codes are as follows:

* EXIT_OK (0): Processing successfully completed.
* EXIT_CONFIG_ERROR (1): Errors occurred loading configuration.
* EXIT_INDEX_ERROR (2): Error occurred during indexing.
* EXIT_NO_DATA_ERROR (3): No sample files or multiplexed data files
  were provided.
* EXIT_DATA_CONFLICT_ERROR (4): Both sample files and multiplexed data
  files were provided.
* EXIT_NO_BARCODES_ERROR (5): Multiplexed data files were provided but
  no barcodes file.
* EXIT_DATA_ERROR (6): No data file was processed successfully.
* EXIT_COLLATION_ERROR (7): Error occurred during TPMs collation.

Commands that are submitted to bash are recorded within a
file specified by a cmd_file configuration parameter.

If --dry-run is provided then the commands submitted to bash will not
be executed. This can be useful for seeing what commands will be run
without actually running them.
"""

import collections
from datetime import datetime
import errno
import logging
import os
import os.path
import sys
import yaml
from riboviz import process_utils
from riboviz import logging_utils
from riboviz.utils import value_in_dict


EXIT_OK = 0
""" Processing successfully completed. """
EXIT_CONFIG_ERROR = 1
""" Errors occurred loading configuration. """
EXIT_INDEX_ERROR = 2
""" Error occurred during indexing. """
EXIT_NO_DATA_ERROR = 3
""" No sample files or multiplexed data files were provided. """
EXIT_DATA_CONFLICT_ERROR = 4
""" Both sample files and multiplexed data files were provided. """
EXIT_NO_BARCODES_ERROR = 5
""" Multiplexed data files were provided but no barcodes file. """
EXIT_DATA_ERROR = 6
""" No data file was processed successfully. """
EXIT_COLLATION_ERROR = 7
""" Error occurred during TPMs collation. """


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


def demultiplex_fastq(fastq, barcodes, deplex_dir, log_file,
                      cmd_config):
    """
    Demultiplex reads.

    :param fastq: FASTQ file to demultiplex
    :type fastq: str or unicode
    :param barcodes: Sample sheet filename, tab-delimited
    text format with SampleID and TagRead columns, where TagReads are
    the barcodes to use to demultiplex fastq
    :type barcodes: str or unicode
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
           "-r1", fastq, "-ss", barcodes, "-o", deplex_dir,
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


def process_sample(sample, fastq, r_rna_index, orf_index, is_trimmed,
                   config, py_scripts, r_scripts, output_config,
                   cmd_config):
    """
    Process a single FASTQ sample file.

    :param sample: Sample name
    :type sample: str or unicode
    :param fastq: Sample FASTQ file
    :type fastq: str or unicode
    :param r_rna_index: Prefix of rRNA HT2 index files
    :type r_rna_index: str or unicode
    :param orf_index: Prefix of ORF HT2 index files
    :type orf_index: str or unicode
    :param is_trimmed: Have adapters been cut and barcodes and UMIs
    extracted already?
    :type are_trimmed: bool
    :param config: RiboViz configuration
    :type config: dict
    :param output_config: Output directories
    :type output_config: OutputConfigTuple
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :param py_scripts: Python scripts directory
    :type py_scripts: str or unicode
    :param r_scripts: R scripts directory
    :type r_scripts: str or unicode
    :raise FileNotFoundError: if fastq, other files or a third-party
    tool cannot be found
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    :raise KeyError: if config is missing required configuration
    """
    LOGGER.info("Processing sample: %s", sample)
    step = 1
    if not os.path.exists(fastq):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                fastq)
    LOGGER.info("Processing file: %s", fastq)

    if value_in_dict("nprocesses", config):
        nprocesses = int(config["nprocesses"])
    else:
        nprocesses = 1

    if is_trimmed:
        LOGGER.info("Skipping adaptor trimming and barcode/UMI extraction")
        trim_fq = fastq
    else:
        log_file = get_sample_log_file(
            output_config.logs_dir, sample, "cutadapt", step)
        trim_fq = os.path.join(output_config.tmp_dir, sample + "_trim.fq")
        cut_adapters(config["adapters"], fastq, trim_fq, log_file,
                     cmd_config)
        step += 1

        is_extract_umis = value_in_dict("extract_umis", config)
        if is_extract_umis:
            extract_trim_fq = os.path.join(
                output_config.tmp_dir, sample + "_extract_trim.fq")
            log_file = get_sample_log_file(
                output_config.logs_dir, sample, "umi_tools_extract",
                step)
            extract_barcodes_umis(
                trim_fq, extract_trim_fq, config["umi_regexp"],
                log_file, cmd_config)
            trim_fq = extract_trim_fq
            step += 1

    non_r_rna_trim_fq = os.path.join(
        output_config.tmp_dir, sample + "_nonrRNA.fq")
    r_rna_map_sam = os.path.join(
        output_config.tmp_dir, sample + "_rRNA_map.sam")
    log_file = get_sample_log_file(
        output_config.logs_dir, sample, "hisat2_rrna", step)
    map_to_r_rna(trim_fq, r_rna_index, r_rna_map_sam,
                 non_r_rna_trim_fq, nprocesses, log_file, cmd_config)
    step += 1

    orf_map_sam = os.path.join(
        output_config.tmp_dir, sample + "_orf_map.sam")
    unaligned_fq = os.path.join(
        output_config.tmp_dir, sample + "_unaligned.fq")
    log_file = get_sample_log_file(
        output_config.logs_dir, sample, "hisat2_orf", step)
    map_to_orf(non_r_rna_trim_fq, orf_index, orf_map_sam,
               unaligned_fq, nprocesses, log_file, cmd_config)
    step += 1

    orf_map_sam_clean = os.path.join(
        output_config.tmp_dir, sample + "_orf_map_clean.sam")
    log_file = get_sample_log_file(
        output_config.logs_dir, sample, "trim_5p_mismatch", step)
    trim_5p_mismatches(orf_map_sam, orf_map_sam_clean, py_scripts,
                       log_file, cmd_config)
    step += 1

    log_file = get_sample_log_file(
        output_config.logs_dir, sample, "samtools_view_sort", step)
    sample_out_prefix = os.path.join(output_config.out_dir, sample)
    sample_out_bam = sample_out_prefix + ".bam"
    sort_bam(orf_map_sam_clean, sample_out_bam, nprocesses, log_file,
             cmd_config)
    step += 1

    log_file = get_sample_log_file(
        output_config.logs_dir, sample, "samtools_index", step)
    index_bam(sample_out_bam, log_file, cmd_config)
    step += 1

    is_dedup_umis = value_in_dict("dedup_umis", config)
    if is_dedup_umis:
        if not is_extract_umis:
            LOGGER.warning(
                "WARNING: dedup_umis was TRUE but extract_umis was FALSE.")
        is_group_umis = value_in_dict("group_umis", config)
        if is_group_umis:
            umi_groups = os.path.join(
                output_config.tmp_dir, sample + "_pre_dedup_groups.tsv")
            log_file = get_sample_log_file(
                output_config.logs_dir, sample, "umi_tools_group",
                step)
            group_umis(sample_out_bam, umi_groups, log_file, cmd_config)
            step += 1

        sample_dedup_bam = sample_out_prefix + "_dedup.bam"
        log_file = get_sample_log_file(
            output_config.logs_dir, sample, "umi_tools_dedup", step)
        dedup_stats_prefix = os.path.join(
            output_config.tmp_dir, sample + "_dedup_stats")
        deduplicate_umis(sample_out_bam, sample_dedup_bam,
                         dedup_stats_prefix, log_file, cmd_config)
        step += 1

        log_file = get_sample_log_file(
            output_config.logs_dir, sample, "samtools_index", step)
        index_bam(sample_dedup_bam, log_file, cmd_config)
        sample_out_bam = sample_dedup_bam
        step += 1

        if is_group_umis:
            umi_groups = os.path.join(
                output_config.tmp_dir, sample + "_post_dedup_groups.tsv")
            log_file = get_sample_log_file(
                output_config.logs_dir, sample, "umi_tools_group",
                step)
            group_umis(sample_out_bam, umi_groups, log_file, cmd_config)
            step += 1

    is_make_bedgraph = value_in_dict("make_bedgraph", config)
    if is_make_bedgraph:
        log_file = get_sample_log_file(
            output_config.logs_dir, sample, "bedtools_genome_cov_plus", step)
        plus_bedgraph = sample_out_prefix + "_plus.bedgraph"
        make_bedgraph(sample_out_bam, plus_bedgraph, True, log_file,
                      cmd_config)
        step += 1

        log_file = get_sample_log_file(
            output_config.logs_dir, sample, "bedtools_genome_cov_minus", step)
        minus_bedgraph = sample_out_prefix + "_minus.bedgraph"
        make_bedgraph(sample_out_bam, minus_bedgraph, False, log_file,
                      cmd_config)
        step += 1

    orf_gff_file = config["orf_gff_file"]
    if not os.path.exists(orf_gff_file):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                orf_gff_file)
    log_file = get_sample_log_file(
        output_config.logs_dir, sample, "bam_to_h5", step)
    sample_out_h5 = sample_out_prefix + ".h5"
    bam_to_h5(sample_out_bam, sample_out_h5, orf_gff_file, config,
              nprocesses, r_scripts, log_file, cmd_config)
    step += 1

    log_file = get_sample_log_file(
        output_config.logs_dir, sample, "generate_stats_figs", step)
    generate_stats_figs(sample_out_h5, output_config.out_dir,
                        sample_out_prefix, config, nprocesses,
                        r_scripts, log_file, cmd_config)

    LOGGER.info("Finished processing sample: %s", fastq)


def process_samples(samples, in_dir, r_rna_index, orf_index, config,
                    py_scripts, r_scripts, output_config,
                    cmd_config):
    """
    Process FASTQ sample files. Any exceptions in the processing of
    any sample are logged but are not thrown from this function.

    :param samples: Sample names and files
    :type samples: dict
    :param in_dir: Directory with sample files
    :type in_dir: str or unicode
    :param r_rna_index: Prefix of rRNA HT2 index files
    :type r_rna_index: str or unicode
    :param orf_index: Prefix of ORF HT2 index files
    :type orf_index: str or unicode
    :param config: RiboViz configuration
    :type config: dict
    :param output_config: Output directories
    :type output_config: OutputConfigTuple
    :param cmd_config: Command-related configuration
    :type cmd_config: CmdConfigTuple
    :param py_scripts: Python scripts directory
    :type py_scripts: str or unicode
    :param r_scripts: R scripts directory
    :type r_scripts: str or unicode
    :return: names of successfully-processed samples
    :rtype: list(str or unicode)
    """
    LOGGER.info("Processing samples")
    successes = []
    num_samples = len(samples)
    for sample in list(samples.keys()):
        try:
            fastq = os.path.join(in_dir, samples[sample])
            process_sample(sample,
                           fastq,
                           r_rna_index,
                           orf_index,
                           False,
                           config,
                           py_scripts,
                           r_scripts,
                           output_config,
                           cmd_config)
            successes.append(sample)
        except FileNotFoundError as e:
            logging.error("File not found: %s", e.filename)
        except Exception:
            logging.error("Problem processing sample: %s", sample)
            exc_type, _, _ = sys.exc_info()
            logging.exception(exc_type.__name__)
        LOGGER.info("Finished processing %d samples, %d failed",
                    num_samples,
                    num_samples - len(successes))
    return successes


def prep_riboviz(py_scripts, r_scripts, config_yaml, is_dry_run=False):
    """
    Run the RiboViz workflow.

    Exit codes are as follows:

    * EXIT_OK (0): Processing successfully completed.
    * EXIT_CONFIG_ERROR (1): Errors occurred loading configuration.
    * EXIT_INDEX_ERROR (2): Error occurred during indexing.
    * EXIT_NO_DATA_ERROR (3): No sample files or multiplexed data
      files were provided.
    * EXIT_DATA_CONFLICT_ERROR (4): Both sample files and multiplexed
      data files were provided.
    * EXIT_NO_BARCODES_ERROR (5): Multiplexed data files were provided
      but no barcodes file.
    * EXIT_DATA_ERROR (6): No data file was processed successfully.
    * EXIT_COLLATION_ERROR (7): Error occurred during TPMs collation.

    :param py_scripts: Python scripts directory
    :type py_scripts: str or unicode
    :param r_scripts: R scripts directory
    :type r_scripts: str or unicode
    :param config_yaml: YAML configuration file path
    :type config_yaml: str or unicode
    :param is_dry_run: Don't execute workflow commands (useful for
    seeing what commands would be executed)
    :type is_dry_run: bool
    :return: exit code
    :rtype: int
    """
    LOGGER.info("Running under Python: %s", sys.version)
    LOGGER.info("Configuration file: %s", config_yaml)

    LOGGER.info("Load configuration: %s", config_yaml)
    try:
        with open(config_yaml, 'r') as f:
            config = yaml.load(f, yaml.SafeLoader)
    except FileNotFoundError as e:
        logging.error("File not found: %s", e.filename)
        return EXIT_CONFIG_ERROR
    except Exception:
        logging.error("Problem reading: %s", config_yaml)
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_CONFIG_ERROR

    if value_in_dict("cmd_file", config):
        cmd_file = config["cmd_file"]
    else:
        cmd_file = "run_riboviz_vignette.sh"
    LOGGER.info("Command file: %s", cmd_file)
    if os.path.exists(cmd_file):
        os.remove(cmd_file)
    cmd_config = CmdConfigTuple(cmd_file, is_dry_run)

    try:
        in_dir = config["dir_in"]
        output_config = setup_output_directories(config, cmd_config)
    except KeyError as e:
        logging.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except Exception:
        logging.error(("Problem configuring directories"))
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_CONFIG_ERROR

    LOGGER.info("Build indices for alignment, if necessary/requested")
    try:
        r_rna_fasta = config["rRNA_fasta"]
        orf_fasta = config["orf_fasta"]
        r_rna_index = os.path.join(
            output_config.index_dir, config["rRNA_index"])
        orf_index = os.path.join(
            output_config.index_dir, config["orf_index"])

        is_build_indices = value_in_dict("build_indices", config)
        if is_build_indices:
            log_file = os.path.join(
                output_config.logs_dir, "hisat2_build_r_rna.log")
            build_indices(r_rna_fasta, r_rna_index, log_file,
                          cmd_config)
            log_file = os.path.join(
                output_config.logs_dir, "hisat2_build_orf.log")
            build_indices(orf_fasta, orf_index, log_file, cmd_config)
    except KeyError as e:
        logging.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except FileNotFoundError as e:
        logging.error("File not found: %s", e.filename)
        return EXIT_INDEX_ERROR
    except Exception:
        logging.error("Problem creating indices")
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_INDEX_ERROR

    is_sample_files = value_in_dict("fq_files", config)
    is_multiplex_files = value_in_dict("multiplex_fq_files", config)
    is_barcodes = value_in_dict("barcodes", config)
    if not is_sample_files and not is_multiplex_files:
        LOGGER.error("No sample files (fq_files) or multiplexed files (multiplex_fq_files) are specified.")
        return EXIT_NO_DATA_ERROR
    elif is_sample_files and is_multiplex_files:
        LOGGER.error("Both sample files (fq_files) and multiplexed files (multiplex_fq_files) were specified.")
        return EXIT_DATA_CONFLICT_ERROR
    elif is_multiplex_files and not is_barcodes:
        LOGGER.error("Multiplexed files (multiplex_fq_files) are specified but no barcodes (barcodes) file.")
        return EXIT_NO_BARCODES_ERROR

    if is_sample_files:
        # Process sample files
        samples = config["fq_files"]
        processed_samples = process_samples(
            samples, in_dir, r_rna_index, orf_index, config,
            py_scripts, r_scripts, output_config, cmd_config)
        if not processed_samples:
            return EXIT_DATA_ERROR
        try:
            log_file = os.path.join(
                output_config.logs_dir, "collate_tpms.log")
            collate_tpms(output_config.out_dir, processed_samples,
                         r_scripts, log_file, cmd_config)
        except Exception:
            logging.error(("Problem collating TPMs"))
            exc_type, _, _ = sys.exc_info()
            logging.exception(exc_type.__name__)
            return EXIT_COLLATION_ERROR
    else:
        # Process multiplexed files
        LOGGER.info("WIP: multiplexed file processing")

        try:
            barcodes = os.path.join(in_dir, config["barcodes"])
            if not os.path.exists(barcodes):
                raise FileNotFoundError(errno.ENOENT,
                                        os.strerror(errno.ENOENT),
                                        barcodes)
        except FileNotFoundError as e:
            logging.error("File not found: %s", e.filename)
            return EXIT_CONFIG_ERROR
        except Exception:
            logging.error("Problem processing: %s", barcodes)
            exc_type, _, _ = sys.exc_info()
            logging.exception(exc_type.__name__)
            return EXIT_CONFIG_ERROR

        multiplex_files = config["multiplex_fq_files"]
        # WIP: take first file only
        multiplex_file = multiplex_files[0]
        multiplex_name = os.path.splitext(os.path.basename(multiplex_file))[0]
        multiplex_file = os.path.join(in_dir, multiplex_file)
        LOGGER.info("Processing file: %s", multiplex_file)
        try:
            if not os.path.exists(multiplex_file):
                raise FileNotFoundError(errno.ENOENT,
                                        os.strerror(errno.ENOENT),
                                        multiplex_file)
            trim_fq = os.path.join(
                output_config.tmp_dir, multiplex_name + "_trim.fq")
            log_file = os.path.join(
                output_config.logs_dir, "cutadapt.log")
            cut_adapters(config["adapters"], multiplex_file, trim_fq,
                         log_file, cmd_config)

            extract_trim_fq = os.path.join(
                output_config.tmp_dir, multiplex_name + "_extract_trim.fq")
            log_file = os.path.join(
                output_config.logs_dir, "umi_tools_extract.log")
            extract_barcodes_umis(
                trim_fq, extract_trim_fq, config["umi_regexp"], log_file,
                cmd_config)

            deplex_dir = os.path.join(
                output_config.tmp_dir, multiplex_name + "_deplex")
            log_file = os.path.join(
                output_config.logs_dir, "demultiplex_fastq.log")
            demultiplex_fastq(extract_trim_fq, barcodes, deplex_dir,
                              log_file, cmd_config)

        except FileNotFoundError as e:
            logging.error("File not found: %s", e.filename)
            return EXIT_DATA_ERROR
        except Exception:
            logging.error("Problem processing: %s", multiplex_file)
            exc_type, _, _ = sys.exc_info()
            logging.exception(exc_type.__name__)
            return EXIT_DATA_ERROR

    LOGGER.info("Completed")
    return EXIT_OK


if __name__ == "__main__":
    # Assume RiboViz Python scripts are peers in same directory.
    py_scripts_arg = os.path.dirname(os.path.realpath(__file__))
    sys.argv.pop(0)  # Remove program.
    is_dry_run_arg = (sys.argv[0] == "--dry-run")
    if is_dry_run_arg:
        sys.argv.pop(0)
    r_scripts_arg = sys.argv[0]
    config_yaml_arg = sys.argv[1]
    exit_code = prep_riboviz(py_scripts_arg,
                             r_scripts_arg,
                             config_yaml_arg,
                             is_dry_run_arg)
    sys.exit(exit_code)
