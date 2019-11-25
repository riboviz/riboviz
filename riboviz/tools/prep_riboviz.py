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
* EXIT_NO_SAMPLES_ERROR (3): No samples were provided.
* EXIT_SAMPLES_ERROR (4): No sample was processed successfully.
* EXIT_COLLATION_ERROR (5): Error occurred during TPMs collation.

Commands that are submitted to bash are recorded within a
file specified by a cmd_file configuration parameter.

If --dry-run is provided then the commands submitted to bash will not
be executed. This can be useful for seeing what commands will be run
without actually running them.
"""

from datetime import datetime
import errno
import logging
import os
import os.path
import sys
import yaml
from riboviz import process_utils
from riboviz import logging_utils


EXIT_OK = 0
""" Processing successfully completed. """
EXIT_CONFIG_ERROR = 1
""" Errors occurred loading configuration. """
EXIT_INDEX_ERROR = 2
""" Error occurred during indexing. """
EXIT_NO_SAMPLES_ERROR = 3
""" No samples were provided. """
EXIT_SAMPLES_ERROR = 4
""" No sample was processed successfully. """
EXIT_COLLATION_ERROR = 5
""" Error occurred during TPMs collation. """


logging_utils.configure_logging()
LOGGER = logging.getLogger(__name__)
""" Logger """


def build_indices(fasta, ht_prefix, file_type, logs_dir, cmd_file,
                  dry_run=False):
    """
    Build indices for alignment via invocation of hisat2-build.
    Index files have name <ht_prefix>.<N>.ht2.

    :param fasta: FASTA file to be indexed
    :type fasta: str or unicode
    :param ht_prefix: Prefix of HT2 index files
    :type ht_prefix: str or unicode
    :param logs_dir Log files directory
    :type logs_dir: str or unicode
    :param file_type: Type of file being indexed, used for logging
    :type file_type: str or unicode
    :param cmd_file: File to log command to, if not None
    :type cmd_file: str or unicode
    :param dry_run: Don't execute workflow commands (useful for seeing
    what commands would be executed)
    :type dry_run: bool
    :raise FileNotFoundError: if hisat2-build cannot be found
    :raise AssertionError: if hisat2-build returns non-zero exit
    code
    """
    log_file = os.path.join(logs_dir, "hisat2_build_%s.log" % file_type)
    LOGGER.info("Create indices: %s. Log: %s", file_type, log_file)
    cmd = ["hisat2-build", fasta, ht_prefix]
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)


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


def group_umis(sample_bam,
               tmp_dir,
               sample,
               tag,
               step,
               logs_dir,
               cmd_file,
               dry_run):
    """
    Run "umi_tools group" on a BAM file for a sample.

    :param sample_bam: Sample BAM file
    :type sample_bam: str or unicode
    :param tmp_dir Temporary files directory
    :type tmp_dir: str or unicode
    :param sample_prefix: Sample name.
    :type sample: str or unicode
    :param tag: Tag used as part of groups file name to distinguish it
    from other groups files for the sample.
    :type tag: str or unicode
    :param step: Current step in workflow
    :type step: str or unicode
    :param logs_dir Log files directory
    :type logs_dir: str or unicode
    :param cmd_file: File to log command to, if not None
    :type cmd_file: str or unicode
    :param dry_run: Don't execute workflow commands (useful for seeing
    what commands would be executed)
    :type dry_run: bool
    :raise OSError: if a third-party tool cannot be found
    :raise FileNotFoundError: if umi_tools or a third-party tool
    cannot be found
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    """
    sample_umi_groups = os.path.join(
        tmp_dir, sample + "_" + tag + "_groups.tsv")
    log_file = get_sample_log_file(logs_dir,
                                   sample,
                                   "umi_tools_group",
                                   step)
    LOGGER.info("Identify UMI groups. Log: %s", log_file)
    cmd = ["umi_tools", "group", "-I", sample_bam,
           "--group-out", sample_umi_groups]
    process_utils.run_logged_command(cmd,
                                     log_file,
                                     cmd_file,
                                     dry_run)


def process_sample(sample,
                   fastq,
                   r_rna_index,
                   orf_index,
                   config,
                   py_scripts,
                   r_scripts,
                   tmp_dir,
                   out_dir,
                   logs_dir,
                   cmd_file,
                   dry_run=False):
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
    :param config: RiboViz configuration
    :type config: dict
    :param py_scripts: Directory with RiboViz Python scripts
    :type py_scripts: str or unicode
    :param r_scripts:  Directory with RiboViz R scripts
    :type r_scripts: str or unicode
    :param tmp_dir Temporary files directory
    :type tmp_dir: str or unicode
    :param out_dir Output files directory
    :type out_dir: str or unicode
    :param logs_dir Log files directory
    :type logs_dir: str or unicode
    :param cmd_file: File to log command to, if not None
    :type cmd_file: str or unicode
    :param dry_run: Don't execute workflow commands (useful for seeing
    what commands would be executed)
    :type dry_run: bool
    :raise FileNotFoundError: if fastq, other files or a third-party
    tool cannot be found
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    :raise KeyError: if config is missing required configuration
    """
    LOGGER.info("Processing sample: %s", sample)
    step = 1

    if not os.path.exists(fastq):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                fastq)
    LOGGER.info("Processing file: %s", fastq)

    if "nprocesses" not in config:
        nprocesses = 1
    else:
        nprocesses = config["nprocesses"]

    # Trimmed reads.
    log_file = get_sample_log_file(logs_dir, sample, "cutadapt", step)
    LOGGER.info("Cut Illumina adapters. Log: %s", log_file)
    trim_fq = os.path.join(tmp_dir, sample + "_trim.fq")
    cmd = ["cutadapt", "--trim-n", "-O", "1", "-m", "5",
           "-a", config["adapters"], "-o", trim_fq, fastq]
    # cutadapt allows all available processors to be requested.
    cmd += ["-j", str(0)]
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
    step += 1

    extract_umis = "extract_umis" in config and config["extract_umis"]
    if extract_umis:
        extract_trim_fq = os.path.join(tmp_dir, sample + "_extract_trim.fq")
        log_file = get_sample_log_file(logs_dir,
                                       sample,
                                       "umi_tools_extract",
                                       step)
        LOGGER.info("Extract UMIs. Log: %s", log_file)
        pattern_parameter = "--bc-pattern=" + config["umi_regexp"]
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
        cmd = ["umi_tools", "extract", "-I", trim_fq,
               pattern_parameter,
               "--extract-method=regex",
               "-S", extract_trim_fq]
        cmd_to_log = ["--bc-pattern=" + "\"" + config["umi_regexp"] + "\""
                      if c == pattern_parameter else c for c in cmd]
        process_utils.run_logged_command(cmd,
                                         log_file,
                                         cmd_file,
                                         dry_run,
                                         cmd_to_log)
        trim_fq = extract_trim_fq
        step += 1

    log_file = get_sample_log_file(logs_dir, sample, "hisat2_rrna", step)
    LOGGER.info("Map reads to rRNA. Log: %s", log_file)
    # Trimmed non-rRNA reads.
    non_r_rna_trim_fq = os.path.join(tmp_dir, sample + "_nonrRNA.fq")
    # rRNA-mapped reads.
    r_rna_map_sam = os.path.join(tmp_dir, sample + "_rRNA_map.sam")
    cmd = ["hisat2", "-p", str(nprocesses), "-N", "1",
           "--un", non_r_rna_trim_fq, "-x", r_rna_index,
           "-S", r_rna_map_sam, "-U", trim_fq]
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
    step += 1

    log_file = get_sample_log_file(logs_dir, sample, "hisat2_orf", step)
    LOGGER.info("Map to ORFs with up to 2 alignments. Log: %s", log_file)
    # ORF-mapped reads.
    orf_map_sam = os.path.join(tmp_dir, sample + "_orf_map.sam")
    # Unaligned reads.
    unaligned_fq = os.path.join(tmp_dir, sample + "_unaligned.fq")
    cmd = ["hisat2", "-p", str(nprocesses), "-k", "2",
           "--no-spliced-alignment", "--rna-strandness",
           "F", "--no-unal", "--un", unaligned_fq,
           "-x", orf_index, "-S", orf_map_sam,
           "-U", non_r_rna_trim_fq]
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
    step += 1

    log_file = get_sample_log_file(logs_dir, sample, "trim_5p_mismatch", step)
    LOGGER.info("Trim 5' mismatched nt and remove reads with >1 mismatch. Log: %s",
                log_file)
    # ORF-mapped reads.
    orf_map_sam_clean = os.path.join(tmp_dir,
                                     sample +
                                     "_orf_map_clean.sam")
    cmd = ["python", os.path.join(py_scripts, "trim_5p_mismatch.py"),
           "-mm", "2", "-in", orf_map_sam,
           "-out", orf_map_sam_clean]
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
    step += 1

    log_file = get_sample_log_file(logs_dir,
                                   sample,
                                   "samtools_view_sort",
                                   step)
    LOGGER.info("Convert SAM to BAM and sort on genome. Log: %s", log_file)
    sample_out_prefix = os.path.join(out_dir, sample)
    sample_out_bam = sample_out_prefix + ".bam"
    cmd_view = ["samtools", "view", "-b", orf_map_sam_clean]
    cmd_sort = ["samtools", "sort", "-@", str(nprocesses),
                "-O", "bam", "-o", sample_out_bam, "-"]
    process_utils.run_logged_pipe_command(
        cmd_view,
        cmd_sort,
        log_file,
        cmd_file,
        dry_run)
    step += 1

    log_file = get_sample_log_file(logs_dir, sample, "samtools_index", step)
    LOGGER.info("Index BAM file. Log: %s", log_file)
    cmd = ["samtools", "index", sample_out_bam]
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
    step += 1

    if "dedup_umis" in config and config["dedup_umis"]:
        LOGGER.info("Deduplicate using UMIs. Log: %s", log_file)
        if not extract_umis:
            LOGGER.warning("WARNING: dedup_umis was TRUE but extract_umis was FALSE.")
        if "group_umis" in config and config["group_umis"]:
            group_umis(sample_out_bam,
                       tmp_dir,
                       sample,
                       "pre_dedup",
                       step,
                       logs_dir,
                       cmd_file,
                       dry_run)
            step += 1
        sample_dedup_bam = sample_out_prefix + "_dedup.bam"
        log_file = get_sample_log_file(logs_dir,
                                       sample,
                                       "umi_tools_dedup",
                                       step)
        dedup_stats_prefix = os.path.join(tmp_dir, sample + "_dedup_stats")
        LOGGER.info("Deduplicate. Log: %s", log_file)
        cmd = ["umi_tools", "dedup", "-I", sample_out_bam,
               "-S", sample_dedup_bam,
               "--output-stats=" + dedup_stats_prefix]
        process_utils.run_logged_command(cmd,
                                         log_file,
                                         cmd_file,
                                         dry_run)
        step += 1
        log_file = get_sample_log_file(logs_dir,
                                       sample,
                                       "samtools_index",
                                       step)
        LOGGER.info("Index BAM file. Log: %s", log_file)
        cmd = ["samtools", "index", sample_dedup_bam]
        process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
        sample_out_bam = sample_dedup_bam
        step += 1
        if "group_umis" in config and config["group_umis"]:
            group_umis(sample_dedup_bam,
                       tmp_dir,
                       sample,
                       "post_dedup",
                       step,
                       logs_dir,
                       cmd_file,
                       dry_run)
            step += 1

    if config["make_bedgraph"]:
        LOGGER.info("Record transcriptome coverage as a bedgraph")
        log_file = get_sample_log_file(
            logs_dir, sample, "bedtools_genome_cov_plus", step)
        LOGGER.info("Calculate coverage for plus strand. Log: %s",
                    log_file)
        cmd = ["bedtools", "genomecov", "-ibam", sample_out_bam,
               "-trackline", "-bga", "-5", "-strand", "+"]
        plus_bedgraph = sample_out_prefix + "_plus.bedgraph"
        process_utils.run_logged_redirect_command(
            cmd,
            plus_bedgraph,
            log_file,
            cmd_file,
            dry_run)
        step += 1

        log_file = get_sample_log_file(
            logs_dir, sample, "bedtools_genome_cov_minus", step)
        LOGGER.info("Calculate coverage for minus strand. Log: %s",
                    log_file)
        cmd = ["bedtools", "genomecov", "-ibam", sample_out_bam,
               "-trackline", "-bga", "-5", "-strand", "-"]
        minus_bedgraph = sample_out_prefix + "_minus.bedgraph"
        process_utils.run_logged_redirect_command(
            cmd,
            minus_bedgraph,
            log_file,
            cmd_file,
            dry_run)
        step += 1

    log_file = get_sample_log_file(logs_dir, sample, "bam_to_h5", step)
    LOGGER.info("Make length-sensitive alignments in H5 format. Log: %s",
                log_file)
    sample_out_h5 = sample_out_prefix + ".h5"
    second_id = config["SecondID"]
    if second_id is None:
        second_id = "NULL"
    orf_gff_file = config["orf_gff_file"]
    if not os.path.exists(orf_gff_file):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                orf_gff_file)
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "bam_to_h5.R"),
           "--Ncores=" + str(nprocesses),
           "--MinReadLen=" + str(config["MinReadLen"]),
           "--MaxReadLen=" + str(config["MaxReadLen"]),
           "--Buffer=" + str(config["Buffer"]),
           "--PrimaryID=" + config["PrimaryID"],
           "--SecondID=" + second_id,
           "--dataset=" + config["dataset"],
           "--bamFile=" + sample_out_bam,
           "--hdFile=" + sample_out_h5,
           "--orf_gff_file=" + orf_gff_file,
           "--ribovizGFF=" + str(config["ribovizGFF"]),
           "--StopInCDS=" + str(config["StopInCDS"])]
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
    step += 1

    log_file = get_sample_log_file(logs_dir,
                                   sample,
                                   "generate_stats_figs",
                                   step)
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
           "--hdFile=" + sample_out_h5,
           "--out_prefix=" + sample_out_prefix,
           "--orf_fasta=" + config["orf_fasta"],
           "--rpf=" + str(config["rpf"]),
           "--dir_out=" + out_dir,
           "--do_pos_sp_nt_freq=" + str(config["do_pos_sp_nt_freq"])]

    for flag in ["t_rna", "codon_pos", "features_file", "orf_gff_file"]:
        if flag in config and config[flag] is not None:
            flag_file = config[flag]
            if not os.path.exists(flag_file):
                raise FileNotFoundError(errno.ENOENT,
                                        os.strerror(errno.ENOENT),
                                        flag_file)
            cmd.append("--" + flag + "=" + flag_file)
    for flag in ["count_threshold", "asite_disp_length_file"]:

        if flag in config and config[flag] is not None:
            cmd.append("--" + flag + "=" + str(config[flag]))
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
    LOGGER.info("Finished processing sample: %s", fastq)


def collate_tpms(out_dir, samples, r_scripts, logs_dir, cmd_file,
                 dry_run=False):
    """
    Collate TPMs across sample results.

    :param out_dir Output files directory
    :type out_dir: str or unicode
    :param samples: Sample names
    :type samples: list(str or unicode)
    :param r_scripts:  Directory with RiboViz R scripts
    :type r_scripts: str or unicode
    :param logs_dir Log files directory
    :type logs_dir: str or unicode
    :param cmd_file: File to log command to, if not None
    :type cmd_file: str or unicode
    :param dry_run: Don't execute workflow commands (useful for seeing
    what commands would be executed)
    :type dry_run: bool
    :raise FileNotFoundError: if Rscript cannot be found
    :raise AssertionError: if collate_tpms.R returns non-zero exit
    code
    """
    log_file = os.path.join(logs_dir, "collate_tpms.log")
    LOGGER.info("Collate TPMs across all processed samples. Log: %s",
                log_file)
    cmd = ["Rscript",
           "--vanilla",
           os.path.join(r_scripts, "collate_tpms.R"),
           "--dir_out=" + out_dir]
    cmd += samples
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)


def prep_riboviz(py_scripts, r_scripts, config_yaml, dry_run=False):
    """
    Run the RiboViz workflow.

    Exit codes are as follows:

    * EXIT_OK (0): Processing successfully completed.
    * EXIT_CONFIG_ERROR (1): Errors occurred loading configuration.
    * EXIT_INDEX_ERROR (2): Error occurred during indexing.
    * EXIT_NO_SAMPLES_ERROR (3): No samples were provided.
    * EXIT_SAMPLES_ERROR (4): No sample was processed successfully.
    * EXIT_COLLATION_ERROR (5): Error occurred during TPMs collation.

    :param py_scripts: Directory with RiboViz Python scripts
    :type py_scripts: str or unicode
    :param r_scripts:  Directory with RiboViz R scripts
    :type r_scripts: str or unicode
    :param config_yaml: YAML configuration file path
    :type config_yaml: str or unicode
    :param dry_run: Don't execute workflow commands (useful for seeing
    what commands would be executed)
    :type dry_run: bool
    :return: exit code
    :rtype: int
    """
    LOGGER.info("Running under Python: %s", sys.version)
    LOGGER.info("Configuration file: %s", config_yaml)

    # Extract configuration.
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

    # Set up command file.
    if "cmd_file" in config:
        cmd_file = config["cmd_file"]
    else:
        cmd_file = "run_riboviz_vignette.sh"
    LOGGER.info("Command file: %s", cmd_file)
    if os.path.exists(cmd_file):
        os.remove(cmd_file)

    # Set up directories
    try:
        in_dir = config["dir_in"]
        index_dir = config["dir_index"]
        tmp_dir = config["dir_tmp"]
        out_dir = config["dir_out"]
        base_logs_dir = config["dir_logs"]
        logs_dir = os.path.join(base_logs_dir,
                                datetime.now().strftime('%Y%m%d-%H%M%S'))
        dirs = [index_dir, tmp_dir, out_dir, base_logs_dir, logs_dir]
        for d in dirs:
            with open(cmd_file, "a") as f:
                f.write("mkdir -p %s\n" % d)
        if not dry_run:
            for d in dirs:
                if not os.path.exists(d):
                    os.makedirs(d)
    except KeyError as e:
        logging.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except Exception:
        logging.error(("Problem configuring logs directory"))
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_CONFIG_ERROR

    LOGGER.info("Build indices for alignment, if necessary/requested")
    try:
        r_rna_fasta = config["rRNA_fasta"]
        if not os.path.exists(r_rna_fasta):
            raise FileNotFoundError(errno.ENOENT,
                                    os.strerror(errno.ENOENT),
                                    r_rna_fasta)
        orf_fasta = config["orf_fasta"]
        if not os.path.exists(orf_fasta):
            raise FileNotFoundError(errno.ENOENT,
                                    os.strerror(errno.ENOENT),
                                    orf_fasta)
        r_rna_index = os.path.join(index_dir,
                                   config["rRNA_index"])
        orf_index = os.path.join(index_dir,
                                 config["orf_index"])
        if config["build_indices"]:
            build_indices(r_rna_fasta,
                          r_rna_index,
                          "r_rna",
                          logs_dir,
                          cmd_file,
                          dry_run)
            build_indices(orf_fasta,
                          orf_index,
                          "orf",
                          logs_dir,
                          cmd_file,
                          dry_run)
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

    # Loop over sample fastq.gz files.
    LOGGER.info("Processing samples")
    if ("fq_files" not in config) or \
       (config["fq_files"] is None) or \
       (not config["fq_files"]):
        LOGGER.error("No samples are defined")
        return EXIT_NO_SAMPLES_ERROR
    samples = config["fq_files"]
    num_samples = len(config["fq_files"])
    successes = []
    for sample in list(samples.keys()):
        try:
            fastq = os.path.join(in_dir, samples[sample])
            process_sample(sample,
                           fastq,
                           r_rna_index,
                           orf_index,
                           config,
                           py_scripts,
                           r_scripts,
                           tmp_dir,
                           out_dir,
                           logs_dir,
                           cmd_file,
                           dry_run)
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
    if not successes:
        return EXIT_SAMPLES_ERROR

    # Collate TPMs across sample results.
    try:
        collate_tpms(out_dir, successes, r_scripts, logs_dir, cmd_file,
                     dry_run)
    except Exception:
        logging.error(("Problem collating TPMs"))
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_COLLATION_ERROR

    LOGGER.info("Completed")
    return EXIT_OK


if __name__ == "__main__":
    # Assume RiboViz Python scripts are peers in same directory.
    py_scripts_arg = os.path.dirname(os.path.realpath(__file__))
    sys.argv.pop(0)  # Remove program.
    dry_run = (sys.argv[0] == "--dry-run")
    if dry_run:
        sys.argv.pop(0)
    r_scripts_arg = sys.argv[0]
    config_yaml_arg = sys.argv[1]
    exit_code = prep_riboviz(py_scripts_arg,
                             r_scripts_arg,
                             config_yaml_arg,
                             dry_run)
    sys.exit(exit_code)
