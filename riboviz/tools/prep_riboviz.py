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

* - Read configuration info from YAML file.
* - Build hisat2 indices if requested (config["build_indices"] ==
*   True), in index directory (config["dir_index"]).
* - Process all fastq.gz files (config["dir_in"]).
* - Cut out sequencing library adapters ("CTGTAGGCACC" or
*   config["adapters").
* - Remove rRNA or other contaminating reads by hisat2 alignment to
*   rRNA index file (config["rRNA_index"])
* - Align remaining reads to ORFs or other hisat2 index file
*   (config["orf_index"]).
* - Trim 5' mismatches from reads and remove reads with more than 2
*   mismatches.
* - Parallelize over many processes (config["nprocesses"]):
*   - This value is used to configure hisat2, samtools sort,
*     bam_to_h5.R and generate_stats_figs.R.
*   - For cutadapt and Python 3, the number of available processors
*     on the host will be used.
*   - For cutadapt and Python 2, its default of 1 processor will be
*     used as cutadapt cannot run in parallel under Python 2.
* - Make length-sensitive alignments in compressed h5 format by
*   running "bam_to_h5.R".
* - Generate summary statistics, and analyses and QC plots for both
*   RPF and mRNA datasets, by running "generate_stats_figs.R".
* - Put all intermediate files into a temporary directory
*   (config["dir_tmp"]).
* - When finished, put useful output files into output directory
*   (config["dir_out"]).
* Optionally export bedgraph files for plus and minus strands, if
*   requested (config["make_bedgraph"] == True).
*  - bamfiles (config["dir_out"]/*.bam) are directly usable in genome
*  - browsers, etc.

Exit codes are as follows:

* EXIT_OK (0): Processing successfully completed.
* EXIT_CONFIG_ERROR (1): Errors occurred loading configuration.
* EXIT_INDEX_ERROR (2): Error occurred during indexing.
* EXIT_NO_SAMPLES_ERROR (3): No samples were provided.
* EXIT_SAMPLES_ERROR (4): No sample was processed successfully.
* EXIT_COLLATION_ERROR (5): Error occurred during TPMs collation.

If --dry-run is provided then the workflow commands will not be
executed. This can be useful for seeing what commands will be
run, in advance of running them.
"""

from datetime import datetime
import errno
import logging
import os
import os.path
import sys
import yaml
from riboviz import utils
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
logger = logging.getLogger(__name__)


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
    cmd = ["hisat2-build", fasta, ht_prefix]
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            os.path.join(logs_dir, "hisat2_build_" + file_type + ".log"))


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
    return os.path.join(
        logs_dir,
        sample + ("_%02d_" % index) + step + ".log")


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
    :raise FileNotFoundError: if fastq or a third-party tool cannot be
    found
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    :raise KeyError: if config is missing required configuration
    """
    logger.info("Processing sample: %s", sample)
    step = 1

    if not os.path.exists(fastq):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                fastq)
    logger.info("Processing file: %s", fastq)

    if "nprocesses" not in config:
        nprocesses = 1
    else:
        nprocesses = config["nprocesses"]

    # Cut illumina adapters.

    # Trimmed reads.
    trim_fq = os.path.join(tmp_dir, sample + "_trim.fq")
    cmd = ["cutadapt", "--trim-n", "-O", "1", "-m", "5",
           "-a", config["adapters"], "-o", trim_fq, fastq]
    py_major = sys.version_info.major
    if py_major == 3:
        # cutadapt and Python 3 allows all available processors to
        # be requested.
        cmd += ["-j", str(0)]
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            get_sample_log_file(logs_dir, sample, "cut_adapt", step))
    step += 1

    # Map reads to rRNA.

    # Trimmed non-rRNA reads.
    non_r_rna_trim_fq = os.path.join(tmp_dir, sample + "_nonrRNA.fq")
    # rRNA-mapped reads.
    r_rna_map_sam = os.path.join(tmp_dir, sample + "_rRNA_map.sam")
    cmd = ["hisat2", "-p", str(nprocesses), "-N", "1",
           "--un", non_r_rna_trim_fq, "-x", r_rna_index,
           "-S", r_rna_map_sam, "-U", trim_fq]
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            get_sample_log_file(logs_dir, sample, "hisat2_rrna", step))
    step += 1

    # Map to ORFs with (mostly) default settings, up to 2 alignments.

    # ORF-mapped reads.
    orf_map_sam = os.path.join(tmp_dir, sample + "_orf_map.sam")
    # Unaligned reads.
    unaligned_fq = os.path.join(tmp_dir, sample + "_unaligned.fq")
    cmd = ["hisat2", "-p", str(nprocesses), "-k", "2",
           "--no-spliced-alignment", "--rna-strandness",
           "F", "--no-unal", "--un", unaligned_fq,
           "-x", orf_index, "-S", orf_map_sam,
           "-U", non_r_rna_trim_fq]
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            get_sample_log_file(logs_dir, sample, "hisat2_orf", step))
    step += 1

    # Trim 5' mismatched nt and remove reads with >1 mismatch.

    # ORF-mapped reads.
    orf_map_sam_clean = os.path.join(tmp_dir,
                                     sample +
                                     "_orf_map_clean.sam")
    cmd = ["python", os.path.join(py_scripts, "trim_5p_mismatch.py"),
           "-mm", "2", "-in", orf_map_sam,
           "-out", orf_map_sam_clean]
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            get_sample_log_file(logs_dir, sample, "trim_5p_mismatch", step))
    step += 1

    # Convert SAM (text) output to BAM (compressed binary) and sort on
    # genome.

    sample_out_prefix = os.path.join(out_dir, sample)
    sample_out_bam = sample_out_prefix + ".bam"
    cmd_view = ["samtools", "view", "-b", orf_map_sam_clean]
    cmd_sort = ["samtools", "sort", "-@", str(nprocesses),
                "-O", "bam", "-o", sample_out_bam, "-"]
    logger.info("Running: %s | %s",
                utils.list_to_str(cmd_view),
                utils.list_to_str(cmd_sort))
    process_utils.log_pipe_command(cmd_view, cmd_sort, cmd_file)
    if not dry_run:
        process_utils.run_logged_pipe_command(
            cmd_view,
            cmd_sort,
            get_sample_log_file(logs_dir, sample, "samtools_view_sort", step))
    step += 1

    # Index BAM file.
    cmd = ["samtools", "index", sample_out_bam]
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            get_sample_log_file(logs_dir, sample, "samtools_index", step))
    step += 1

    if config["make_bedgraph"]:
        # Record transcriptome coverage as a bedgraph.
        # Calculate transcriptome coverage for plus strand.
        cmd = ["bedtools", "genomecov", "-ibam", sample_out_bam,
               "-trackline", "-bga", "-5", "-strand", "+"]
        plus_bedgraph = sample_out_prefix + "_plus.bedgraph"
        logger.info("Running: %s > %s",
                    utils.list_to_str(cmd),
                    plus_bedgraph)
        process_utils.log_redirect_command(cmd, plus_bedgraph, cmd_file)
        if not dry_run:
            process_utils.run_logged_redirect_command(
                cmd,
                plus_bedgraph,
                get_sample_log_file(
                    logs_dir, sample, "bedtools_genome_cov_plus", step))
        step += 1
        # Calculate transcriptome coverage for minus strand.
        cmd = ["bedtools", "genomecov", "-ibam", sample_out_bam,
               "-trackline", "-bga", "-5", "-strand", "-"]
        minus_bedgraph = sample_out_prefix + "_minus.bedgraph"
        logger.info("Running: %s > %s",
                    utils.list_to_str(cmd),
                    minus_bedgraph)
        process_utils.log_redirect_command(cmd, minus_bedgraph, cmd_file)
        if not dry_run:
            process_utils.run_logged_redirect_command(
                cmd,
                minus_bedgraph,
                get_sample_log_file(
                    logs_dir, sample, "bedtools_genome_cov_minus", step))
        step += 1
        logger.info("bedgraphs made on plus and minus strands")

    # Make length-sensitive alignments in H5 format.

    sample_out_h5 = sample_out_prefix + ".h5"
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
           "--bamFile=" + sample_out_bam,
           "--hdFile=" + sample_out_h5,
           "--orf_gff_file=" + config["orf_gff_file"],
           "--ribovizGFF=" + str(config["ribovizGFF"]),
           "--StopInCDS=" + str(config["StopInCDS"])]
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            get_sample_log_file(logs_dir, sample, "bam_to_h5", step))
    step += 1

    # Create summary statistics and analyses plots.

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
           "--orf_gff_file=" + config["orf_gff_file"],
           "--dir_out=" + out_dir,
           "--t_rna=" + config["t_rna"],
           "--codon_pos=" + config["codon_pos"],
           "--features_file=" + config["features_file"],
           "--do_pos_sp_nt_freq=" + str(config["do_pos_sp_nt_freq"])]
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            get_sample_log_file(logs_dir, sample, "generate_stats_figs", step))
    logger.info("Finished processing sample: %s", fastq)


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
    logger.info("Collating TPMs across all processed samples")
    cmd = ["Rscript",
           "--vanilla",
           os.path.join(r_scripts, "collate_tpms.R"),
           "--dir_out=" + out_dir]
    cmd += samples
    logger.info("Running: %s", utils.list_to_str(cmd))
    process_utils.log_command(cmd, cmd_file)
    if not dry_run:
        process_utils.run_logged_command(
            cmd,
            os.path.join(logs_dir, "collate_tpms.log"))


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
    logger.info("Running under Python: %s ", sys.version)

    cmd_file = "prep_riboviz.sh"
    if os.path.exists(cmd_file):
        os.remove(cmd_file)

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

    try:
        base_logs_dir = config["dir_logs"]
        logs_dir = os.path.join(base_logs_dir,
                                datetime.now().strftime('%Y%m%d-%H%M%S'))
        if not dry_run:
            if not os.path.exists(base_logs_dir):
                os.makedirs(base_logs_dir)
            if not os.path.exists(logs_dir):
                os.makedirs(logs_dir)
    except KeyError as e:
        logging.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except Exception:
        logging.error(("Problem configuring logs directory"))
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_CONFIG_ERROR

    # Build indices for alignment, if necessary/requested.
    try:
        index_dir = config["dir_index"]
        r_rna_index = os.path.join(index_dir,
                                   config["rRNA_index"])
        orf_index = os.path.join(index_dir,
                                 config["orf_index"])
        if config["build_indices"]:
            if not os.path.exists(index_dir) and not dry_run:
                os.makedirs(index_dir)
            build_indices(config["rRNA_fasta"],
                          r_rna_index,
                          "r_rna",
                          logs_dir,
                          cmd_file,
                          dry_run)
            logger.info("rRNA index built")
            build_indices(config["orf_fasta"],
                          orf_index,
                          "orf",
                          logs_dir,
                          cmd_file,
                          dry_run)
            logger.info("ORF index built")
    except KeyError as e:
        logging.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except Exception:
        logging.error("Problem creating indices")
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_INDEX_ERROR

    # Loop over sample fastq.gz files.
    logger.info("Processing samples")
    try:
        in_dir = config["dir_in"]
        tmp_dir = config["dir_tmp"]
        out_dir = config["dir_out"]
    except KeyError as e:
        logging.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    if not os.path.exists(tmp_dir) and not dry_run:
        os.makedirs(tmp_dir)
    if not os.path.exists(out_dir) and not dry_run:
        os.makedirs(out_dir)
    if ("fq_files" not in config) or \
       (config["fq_files"] is None) or \
       (not config["fq_files"]):
        logger.error("No samples are defined")
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
    logger.info("Finished processing %d samples, %d failed",
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

    logger.info("Completed")
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
