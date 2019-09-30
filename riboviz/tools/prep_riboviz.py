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
    :raise OSError: if hisat2-build cannot be found (Python 2)
    :raise FileNotFoundError: if hisat2-build cannot be found (Python 3)
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
    :raise IOError: if fastq cannot be found (Python 2)
    :raise OSError: if a third-party tool cannot be found (Python 2)
    :raise FileNotFoundError: if fastq or a third-party tool cannot be
    found (Python 3)
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    :raise KeyError: if config is missing required configuration
    """
    LOGGER.info("Processing sample: %s", sample)
    step = 1

    if not os.path.exists(fastq):
        raise IOError(errno.ENOENT,
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
    py_major = sys.version_info.major
    if py_major == 3:
        # cutadapt and Python 3 allows all available processors to
        # be requested.
        cmd += ["-j", str(0)]
    process_utils.run_logged_command(cmd, log_file, cmd_file, dry_run)
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
           "--orf_gff_file=" + config["orf_gff_file"],
           "--dir_out=" + out_dir,
           "--t_rna=" + config["t_rna"],
           "--codon_pos=" + config["codon_pos"],
           "--features_file=" + config["features_file"],
           "--do_pos_sp_nt_freq=" + str(config["do_pos_sp_nt_freq"])]
    if "count_threshold" in config and config["count_threshold"] is not None:
        cmd.append("--count_threshold=" + str(config["count_threshold"]))
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
    :raise OSError: if Rscript cannot be found (Python 2)
    :raise FileNotFoundError: if Rscript cannot be found (Python 3)
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
    except IOError as e:
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
        r_rna_index = os.path.join(index_dir,
                                   config["rRNA_index"])
        orf_index = os.path.join(index_dir,
                                 config["orf_index"])
        if config["build_indices"]:
            build_indices(config["rRNA_fasta"],
                          r_rna_index,
                          "r_rna",
                          logs_dir,
                          cmd_file,
                          dry_run)
            build_indices(config["orf_fasta"],
                          orf_index,
                          "orf",
                          logs_dir,
                          cmd_file,
                          dry_run)
    except KeyError as e:
        logging.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
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
        except IOError as e:
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
