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
* EXIT_FILE_NOT_FOUND_ERROR (1): A file does not seem to exist.
* EXIT_CONFIG_ERROR (2): Errors occurred loading or accessing
  configuration e.g. missing configuration parameters, inconsistent
  configuration parameters.
* EXIT_PROCESSING_ERROR (3): Error occurred during processing.

Commands that are submitted to bash are recorded within a
file specified by a cmd_file configuration parameter.

If --dry-run is provided then the commands submitted to bash will not
be executed. This can be useful for seeing what commands will be run
without actually running them.
"""

import errno
import logging
import os
import os.path
import sys
import yaml
from riboviz import logging_utils
from riboviz import params
from riboviz import workflow
from riboviz.tools.demultiplex_fastq import NUM_READS_FILE
from riboviz.utils import value_in_dict
from riboviz.utils import get_fastq_filename
from riboviz.utils import load_deplexed_sample_sheet
from riboviz.utils import get_non_zero_deplexed_samples


EXIT_OK = 0
"""Processing successfully completed. """
EXIT_FILE_NOT_FOUND_ERROR = 1
"""
A file does not seem to exist.
"""
EXIT_CONFIG_ERROR = 2
"""
Errors occurred loading or accessing configuration e.g. missing
configuration parameters, inconsistent configuration parameters.
"""
EXIT_PROCESSING_ERROR = 3
"""
Error occurred during processing.
"""


logging_utils.configure_logging()
LOGGER = logging.getLogger(__name__)
""" Logger """


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

    if value_in_dict(params.NPROCESSES, config):
        nprocesses = int(config[params.NPROCESSES])
    else:
        nprocesses = 1

    is_extract_umis = value_in_dict(params.EXTRACT_UMIS, config)
    if is_trimmed:
        LOGGER.info("Skipping adaptor trimming and barcode/UMI extraction")
        trim_fq = fastq
    else:
        log_file = workflow.get_sample_log_file(
            output_config.logs_dir, sample, "cutadapt", step)
        trim_fq = os.path.join(output_config.tmp_dir, sample + "_trim.fq")
        workflow.cut_adapters(
            config[params.ADAPTERS], fastq, trim_fq, log_file, cmd_config)
        step += 1

        if is_extract_umis:
            extract_trim_fq = os.path.join(
                output_config.tmp_dir, sample + "_extract_trim.fq")
            log_file = workflow.get_sample_log_file(
                output_config.logs_dir, sample, "umi_tools_extract",
                step)
            workflow.extract_barcodes_umis(
                trim_fq, extract_trim_fq, config[params.UMI_REGEXP],
                log_file, cmd_config)
            trim_fq = extract_trim_fq
            step += 1

    non_r_rna_trim_fq = os.path.join(
        output_config.tmp_dir, sample + "_nonrRNA.fq")
    r_rna_map_sam = os.path.join(
        output_config.tmp_dir, sample + "_rRNA_map.sam")
    log_file = workflow.get_sample_log_file(
        output_config.logs_dir, sample, "hisat2_rrna", step)
    workflow.map_to_r_rna(
        trim_fq, r_rna_index, r_rna_map_sam, non_r_rna_trim_fq,
        nprocesses, log_file, cmd_config)
    step += 1

    orf_map_sam = os.path.join(
        output_config.tmp_dir, sample + "_orf_map.sam")
    unaligned_fq = os.path.join(
        output_config.tmp_dir, sample + "_unaligned.fq")
    log_file = workflow.get_sample_log_file(
        output_config.logs_dir, sample, "hisat2_orf", step)
    workflow.map_to_orf(
        non_r_rna_trim_fq, orf_index, orf_map_sam, unaligned_fq,
        nprocesses, log_file, cmd_config)
    step += 1

    orf_map_sam_clean = os.path.join(
        output_config.tmp_dir, sample + "_orf_map_clean.sam")
    log_file = workflow.get_sample_log_file(
        output_config.logs_dir, sample, "trim_5p_mismatch", step)
    workflow.trim_5p_mismatches(
        orf_map_sam, orf_map_sam_clean, py_scripts, log_file,
        cmd_config)
    step += 1

    log_file = workflow.get_sample_log_file(
        output_config.logs_dir, sample, "samtools_view_sort", step)
    sample_out_prefix = os.path.join(output_config.out_dir, sample)
    sample_out_bam = sample_out_prefix + ".bam"
    workflow.sort_bam(
        orf_map_sam_clean, sample_out_bam, nprocesses, log_file,
        cmd_config)
    step += 1

    log_file = workflow.get_sample_log_file(
        output_config.logs_dir, sample, "samtools_index", step)
    workflow.index_bam(sample_out_bam, log_file, cmd_config)
    step += 1

    is_dedup_umis = value_in_dict(params.DEDUP_UMIS, config)
    if is_dedup_umis:
        if not is_extract_umis:
            LOGGER.warning(
                "WARNING: dedup_umis was TRUE but extract_umis was FALSE.")
        is_group_umis = value_in_dict(params.GROUP_UMIS, config)
        if is_group_umis:
            umi_groups = os.path.join(
                output_config.tmp_dir, sample + "_pre_dedup_groups.tsv")
            log_file = workflow.get_sample_log_file(
                output_config.logs_dir, sample, "umi_tools_group",
                step)
            workflow.group_umis(
                sample_out_bam, umi_groups, log_file, cmd_config)
            step += 1

        sample_dedup_bam = sample_out_prefix + "_dedup.bam"
        log_file = workflow.get_sample_log_file(
            output_config.logs_dir, sample, "umi_tools_dedup", step)
        dedup_stats_prefix = os.path.join(
            output_config.tmp_dir, sample + "_dedup_stats")
        workflow.deduplicate_umis(
            sample_out_bam, sample_dedup_bam, dedup_stats_prefix,
            log_file, cmd_config)
        step += 1

        log_file = workflow.get_sample_log_file(
            output_config.logs_dir, sample, "samtools_index", step)
        workflow.index_bam(sample_dedup_bam, log_file, cmd_config)
        sample_out_bam = sample_dedup_bam
        step += 1

        if is_group_umis:
            umi_groups = os.path.join(
                output_config.tmp_dir, sample + "_post_dedup_groups.tsv")
            log_file = workflow.get_sample_log_file(
                output_config.logs_dir, sample, "umi_tools_group",
                step)
            workflow.group_umis(
                sample_out_bam, umi_groups, log_file, cmd_config)
            step += 1

    is_make_bedgraph = value_in_dict(params.MAKE_BEDGRAPH, config)
    if is_make_bedgraph:
        log_file = workflow.get_sample_log_file(
            output_config.logs_dir, sample, "bedtools_genome_cov_plus", step)
        plus_bedgraph = sample_out_prefix + "_plus.bedgraph"
        workflow.make_bedgraph(
            sample_out_bam, plus_bedgraph, True, log_file, cmd_config)
        step += 1

        log_file = workflow.get_sample_log_file(
            output_config.logs_dir, sample, "bedtools_genome_cov_minus", step)
        minus_bedgraph = sample_out_prefix + "_minus.bedgraph"
        workflow.make_bedgraph(
            sample_out_bam, minus_bedgraph, False, log_file, cmd_config)
        step += 1

    orf_gff_file = config[params.ORF_GFF_FILE]
    if not os.path.exists(orf_gff_file):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                orf_gff_file)
    log_file = workflow.get_sample_log_file(
        output_config.logs_dir, sample, "bam_to_h5", step)
    sample_out_h5 = sample_out_prefix + ".h5"
    workflow.bam_to_h5(
        sample_out_bam, sample_out_h5, orf_gff_file, config,
        nprocesses, r_scripts, log_file, cmd_config)
    step += 1

    log_file = workflow.get_sample_log_file(
        output_config.logs_dir, sample, "generate_stats_figs", step)
    workflow.generate_stats_figs(
        sample_out_h5, output_config.out_dir, sample_out_prefix,
        config, nprocesses, r_scripts, log_file, cmd_config)

    LOGGER.info("Finished processing sample: %s", fastq)


def process_samples(samples, in_dir, r_rna_index, orf_index,
                    is_trimmed, config, py_scripts, r_scripts,
                    output_config, cmd_config):
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
                           is_trimmed,
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
    * EXIT_FILE_NOT_FOUND_ERROR (1): A file does not seem to exist.
    * EXIT_CONFIG_ERROR (2): Errors occurred loading or accessing
      configuration e.g. missing configuration parameters,
      inconsistent configuration parameters.
    * EXIT_PROCESSING_ERROR (3): Error occurred during processing.

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
        return EXIT_FILE_NOT_FOUND_ERROR
    except Exception:
        logging.error("Problem reading: %s", config_yaml)
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_CONFIG_ERROR

    if value_in_dict(params.CMD_FILE, config):
        cmd_file = config[params.CMD_FILE]
    else:
        cmd_file = "run_riboviz_vignette.sh"
    LOGGER.info("Command file: %s", cmd_file)
    if os.path.exists(cmd_file):
        os.remove(cmd_file)
    cmd_config = workflow.CmdConfigTuple(cmd_file, is_dry_run)

    try:
        in_dir = config[params.INPUT_DIR]
        output_config = workflow.setup_output_directories(config, cmd_config)
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
        r_rna_index = os.path.join(
            output_config.index_dir, config[params.R_RNA_INDEX_PREFIX])
        orf_index = os.path.join(
            output_config.index_dir, config[params.ORF_INDEX_PREFIX])
        is_build_indices = value_in_dict(params.BUILD_INDICES, config)
        if is_build_indices:
            r_rna_fasta = config[params.R_RNA_FASTA_FILE]
            log_file = os.path.join(
                output_config.logs_dir, "hisat2_build_r_rna.log")
            workflow.build_indices(r_rna_fasta, r_rna_index, log_file,
                                   cmd_config)
            orf_fasta = config[params.ORF_FASTA_FILE]
            log_file = os.path.join(
                output_config.logs_dir, "hisat2_build_orf.log")
            workflow.build_indices(orf_fasta, orf_index, log_file,
                                   cmd_config)
    except KeyError as e:
        logging.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except FileNotFoundError as e:
        logging.error("File not found: %s", e.filename)
        return EXIT_FILE_NOT_FOUND_ERROR
    except Exception:
        logging.error("Problem creating indices")
        exc_type, _, _ = sys.exc_info()
        logging.exception(exc_type.__name__)
        return EXIT_PROCESSING_ERROR

    is_sample_files = value_in_dict(params.FQ_FILES, config)
    is_multiplex_files = value_in_dict(params.MULTIPLEX_FQ_FILES, config)
    is_sample_sheet_file = value_in_dict(params.SAMPLE_SHEET_FILE, config)
    if not is_sample_files and not is_multiplex_files:
        LOGGER.error("No sample files (fq_files) or multiplexed files (multiplex_fq_files) are specified.")
        return EXIT_CONFIG_ERROR
    elif is_sample_files and is_multiplex_files:
        LOGGER.error("Both sample files (fq_files) and multiplexed files (multiplex_fq_files) were specified.")
        return EXIT_CONFIG_ERROR
    elif is_multiplex_files and not is_sample_sheet_file:
        LOGGER.error("Multiplexed files (multiplex_fq_files) are specified but no sample sheet (sample_sheet) file.")
        return EXIT_CONFIG_ERROR

    if is_sample_files:
        samples = config[params.FQ_FILES]
        processed_samples = process_samples(
            samples, in_dir, r_rna_index, orf_index, False, config,
            py_scripts, r_scripts, output_config, cmd_config)
        if not processed_samples:
            return EXIT_PROCESSING_ERROR
        try:
            log_file = os.path.join(
                output_config.logs_dir, "collate_tpms.log")
            workflow.collate_tpms(
                output_config.out_dir, processed_samples, r_scripts,
                log_file, cmd_config)
        except Exception:
            logging.error(("Problem collating TPMs"))
            exc_type, _, _ = sys.exc_info()
            logging.exception(exc_type.__name__)
            return EXIT_PROCESSING_ERROR
    else:
        LOGGER.info("WIP: multiplexed file processing")
        try:
            sample_sheet_file = os.path.join(in_dir,
                                             config[params.SAMPLE_SHEET_FILE])
            if not os.path.exists(sample_sheet_file):
                raise FileNotFoundError(errno.ENOENT,
                                        os.strerror(errno.ENOENT),
                                        sample_sheet_file)
        except FileNotFoundError as e:
            logging.error("File not found: %s", e.filename)
            return EXIT_FILE_NOT_FOUND_ERROR
        except Exception:
            logging.error("Problem processing: %s", sample_sheet_file)
            exc_type, _, _ = sys.exc_info()
            logging.exception(exc_type.__name__)
            return EXIT_PROCESSING_ERROR

        multiplex_files = config[params.MULTIPLEX_FQ_FILES]
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
            workflow.cut_adapters(
                config[params.ADAPTERS], multiplex_file, trim_fq, log_file,
                cmd_config)

            extract_trim_fq = os.path.join(
                output_config.tmp_dir, multiplex_name + "_extract_trim.fq")
            log_file = os.path.join(
                output_config.logs_dir, "umi_tools_extract.log")
            workflow.extract_barcodes_umis(
                trim_fq, extract_trim_fq, config[params.UMI_REGEXP], log_file,
                cmd_config)

            deplex_dir = os.path.join(
                output_config.tmp_dir, multiplex_name + "_deplex")
            log_file = os.path.join(
                output_config.logs_dir, "demultiplex_fastq.log")
            workflow.demultiplex_fastq(extract_trim_fq, sample_sheet_file,
                                       deplex_dir, log_file, cmd_config)

            # TODO this won't work for dry run as it does not exist!
            num_reads_file = os.path.join(deplex_dir, NUM_READS_FILE)
            num_reads = load_deplexed_sample_sheet(num_reads_file)
            samples = get_non_zero_deplexed_samples(num_reads)
            print(samples)  # TODO remove
            if not samples:
                LOGGER.error("No sample files with any reads resulted from demultiplexing.")
                return EXIT_PROCESSING_ERROR

            # Use get_fastq_filename to deduce sample-specific files
            # output by demultiplex_fastq.py - demultiplex_fastq.py
            # uses this function too.
            sample_files = {sample: get_fastq_filename(sample)
                            for sample in samples}
            print(sample_files)  # TODO remove
            processed_samples = process_samples(
                sample_files, deplex_dir, r_rna_index, orf_index,
                True, config, py_scripts, r_scripts, output_config,
                cmd_config)
            if not processed_samples:
                return EXIT_PROCESSING_ERROR
        except FileNotFoundError as e:
            logging.error("File not found: %s", e.filename)
            return EXIT_FILE_NOT_FOUND_ERROR
        except Exception:
            logging.error("Problem processing: %s", multiplex_file)
            exc_type, _, _ = sys.exc_info()
            logging.exception(exc_type.__name__)
            return EXIT_PROCESSING_ERROR

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
