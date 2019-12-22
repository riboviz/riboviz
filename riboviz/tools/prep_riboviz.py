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

The script can process sample data or multiplexed sample data
(relevant configuration parameters are shown in brackets).

Process ribosome profiling data:

* Reads configuration information from YAML configuration file.
* Builds hisat2 indices if requested ("build_indices: TRUE") using
  "hisat2 build" and saves these into an index directory
  ("dir_index").
* Processes each sample fastq[.gz] file (sample IDs and files are
  listed in "fq_files" and are assumed to be relative to "dir_in") in
  turn:
  - Cuts out sequencing library adapters ("adapters", default
    "CTGTAGGCACC") using "cutadapt".
  - Extracts UMIs using "umi_tools extract", if requested
    ("extract_umis: TRUE"), using a UMI-tools-compliant
    regular expression pattern ("umi_regexp").
  - Removes rRNA or other contaminating reads by alignment to
    rRNA index files ("rrna_index_prefix") using "hisat2".
  - Aligns remaining reads to ORFs index files
    ("orf_index_prefix"). using "hisat2".
  - Trims 5' mismatches from reads and remove reads with more than 2
    mismatches using "trim_5p_mismatch.py".
  - Outputs UMI groups pre-deduplication using "umi_tools group" if
    requested ("dedup_umis: TRUE" and "group_umis: TRUE")
  - Deduplicates UMIs using "umi_tools dedup", if requested
    ("dedup_umis: TRUE")
  - Outputs UMI groups post-deduplication using "umi_tools group" if
    requested ("dedup_umis: TRUE" and "group_umis:TRUE")
  - Exports bedgraph files for plus and minus strands, if requested
    ("make_bedgraph: TRUE") using "bedtools genomecov".
  - Writes intermediate files produced above into a sample-specific
    directory under the temporary directory  ("dir_tmp").
  - Makes length-sensitive alignments in compressed h5 format using
    "bam_to_h5.R".
  - Generates summary statistics, and analyses and QC plots for both
    RPF and mRNA datasets using "generate_stats_figs.R". This
    includes estimated read counts, reads per base, and transcripts
    per million for each ORF in each sample.
  - Writes output files produced above into an sample-specific
    directory under the output directory ("dir_out").
* Collates TPMs across all processed fastq[.gz] files, using
  "collate_tpms.R" and writes into output directory ("dir_out").

Process multiplexed ribosome profiling data:

* Reads configuration information from YAML configuration file.
* Builds hisat2 indices if requested ("build_indices: TRUE") using
  "hisat2 build" and saves these into an index directory
  ("dir_index").
* Reads fastq[.gz] file (the file is listed in "multiplex_fq_files"
  and is assumed to be relative to "dir_in").
* Cuts out sequencing library adapters ("adapters", default
  "CTGTAGGCACC") using "cutadapt".
* Extracts barcodes and UMIs using "umi_tools extract", if requested
  ("extract_umis: TRUE"), using a UMI-tools-compliant
  regular expression pattern ("umi_regexp").
* Demultiplexes fastq[.gz] file with reference to a sample sheet
  ("sample_sheet"), using "demultiplex_fastq.py".
* Processes each demultiplexed fastq[.gz], which has one or more
  reads, in turn:
  - Removes rRNA or other contaminating reads by alignment to rRNA
    index files ("rrna_index_prefix") using "hisat2".
  - Aligns remaining reads to ORFs index files ("orf_index_prefix")
    using "hisat2".
  - Trims 5' mismatches from reads and remove reads with more than 2
    mismatches using "trim_5p_mismatch.py".
  - Outputs UMI groups pre-deduplication using "umi_tools group" if
    requested ("dedup_umis: TRUE" and "group_umis: TRUE").
  - Deduplicates UMIs using "umi_tools dedup", if requested
    ("dedup_umis: TRUE").
  - Outputs UMI groups post-deduplication using "umi_tools group" if
    requested ("dedup_umis: TRUE" and "group_umis: TRUE")
  - Exports bedgraph files for plus and minus strands, if requested
    ("make_bedgraph: TRUE") using "bedtools genomecov".
  - Writes intermediate files produced above into a sample-specific
    directory under the temporary directory  ("dir_tmp").
  - Makes length-sensitive alignments in compressed h5 format using
    "bam_to_h5.R".
  - Generates summary statistics, and analyses and QC plots for both
    RPF and mRNA datasets using "generate_stats_figs.R". This
    includes estimated read counts, reads per base, and transcripts
    per million for each ORF in each sample.
  - Writes output files produced above into an sample-specific
    directory under the output directory ("dir_out").
* Collates TPMs across all demultiplexed fastq[.gz] files, using
  "collate_tpms.R" and writes into output directory ("dir_out").

The script can parallelize parts of its operation over many
processes ("num_processes"):

* This value is used to configure "hisat2", "samtools sort",
  "bam_to_h5.R" and "generate_stats_figs.R".
* For "cutadapt", the number of available processors on the host will
  be used.

`prep_riboviz.py` returns the following exit codes:

* 0: Processing successfully completed.
* 1: A file does not seem to exist.
* 2: Errors occurred loading or accessing configuration e.g. missing configuration parameters, inconsistent configuration parameters.
* 3: Error occurred during processing.

Commands that are submitted to bash are recorded within a file
specified by a "cmd_file" configuration parameter.

If "--dry-run" is provided then the commands submitted to bash will
not be executed. This can be useful for both seeing what commands will
be run, and validating the configuration, without actually running
the commands.
"""

from datetime import datetime
import errno
import logging
import os
import os.path
import sys
import yaml
from riboviz import logging_utils
from riboviz import params
from riboviz import provenance
from riboviz import workflow
from riboviz.tools.demultiplex_fastq import NUM_READS_FILE
from riboviz.utils import SAMPLE_ID
from riboviz.utils import value_in_dict
from riboviz.utils import get_fastq_filename
from riboviz.utils import load_sample_sheet
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


def process_sample(sample, fastq, r_rna_index, orf_index,
                   is_trimmed, config, tmp_dir, out_dir, run_config):
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
    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    :param out_dir: Output directory
    :type out_dir: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if fastq, other files or a third-party
    tool cannot be found
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    :raise KeyError: if config is missing required configuration
    """
    LOGGER.info("Processing sample: %s", sample)
    step = 1
    LOGGER.info("Processing file: %s", fastq)
    is_extract_umis = value_in_dict(params.EXTRACT_UMIS, config)
    if is_trimmed:
        LOGGER.info("Skipping adaptor trimming and barcode/UMI extraction")
        trim_fq = fastq
    else:
        log_file = workflow.get_log_file(
            run_config.logs_dir, "cutadapt", step)
        trim_fq = os.path.join(tmp_dir, "trim.fq")
        workflow.cut_adapters(
            config[params.ADAPTERS], fastq, trim_fq, log_file, run_config)
        step += 1

        if is_extract_umis:
            extract_trim_fq = os.path.join(tmp_dir, "extract_trim.fq")
            log_file = workflow.get_log_file(
                run_config.logs_dir, "umi_tools_extract", step)
            workflow.extract_barcodes_umis(
                trim_fq, extract_trim_fq, config[params.UMI_REGEXP],
                log_file, run_config)
            trim_fq = extract_trim_fq
            step += 1

    non_r_rna_trim_fq = os.path.join(tmp_dir, "nonrRNA.fq")
    r_rna_map_sam = os.path.join(tmp_dir, "rRNA_map.sam")
    log_file = workflow.get_log_file(
        run_config.logs_dir, "hisat2_rrna", step)
    workflow.map_to_r_rna(
        trim_fq, r_rna_index, r_rna_map_sam, non_r_rna_trim_fq,
        log_file, run_config)
    step += 1

    orf_map_sam = os.path.join(tmp_dir, "orf_map.sam")
    unaligned_fq = os.path.join(tmp_dir, "unaligned.fq")
    log_file = workflow.get_log_file(
        run_config.logs_dir, "hisat2_orf", step)
    workflow.map_to_orf(
        non_r_rna_trim_fq, orf_index, orf_map_sam, unaligned_fq,
        log_file, run_config)
    step += 1

    orf_map_sam_clean = os.path.join(tmp_dir, "orf_map_clean.sam")
    log_file = workflow.get_log_file(
        run_config.logs_dir, "trim_5p_mismatch", step)
    workflow.trim_5p_mismatches(
        orf_map_sam, orf_map_sam_clean, log_file, run_config)
    step += 1

    log_file = workflow.get_log_file(
        run_config.logs_dir, "samtools_view_sort", step)
    sample_out_prefix = os.path.join(out_dir, sample)

    is_dedup_umis = value_in_dict(params.DEDUP_UMIS, config)
    if is_dedup_umis:
        # Create BAM file in temporary directory
        sample_bam = os.path.join(tmp_dir, "pre_dedup.bam")
    else:
        sample_bam = sample_out_prefix + ".bam"
        sample_out_bam = sample_bam

    workflow.sort_bam(
        orf_map_sam_clean, sample_bam, log_file, run_config)
    step += 1
    log_file = workflow.get_log_file(
        run_config.logs_dir, "samtools_index", step)
    workflow.index_bam(sample_bam, log_file, run_config)
    step += 1

    if is_dedup_umis:
        if not is_extract_umis:
            LOGGER.warning(
                "WARNING: dedup_umis was TRUE but extract_umis was FALSE.")
        is_group_umis = value_in_dict(params.GROUP_UMIS, config)
        if is_group_umis:
            umi_groups = os.path.join(tmp_dir, "pre_dedup_groups.tsv")
            log_file = workflow.get_log_file(
                run_config.logs_dir, "umi_tools_group", step)
            workflow.group_umis(
                sample_bam, umi_groups, log_file, run_config)
            step += 1

        sample_out_bam = sample_out_prefix + ".bam"
        log_file = workflow.get_log_file(
            run_config.logs_dir, "umi_tools_dedup", step)
        dedup_stats_prefix = os.path.join(tmp_dir, "dedup_stats")
        workflow.deduplicate_umis(
            sample_bam, sample_out_bam, dedup_stats_prefix,
            log_file, run_config)
        step += 1

        log_file = workflow.get_log_file(
            run_config.logs_dir, "samtools_index", step)
        workflow.index_bam(sample_out_bam, log_file, run_config)
        step += 1

        if is_group_umis:
            umi_groups = os.path.join(tmp_dir, "post_dedup_groups.tsv")
            log_file = workflow.get_log_file(
                run_config.logs_dir, "umi_tools_group", step)
            workflow.group_umis(
                sample_out_bam, umi_groups, log_file, run_config)
            step += 1

    is_make_bedgraph = value_in_dict(params.MAKE_BEDGRAPH, config)
    if is_make_bedgraph:
        log_file = workflow.get_log_file(
            run_config.logs_dir, "bedtools_genome_cov_plus", step)
        plus_bedgraph = os.path.join(out_dir, "plus.bedgraph")
        workflow.make_bedgraph(
            sample_out_bam, plus_bedgraph, True, log_file, run_config)
        step += 1

        log_file = workflow.get_log_file(
            run_config.logs_dir, "bedtools_genome_cov_minus", step)
        minus_bedgraph = os.path.join(out_dir, "minus.bedgraph")
        workflow.make_bedgraph(
            sample_out_bam, minus_bedgraph, False, log_file, run_config)
        step += 1

    orf_gff_file = config[params.ORF_GFF_FILE]
    log_file = workflow.get_log_file(
        run_config.logs_dir, "bam_to_h5", step)
    sample_out_h5 = sample_out_prefix + ".h5"
    workflow.bam_to_h5(sample_out_bam, sample_out_h5, orf_gff_file,
                       config, log_file, run_config)
    step += 1


    log_file = workflow.get_log_file(
        run_config.logs_dir, "generate_stats_figs", step)
    workflow.generate_stats_figs(
        sample_out_h5, out_dir, config, log_file, run_config)

    LOGGER.info("Finished processing sample: %s", fastq)


def process_samples(samples, in_dir, r_rna_index, orf_index,
                    is_trimmed, config, tmp_dir, out_dir, run_config,
                    check_samples_exist=True):
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
    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    :param out_dir: Output directory
    :type out_dir: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :param check_samples_exist: If run_config.is_dry_run then should a
    check be made for the existence of sample files?
    :type check_samples_exist: bool
    :type are_trimmed: bool
    :return: names of successfully-processed samples
    :rtype: list(str or unicode)
    """
    LOGGER.info("Processing samples")
    successes = []
    num_samples = len(samples)
    for sample in list(samples.keys()):
        try:
            fastq = os.path.join(in_dir, samples[sample])
            if check_samples_exist:
                if not os.path.exists(fastq):
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), fastq)
            sample_tmp_dir = os.path.join(tmp_dir, sample)
            sample_out_dir = os.path.join(out_dir, sample)
            sample_logs_dir = os.path.join(run_config.logs_dir, sample)
            sample_run_config = workflow.RunConfigTuple(
                run_config.py_scripts,
                run_config.r_scripts,
                run_config.cmd_file,
                run_config.is_dry_run,
                sample_logs_dir,
                run_config.nprocesses)
            for directory in [sample_tmp_dir, sample_out_dir, sample_logs_dir]:
                workflow.create_directory(directory,
                                          run_config.cmd_file,
                                          run_config.is_dry_run)
            process_sample(
                sample, fastq, r_rna_index, orf_index, is_trimmed,
                config, sample_tmp_dir, sample_out_dir, sample_run_config)
            successes.append(sample)
        except FileNotFoundError as e:
            LOGGER.error("File not found: %s", e.filename)
        except Exception:
            LOGGER.error("Problem processing sample: %s", sample)
            exc_type, _, _ = sys.exc_info()
            LOGGER.exception(exc_type.__name__)
    num_failed = num_samples - len(successes)
    LOGGER.info("Finished processing %d samples, %d failed",
                num_samples, num_failed)
    return successes


def run_workflow(py_scripts, r_scripts, config_yaml, is_dry_run=False):
    """
    Run the RiboViz workflow.

    :param py_scripts: Python scripts directory
    :type py_scripts: str or unicode
    :param r_scripts: R scripts directory
    :type r_scripts: str or unicode
    :param config_yaml: YAML configuration file path
    :type config_yaml: str or unicode
    :param is_dry_run: Don't execute workflow commands (useful for
    seeing what commands would be executed)
    :type is_dry_run: bool
    :raise FileNotFoundError: if a file cannot be found
    :raise KeyError: if a configuration parameter is missing
    :raise ValueError: if a configuration parameter has an invalid
    value
    :raise TypeError: if a configuration parameter has an invalid type
    :raise Exception: if any other error arises
    """
    LOGGER.info("Running under Python: %s", sys.version)
    LOGGER.info("Configuration file: %s", config_yaml)

    with open(config_yaml, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)

    if value_in_dict(params.CMD_FILE, config):
        cmd_file = config[params.CMD_FILE]
    else:
        cmd_file = "run_riboviz_vignette.sh"
    if os.path.exists(cmd_file):
        os.remove(cmd_file)
    LOGGER.info("Command file: %s", cmd_file)

    if value_in_dict(params.NUM_PROCESSES, config):
        nprocesses = int(config[params.NUM_PROCESSES])
    else:
        nprocesses = 1
    LOGGER.info("Number of processes: %d", nprocesses)

    index_dir = config[params.INDEX_DIR]
    tmp_dir = config[params.TMP_DIR]
    out_dir = config[params.OUTPUT_DIR]
    base_logs_dir = config[params.LOGS_DIR]
    logs_dir = os.path.join(
        base_logs_dir, datetime.now().strftime('%Y%m%d-%H%M%S'))
    dirs = [index_dir, tmp_dir, out_dir, logs_dir]
    for directory in dirs:
        workflow.create_directory(directory, cmd_file, is_dry_run)

    run_config = workflow.RunConfigTuple(
        py_scripts, r_scripts, cmd_file, is_dry_run, logs_dir,
        nprocesses)

    in_dir = config[params.INPUT_DIR]
    LOGGER.info("Build indices for alignment, if necessary/requested")
    r_rna_index = os.path.join(
        index_dir, config[params.RRNA_INDEX_PREFIX])
    orf_index = os.path.join(
        index_dir, config[params.ORF_INDEX_PREFIX])
    is_build_indices = value_in_dict(params.BUILD_INDICES, config)
    if is_build_indices:
        r_rna_fasta = config[params.RRNA_FASTA_FILE]
        log_file = os.path.join(logs_dir, "hisat2_build_r_rna.log")
        workflow.build_indices(r_rna_fasta, r_rna_index, log_file,
                               run_config)
        orf_fasta = config[params.ORF_FASTA_FILE]
        log_file = os.path.join(logs_dir, "hisat2_build_orf.log")
        workflow.build_indices(orf_fasta, orf_index, log_file,
                               run_config)

    is_sample_files = value_in_dict(params.FQ_FILES, config)
    is_multiplex_files = value_in_dict(params.MULTIPLEX_FQ_FILES, config)
    is_sample_sheet_file = value_in_dict(params.SAMPLE_SHEET, config)
    if not is_sample_files and not is_multiplex_files:
        raise ValueError(
            "No sample files ({}) or multiplexed files ({}) are specified".format(
                params.FQ_FILES, params.MULTIPLEX_FQ_FILES))
    elif is_sample_files and is_multiplex_files:
        raise ValueError(
            "Both sample files ({}) and multiplexed files ({}) are specified".format(
                params.FQ_FILES, params.MULTIPLEX_FQ_FILES))
    elif is_multiplex_files and not is_sample_sheet_file:
        raise ValueError(
            "Multiplexed files ({}) are specified but no sample sheet ({})".format(
                params.MULTIPLEX_FQ_FILES, params.SAMPLE_SHEET))

    if is_sample_files:
        samples = config[params.FQ_FILES]
        processed_samples = process_samples(
            samples, in_dir, r_rna_index, orf_index,
            False, config, tmp_dir, out_dir, run_config)
        if not processed_samples:
            raise Exception("No samples were processed successfully")

    else:
        sample_sheet_file = os.path.join(in_dir,
                                         config[params.SAMPLE_SHEET])
        if not os.path.exists(sample_sheet_file):
            raise FileNotFoundError(errno.ENOENT,
                                    os.strerror(errno.ENOENT),
                                    sample_sheet_file)
        multiplex_files = config[params.MULTIPLEX_FQ_FILES]
        multiplex_file = multiplex_files[0]
        multiplex_name, ext = os.path.splitext(
            os.path.basename(multiplex_file))
        if ext.lower() == ".gz":
            # Trim off .fastq
            multiplex_name = os.path.splitext(os.path.basename(multiplex_name))[0]
        multiplex_file = os.path.join(in_dir, multiplex_file)
        LOGGER.info("Processing file: %s", multiplex_file)
        if not os.path.exists(multiplex_file):
            raise FileNotFoundError(errno.ENOENT,
                                    os.strerror(errno.ENOENT),
                                    multiplex_file)
        trim_fq = os.path.join(
            tmp_dir, multiplex_name + "_trim.fq")
        log_file = os.path.join(logs_dir, "cutadapt.log")
        workflow.cut_adapters(
            config[params.ADAPTERS], multiplex_file, trim_fq, log_file,
            run_config)

        extract_trim_fq = os.path.join(
            tmp_dir, multiplex_name + "_extract_trim.fq")
        log_file = os.path.join(logs_dir, "umi_tools_extract.log")
        workflow.extract_barcodes_umis(
            trim_fq, extract_trim_fq, config[params.UMI_REGEXP], log_file,
            run_config)

        deplex_dir = os.path.join(tmp_dir, multiplex_name + "_deplex")
        log_file = os.path.join(logs_dir, "demultiplex_fastq.log")
        workflow.demultiplex_fastq(extract_trim_fq, sample_sheet_file,
                                   deplex_dir, log_file, run_config)

        if not is_dry_run:
            num_reads_file = os.path.join(deplex_dir, NUM_READS_FILE)
            num_reads = load_deplexed_sample_sheet(num_reads_file)
            samples = get_non_zero_deplexed_samples(num_reads)
            if not samples:
                raise Exception(
                    "No non-empty sample files were produced by demultiplexing")
        else:
            # If doing a dry run, use the sample_sheet_file to get
            # the names of the samples.
            sample_sheet = load_sample_sheet(sample_sheet_file)
            samples = list(sample_sheet[SAMPLE_ID])
            if not samples:
                raise Exception(
                    "No samples are specified in {}".format(sample_sheet_file))
        # Use get_fastq_filename to deduce sample-specific files
        # output by demultiplex_fastq.py - demultiplex_fastq.py
        # uses this function too.
        sample_files = {sample: get_fastq_filename(sample)
                        for sample in samples}
        processed_samples = process_samples(
            sample_files, deplex_dir, r_rna_index, orf_index,
            True, config, tmp_dir, out_dir, run_config, False)
        if not processed_samples:
            raise Exception("No samples were processed successfully")

    log_file = os.path.join(logs_dir, "collate_tpms.log")
    workflow.collate_tpms(
        out_dir, processed_samples, True, log_file, run_config)

    LOGGER.info("Completed")


def prep_riboviz(py_scripts, r_scripts, config_yaml, is_dry_run=False):
    """
    Run the RiboViz workflow. This function invokes run_workflow,
    catchesand logs any exceptions and maps these to exit codes.

    Exit codes are as follows:

    * EXIT_OK (0): Processing successfully completed.
    * EXIT_FILE_NOT_FOUND_ERROR (1): A file does not seem to exist.
    * EXIT_CONFIG_ERROR (2): Errors occurred loading or accessing
      configuration e.g. missing configuration parameters,
      inconsistent configuration parameters.e
    * EXIT_PROCESSING_ERROR (3): Error occurred during processing.

    :param py_scripts: Python scripts directory
    :type py_scripts: str or unicode
    :param r_scripts: R scripts direectory
    :type r_scripts: str or unicode
    :param config_yaml: YAML configuration file path
    :type config_yaml: str or unicodee
    :param is_dry_run: Don't execute weorkflow commands (useful for
    seeing what commands would be execueted)
    :type is_dry_run: boole
    :return: exit code
    :rtype: int
    """
    LOGGER.info(provenance.get_version(__file__))
    try:
        run_workflow(py_scripts, r_scripts, config_yaml, is_dry_run)
    except FileNotFoundError as e:
        LOGGER.error("File not found: %s", e.filename)
        return EXIT_FILE_NOT_FOUND_ERROR
    except KeyError as e:
        LOGGER.error("Missing configuration parameter: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except ValueError as e:
        LOGGER.error("Configuration parameter error: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except TypeError as e:
        LOGGER.error("Configuration parameter error: %s", e.args[0])
        return EXIT_CONFIG_ERROR
    except Exception:
        LOGGER.error("Processing error: %s", config_yaml)
        exc_type, _, _ = sys.exc_info()
        LOGGER.exception(exc_type.__name__)
        return EXIT_PROCESSING_ERROR
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
