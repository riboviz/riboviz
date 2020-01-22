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
* 2: Errors occurred loading or accessing configuration e.g. missing
  configuration parameters, inconsistent configuration parameters.
* 3: Error occurred during processing.

Commands that are submitted to bash are recorded within a file
specified by a "cmd_file" configuration parameter.

Information on the tools used, the files read and written are recorded
within a file specified by a "workflow_record_file" configuration
parameter.

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
from riboviz import demultiplex_fastq
from riboviz import fastq
from riboviz import logging_utils
from riboviz import params
from riboviz import provenance
from riboviz import sample_sheets
from riboviz import workflow
from riboviz import workflow_record
from riboviz.trim_5p_mismatch import TRIM_5P_MISMATCH_FILE
from riboviz.utils import value_in_dict


DEFAULT_CMD_FILE = "run_riboviz_vignette.sh"
""" Default command file """
DEFAULT_WORKFLOW_RECORD_FILE = "workflow_record.tsv"
""" Default workflow record file """

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

DRY_RUN = "--dry-run"
""" Command-line flag for dry run """

EXIT_OK = 0
"""Processing successfully completed. """
EXIT_FILE_NOT_FOUND_ERROR = 1
""" A file does not seem to exist. """
EXIT_CONFIG_ERROR = 2
"""
Errors occurred loading or accessing configuration e.g. missing
configuration parameters, inconsistent configuration parameters.
"""
EXIT_PROCESSING_ERROR = 3
""" Error occurred during processing. """

LOG_FORMAT = "{:02d}_{}"
""" File name format for step-specific log files. """

logging_utils.configure_logging()
LOGGER = logging.getLogger(__name__)
""" Logger """


def process_sample(sample, sample_fastq, index_dir, r_rna_index,
                   orf_index, is_trimmed, config, tmp_dir, out_dir,
                   run_config):
    """
    Process a single FASTQ sample file.

    :param sample: Sample name
    :type sample: str or unicode
    :param sample_fastq: Sample FASTQ file
    :type sample_fastq: str or unicode
    :param index_dir: Index directory
    :type index_dir: str or unicode
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
    :raise FileNotFoundError: if sample_fastq, other files or a third-party
    tool cannot be found
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    :raise KeyError: if config is missing required configuration
    """
    LOGGER.info("Processing sample: %s", sample)
    step = 1
    LOGGER.info("Processing file: %s", sample_fastq)
    is_extract_umis = value_in_dict(params.EXTRACT_UMIS, config)
    if is_trimmed:
        LOGGER.info("Skipping adaptor trimming and barcode/UMI extraction")
        trim_fq = sample_fastq
    else:
        log_file = os.path.join(run_config.logs_dir,
                                LOG_FORMAT.format(step, "cutadapt.log"))
        trim_fq = os.path.join(tmp_dir, ADAPTER_TRIM_FQ)
        workflow.cut_adapters(sample, config[params.ADAPTERS],
                              sample_fastq, trim_fq, log_file,
                              run_config)
        step += 1

        if is_extract_umis:
            extract_trim_fq = os.path.join(tmp_dir, UMI_EXTRACT_FQ)
            log_file = os.path.join(
                run_config.logs_dir,
                LOG_FORMAT.format(step, "umi_tools_extract.log"))
            workflow.extract_barcodes_umis(sample,
                                           trim_fq,
                                           extract_trim_fq,
                                           config[params.UMI_REGEXP],
                                           log_file,
                                           run_config)
            trim_fq = extract_trim_fq
            step += 1

    non_r_rna_trim_fq = os.path.join(tmp_dir, NON_RRNA_FQ)
    r_rna_map_sam = os.path.join(tmp_dir, RRNA_MAP_SAM)
    log_file = os.path.join(run_config.logs_dir,
                            LOG_FORMAT.format(step, "hisat2_rrna.log"))
    workflow.map_to_r_rna(sample, trim_fq, index_dir, r_rna_index,
                          r_rna_map_sam, non_r_rna_trim_fq, log_file,
                          run_config)
    step += 1

    orf_map_sam = os.path.join(tmp_dir, ORF_MAP_SAM)
    unaligned_fq = os.path.join(tmp_dir, UNALIGNED_FQ)
    log_file = os.path.join(run_config.logs_dir,
                            LOG_FORMAT.format(step, "hisat2_orf.log"))
    workflow.map_to_orf(sample, non_r_rna_trim_fq, index_dir,
                        orf_index, orf_map_sam, unaligned_fq,
                        log_file, run_config)
    step += 1

    orf_map_sam_clean = os.path.join(tmp_dir, ORF_MAP_CLEAN_SAM)
    trim_5p_mismatch_tsv = os.path.join(tmp_dir,
                                        TRIM_5P_MISMATCH_FILE)
    log_file = os.path.join(run_config.logs_dir,
                            LOG_FORMAT.format(step, "trim_5p_mismatch.log"))
    workflow.trim_5p_mismatches(sample, orf_map_sam,
                                orf_map_sam_clean,
                                trim_5p_mismatch_tsv, log_file,
                                run_config)
    step += 1

    log_file = os.path.join(run_config.logs_dir,
                            LOG_FORMAT.format(step, "samtools_view_sort.log"))
    sample_out_prefix = os.path.join(out_dir, sample)

    is_dedup_umis = value_in_dict(params.DEDUP_UMIS, config)
    if is_dedup_umis:
        # Create BAM file in temporary directory
        sample_bam = os.path.join(tmp_dir, PRE_DEDUP_BAM)
    else:
        sample_bam = BAM_FORMAT.format(sample_out_prefix)
        sample_out_bam = sample_bam

    workflow.sort_bam(sample, orf_map_sam_clean, sample_bam, log_file,
                      run_config)
    step += 1
    log_file = os.path.join(run_config.logs_dir,
                            LOG_FORMAT.format(step, "samtools_index.log"))
    workflow.index_bam(sample, sample_bam, log_file, run_config)
    step += 1

    if is_dedup_umis:
        if not is_extract_umis:
            LOGGER.warning(
                "WARNING: dedup_umis was TRUE but extract_umis was FALSE.")
        is_group_umis = value_in_dict(params.GROUP_UMIS, config)
        if is_group_umis:
            umi_groups = os.path.join(tmp_dir, PRE_DEDUP_GROUPS_TSV)
            log_file = os.path.join(
                run_config.logs_dir,
                LOG_FORMAT.format(step, "umi_tools_group.log"))
            workflow.group_umis(sample, sample_bam, umi_groups,
                                log_file, run_config)
            step += 1

        sample_out_bam = BAM_FORMAT.format(sample_out_prefix)
        log_file = os.path.join(
            run_config.logs_dir,
            LOG_FORMAT.format(step, "umi_tools_dedup.log"))
        dedup_stats_prefix = os.path.join(tmp_dir, DEDUP_STATS_PREFIX)
        workflow.deduplicate_umis(
            sample, sample_bam, sample_out_bam, dedup_stats_prefix,
            log_file, run_config)
        step += 1

        log_file = os.path.join(
            run_config.logs_dir,
            LOG_FORMAT.format(step, "samtools_index.log"))
        workflow.index_bam(sample, sample_out_bam, log_file,
                           run_config)
        step += 1

        if is_group_umis:
            umi_groups = os.path.join(tmp_dir, POST_DEDUP_GROUPS_TSV)
            log_file = os.path.join(
                run_config.logs_dir,
                LOG_FORMAT.format(step, "umi_tools_group.log"))
            workflow.group_umis(sample, sample_out_bam, umi_groups,
                                log_file, run_config)
            step += 1

    is_make_bedgraph = value_in_dict(params.MAKE_BEDGRAPH, config)
    if is_make_bedgraph:
        log_file = os.path.join(
            run_config.logs_dir,
            LOG_FORMAT.format(step, "bedtools_genome_cov_plus.log"))
        plus_bedgraph = os.path.join(out_dir, PLUS_BEDGRAPH)
        workflow.make_bedgraph(sample, sample_out_bam, plus_bedgraph,
                               True, log_file, run_config)
        step += 1

        log_file = os.path.join(
            run_config.logs_dir,
            LOG_FORMAT.format(step, "bedtools_genome_cov_minus.log"))
        minus_bedgraph = os.path.join(out_dir, MINUS_BEDGRAPH)
        workflow.make_bedgraph(sample, sample_out_bam, minus_bedgraph,
                               False, log_file, run_config)
        step += 1

    orf_gff_file = config[params.ORF_GFF_FILE]
    log_file = os.path.join(run_config.logs_dir,
                            LOG_FORMAT.format(step, "bam_to_h5.log"))
    sample_out_h5 = H5_FORMAT.format(sample_out_prefix)
    workflow.bam_to_h5(sample, sample_out_bam, sample_out_h5,
                       orf_gff_file, config, log_file, run_config)
    step += 1

    log_file = os.path.join(
        run_config.logs_dir,
        LOG_FORMAT.format(step, "generate_stats_figs.log"))
    workflow.generate_stats_figs(sample, sample_out_h5, out_dir,
                                 config, log_file, run_config)

    LOGGER.info("Finished processing sample: %s", sample_fastq)


def process_samples(samples, in_dir, index_dir, r_rna_index,
                    orf_index, is_trimmed, config, tmp_dir, out_dir,
                    run_config, check_samples_exist=True):
    """
    Process FASTQ sample files. Any exceptions in the processing of
    any sample are logged but are not thrown from this function.

    :param samples: Sample names and files
    :type samples: dict
    :param in_dir: Directory with sample files
    :type in_dir: str or unicode
    :param index_dir: Index directory
    :type index_dir: str or unicode
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
            sample_fastq = os.path.join(in_dir, samples[sample])
            if check_samples_exist:
                if not os.path.exists(sample_fastq):
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), sample_fastq)
            sample_tmp_dir = os.path.join(tmp_dir, sample)
            sample_out_dir = os.path.join(out_dir, sample)
            sample_logs_dir = os.path.join(run_config.logs_dir, sample)
            sample_run_config = workflow.RunConfigTuple(
                run_config.r_scripts,
                run_config.cmd_file,
                run_config.workflow_record_file,
                run_config.is_dry_run,
                sample_logs_dir,
                run_config.nprocesses)
            for directory in [sample_tmp_dir, sample_out_dir, sample_logs_dir]:
                workflow.create_directory(directory,
                                          run_config.cmd_file,
                                          run_config.is_dry_run)
            process_sample(sample, sample_fastq, index_dir,
                           r_rna_index, orf_index, is_trimmed, config,
                           sample_tmp_dir, sample_out_dir,
                           sample_run_config)
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


def run_workflow(r_scripts, config_yaml, is_dry_run=False):
    """
    Run the RiboViz workflow.

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
        cmd_file = DEFAULT_CMD_FILE
    if os.path.exists(cmd_file):
        os.remove(cmd_file)
    LOGGER.info("Command file: %s", cmd_file)

    if value_in_dict(params.WORKFLOW_RECORD_FILE, config):
        workflow_record_file = config[params.WORKFLOW_RECORD_FILE]
    else:
        workflow_record_file = DEFAULT_WORKFLOW_RECORD_FILE
    if not is_dry_run:
        if os.path.exists(workflow_record_file):
            os.remove(workflow_record_file)
        LOGGER.info("Workflow record file: %s", workflow_record_file)
        workflow_record.create_record_file(workflow_record_file)

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
        r_scripts,
        cmd_file,
        workflow_record_file,
        is_dry_run,
        logs_dir,
        nprocesses)

    in_dir = config[params.INPUT_DIR]
    LOGGER.info("Build indices for alignment, if necessary/requested")
    r_rna_index = config[params.RRNA_INDEX_PREFIX]
    orf_index = config[params.ORF_INDEX_PREFIX]
    is_build_indices = value_in_dict(params.BUILD_INDICES, config)
    if is_build_indices:
        r_rna_fasta = config[params.RRNA_FASTA_FILE]
        log_file = os.path.join(logs_dir, "hisat2_build_r_rna.log")
        workflow.build_indices(r_rna_fasta, index_dir, r_rna_index,
                               log_file, run_config)
        orf_fasta = config[params.ORF_FASTA_FILE]
        log_file = os.path.join(logs_dir, "hisat2_build_orf.log")
        workflow.build_indices(orf_fasta, index_dir, orf_index,
                               log_file, run_config)

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
            samples, in_dir, index_dir, r_rna_index, orf_index,
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
        multiplex_name = os.path.splitext(
            os.path.basename(fastq.strip_fastq_gz(multiplex_file)))[0]
        multiplex_file = os.path.join(in_dir, multiplex_file)
        LOGGER.info("Processing file: %s", multiplex_file)
        if not os.path.exists(multiplex_file):
            raise FileNotFoundError(errno.ENOENT,
                                    os.strerror(errno.ENOENT),
                                    multiplex_file)
        trim_fq = os.path.join(
            tmp_dir, ADAPTER_TRIM_FQ_FORMAT.format(multiplex_name))
        log_file = os.path.join(logs_dir, "cutadapt.log")
        workflow.cut_adapters(None, config[params.ADAPTERS],
                              multiplex_file, trim_fq, log_file,
                              run_config)

        extract_trim_fq = os.path.join(
            tmp_dir, UMI_EXTRACT_FQ_FORMAT.format(multiplex_name))
        log_file = os.path.join(logs_dir, "umi_tools_extract.log")
        workflow.extract_barcodes_umis(None, trim_fq, extract_trim_fq,
                                       config[params.UMI_REGEXP],
                                       log_file, run_config)

        deplex_dir = os.path.join(
            tmp_dir,
            DEPLEX_DIR_FORMAT.format(multiplex_name))
        log_file = os.path.join(logs_dir, "demultiplex_fastq.log")
        workflow.demultiplex_fastq(extract_trim_fq,
                                   sample_sheet_file, deplex_dir,
                                   log_file, run_config)

        if not is_dry_run:
            num_reads_file = os.path.join(deplex_dir,
                                          demultiplex_fastq.NUM_READS_FILE)
            num_reads = sample_sheets.load_deplexed_sample_sheet(
                num_reads_file)
            samples = sample_sheets.get_non_zero_deplexed_samples(num_reads)
            if not samples:
                raise Exception(
                    "No non-empty sample files were produced by demultiplexing")
        else:
            # If doing a dry run, use the sample_sheet_file to get
            # the names of the samples.
            sample_sheet = sample_sheets.load_sample_sheet(sample_sheet_file)
            samples = list(sample_sheet[sample_sheets.SAMPLE_ID])
            if not samples:
                raise Exception(
                    "No samples are specified in {}".format(sample_sheet_file))
        # Use get_fastq_filename to deduce sample-specific files
        # output by demultiplex_fastq.py - demultiplex_fastq.py
        # uses this function too.
        sample_files = {sample: fastq.get_fastq_filename(sample)
                        for sample in samples}
        processed_samples = process_samples(
            sample_files, deplex_dir, index_dir, r_rna_index,
            orf_index, True, config, tmp_dir, out_dir, run_config,
            False)
        if not processed_samples:
            raise Exception("No samples were processed successfully")

    log_file = os.path.join(logs_dir, "collate_tpms.log")
    workflow.collate_tpms(
        out_dir, processed_samples, True, log_file, run_config)

    LOGGER.info("Completed")


def prep_riboviz(r_scripts, config_yaml, is_dry_run=False):
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
    LOGGER.info(provenance.get_provenance_str(__file__, " "))
    try:
        run_workflow(r_scripts, config_yaml, is_dry_run)
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
    sys.argv.pop(0)  # Remove program.
    is_dry_run_arg = (sys.argv[0] == DRY_RUN)
    if is_dry_run_arg:
        sys.argv.pop(0)
    r_scripts_arg = sys.argv[0]
    config_yaml_arg = sys.argv[1]
    exit_code = prep_riboviz(r_scripts_arg,
                             config_yaml_arg,
                             is_dry_run_arg)
    sys.exit(exit_code)
