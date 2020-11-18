#!/usr/bin/env python
"""
Run the workflow.

Usage::

    python -m riboviz.tools.prep_riboviz [-h] -c CONFIG_FILE [-d]

    -h, --help            show this help message and exit
    -c CONFIG_FILE, --config-file CONFIG_FILE
                          Configuration file
    -d, --dry-run         Dry run?

Example::

    python -m riboviz.tools.prep_riboviz
        -c vignette/vignette_config.yaml

If ``--dry-run`` is provided then the commands submitted to the
operating system, will not be executed. This can be useful for seeing
what commands will be run, validating the configuration, without
actually running the commands, and having a bash script that can be
run directly in future.

The following exit codes are returned:

* 0: Processing successfully completed.
* 1: A file does not seem to exist.
* 2: Errors occurred loading or accessing configuration e.g. missing
  configuration parameters, inconsistent configuration parameters.
* 3: Error occurred during processing.
* 4: User is using Python 2.

Operation is as follows.

:py:func:`invoke_prep_riboviz`:

* Parses command-line options using
  :py:func:`parse_command_line_options`.
* Invokes :py:func:`prep_riboviz` then exits using the exit code
  returned.

:py:func:`prep_riboviz`:

* Checks Python version.
* Writes provenance (file name, Git commit hash and date) to workflow
  log file.
* Invokes :py:func:`run_workflow` to run the workflow.
* Catches and logs errors thrown by :py:func:`run_workflow` and maps
  these to the corresponding exit codes.

:py:func:`run_workflow`:

* Reads configuration information from configuration file.
* Sets up command file into which commands sent to the operating
  system are to be written (a bash script).
* Creates index, temporary, output and logs directories.
* Checks for existence of all non-sample-specific data files to avoid
  having to repeatedly recheck these when processing each sample. This
  also allows the workflow to terminate before processing any samples,
  if any of these files cannot be found.
* Builds HISAT2 indices, if requested, using ``hisat2 build``.
* Checks if non-multiplexed FASTQ sample files or a multiplexed FASTQ
  sample file have been specified.
* If non-multiplexed sample files have been specified:
    - Processes the samples using :py:func:`process_samples`.
* If a multiplexed sample file has been specified:
    - Cuts out sequencing library adapters using ``cutadapt``.
    - Extracts barcodes and UMIs using ``umi_tools extract``, if
      requested, using a UMI-tools-compliant regular expression
      pattern.
    - Demultiplexes the sample file with reference to the sample
      sheet, using :py:mod:`riboviz.tools.demultiplex_fastq`.
    - Parses the number of reads file produced when demultiplexing
      samples to get the names of all samples which had 1 or more
      reads.
    - Determines the names of the corresponding demultiplexed sample
      files.
    - Processes the demultiplexed samples using
      :py:func:`process_samples`.
* Collates TPMs across all processed samples using
  ``collate_tpms.R``.
* Counts the reads at each step using
  :py:mod:`riboviz.tools.count_reads`.

:py:func:`process_samples`:

* Processes a collection of FASTQ sample files. For each sample:
    - Checks that the sample exists (if doing a dry run with
      multiplexed data this check is skipped as the file won't exist).
    - Creates sample-specific temporary, output and logs directories.
    - Processes the sample using :py:func:`process_sample`.
* The number of successfully and unsuccessfully processed samples are
  counted and the number of successfully processed samples is
  returned.

:py:func:`process_sample`:

* Processes a single FASTQ sample file.
* Cuts out sequencing library adapters using ``cutadapt``. If the
  samples arise from demultiplexed FASTQ file then this step in
  skipped as the adapters have already been cut.
* Extracts barcodes and UMIs using ``umi_tools extract``, if
  requested, using a UMI-tools-compliant regular expression
  pattern. If the samples arise from demultiplexed FASTQ file then
  this step is skipped as thew barcodes and UMIs have already been
  extracted.
* Removes rRNA or other contaminating reads by alignment to rRNA index
  files using ``hisat2``.
* Aligns remaining reads to ORFs index files using ``hisat2``.
* Trims 5' mismatches from reads and remove reads with more than 2
  mismatches using :py:mod:`riboviz.tools.trim_5p_mismatch`.
* Sorts resultant BAM file using ``samtools view | samtools sort``.
* Indexes resultant BAM file using ``samtools index``.
* If deduplication has been requested:
    - Outputs UMI groups pre-deduplication using ``umi_tools group``,
      if requested.
    - Deduplicates UMIs using ``umi_tools dedup``, and outputs
      deduplication statistics, if requested.
    - Indexes resultant BAM file using ``samtools index``.
    - Outputs UMI groups post-deduplication using ``umi_tools group``,
      if requested.
* Exports bedgraph files for plus and minus strands, if requested,
  using ``bedtools genomecov``.
* Makes length-sensitive alignments in compressed h5 format using
  ``bam_to_h5.R``.
* Generates summary statistics, and analyses and QC plots for both RPF
  and mRNA datasets using ``generate_stats_figs.R``. This includes
  estimated read counts, reads per base, and transcripts per million
  for each ORF in each sample.
* Writes intermediate files produced above into a sample-specific
  temporary directory.
* Writes output files produced above into a sample-specific output
  directory.

:py:const:`LOGGER` is used to log information in the workflow log
file.

Invocation of specific steps in the workflow is handled by calls to
functions in :py:mod:`riboviz.workflow`.

The script can parallelize parts of its operation over many processes
(``num_processes``):

* This value is used to configure ``hisat2``, ``samtools sort``,
  ``bam_to_h5.R`` and ``generate_stats_figs.R``.
* For ``cutadapt``, the number of available processors on the host will
  be used.
"""
import argparse
from datetime import datetime
import errno
import logging
import os
import os.path
import sys
import yaml
import riboviz
from riboviz import demultiplex_fastq
from riboviz import fastq
from riboviz import h5
from riboviz import logging_utils
from riboviz import params
from riboviz import provenance
from riboviz import sam_bam
from riboviz import sample_sheets
from riboviz import utils
from riboviz import workflow
from riboviz import workflow_files
from riboviz.utils import value_in_dict


DRY_RUN = "--dry-run"
""" Dry run command-line flag. """

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
EXIT_PYTHON_2_ERROR = 4
""" User is using Python 2. """

LOG_FORMAT = "{:02d}_{}"
""" Step-specific log file name format. """

logging_utils.configure_logging()
LOGGER = logging.getLogger(__name__)
""" Logger. """


def process_sample(sample, sample_fastq, index_dir, r_rna_index,
                   orf_index, is_trimmed, config, tmp_dir, out_dir,
                   logs_dir, run_config):
    """
    Process a single FASTQ sample file.

    (relevant configuration parameters are shown in brackets).

    * Processes a single FASTQ sample file.
    * Uses a step counter to number step-specific log files.
    * Cuts out sequencing library adapters (``adapters``) using \
      ``cutadapt`` (via :py:func:`riboviz.workflow.cut_adapters`).
        - If the samples arise from demultiplexed FASTQ file then this
          step is skipped as the adapters have already been cut.
    * Extracts barcodes and UMIs using ``umi_tools extract``, if \
      requested (``extract_umis``), using a UMI-tools-compliant \
      regular expression pattern (``umi_regexp``) (via \
      :py:func:`riboviz.workflow.extract_barcodes_umis`).
        - If the samples arise from demultiplexed FASTQ file then this
          step is skipped as thew barcodes and UMIs have already been
          extracted.
    * Removes rRNA or other contaminating reads by alignment to rRNA
      index files (``rrna_index_prefix``) using ``hisat2`` (via
      :py:func:`riboviz.workflow.map_to_r_rna`).
    * Aligns remaining reads to ORFs index files
      (``orf_index_prefix``). using ``hisat2`` (via
      :py:func:`riboviz.workflow.map_to_orf`).
    * Trims 5' mismatches from reads and remove reads with more than 2
      mismatches using :py:mod:`riboviz.tools.trim_5p_mismatch` (via
      :py:func:`riboviz.workflow.trim_5p_mismatches`).
    * Sorts resultant BAM file using ``samtools view | samtools sort``
      (via :py:func:`riboviz.workflow.sort_bam`).
    * Indexes resultant BAM file using ``samtools index`` (via
      :py:func:`riboviz.workflow.index_bam`).
    * If deduplication has been requested (``dedup_umis``):
        - Outputs UMI groups pre-deduplication using ``umi_tools
          group`` if requested (``group_umis``) (via
          :py:func:`riboviz.workflow.group_umis`).
        - Deduplicates UMIs using ``umi_tools dedup`` (via
          :py:func:`riboviz.workflow.deduplicate_umis`), and
          outputs deduplication statistics, if requested
          (``dedup_stats``).
        - Indexes resultant BAM file using ``samtools index`` (via
          :py:func:`riboviz.workflow.index_bam`).
        - Outputs UMI groups post-deduplication using ``umi_tools
          group``, if requested (``group_umis``) (via
          :py:func:`riboviz.workflow.group_umis`).
    * Exports bedgraph files for plus and minus strands, if requested
      (``make_bedgraph``), using ``bedtools genomecov``. (via
      :py:func:`riboviz.workflow.make_bedgraph`).
    * Makes length-sensitive alignments in compressed h5 format using
      ``bam_to_h5.R``. (via :py:func:`riboviz.workflow.bam_to_h5`).
    * Generates summary statistics, and analyses and QC plots for both
      RPF and mRNA datasets using ``generate_stats_figs.R``. This
      includes estimated read counts, reads per base, and transcripts
      per million for each ORF in each sample. (via
      :py:func:`riboviz.workflow.generate_stats_figs`).
    * Writes intermediate files produced above into a sample-specific
      temporary directory.
    * Writes output files produced above into a sample-specific output
      directory.

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
    :param is_trimmed: Have adapters been cut and barcodes \
    and UMIs extracted?
    :type is_trimmed: bool
    :param config: Workflow configuration
    :type config: dict
    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    :param out_dir: Output directory
    :type out_dir: str or unicode
    :param logs_dir: Logs directory
    :type logs_dir: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :raise FileNotFoundError: if ``sample_fastq`` or a third-party \
    tool cannot be found
    :raise AssertionError: if invocation of a third-party tool \
    returns non-zero exit code
    :raise KeyError: if ``config`` is missing required configuration
    """
    LOGGER.info("Processing sample: %s", sample)
    step = 1
    LOGGER.info("Processing file: %s", sample_fastq)
    is_extract_umis = value_in_dict(params.EXTRACT_UMIS, config)
    if is_trimmed:
        LOGGER.info("Skipping adaptor trimming and barcode/UMI extraction")
        trim_fq = sample_fastq
    else:
        log_file = os.path.join(logs_dir,
                                LOG_FORMAT.format(step, "cutadapt.log"))
        trim_fq = os.path.join(tmp_dir, workflow_files.ADAPTER_TRIM_FQ)
        workflow.cut_adapters(config[params.ADAPTERS],
                              sample_fastq, trim_fq, log_file,
                              run_config)
        step += 1

        if is_extract_umis:
            extract_trim_fq = os.path.join(tmp_dir,
                                           workflow_files.UMI_EXTRACT_FQ)
            log_file = os.path.join(
                logs_dir,
                LOG_FORMAT.format(step, "umi_tools_extract.log"))
            workflow.extract_barcodes_umis(trim_fq,
                                           extract_trim_fq,
                                           config[params.UMI_REGEXP],
                                           log_file,
                                           run_config)
            trim_fq = extract_trim_fq
            step += 1

    non_r_rna_trim_fq = os.path.join(tmp_dir, workflow_files.NON_RRNA_FQ)
    r_rna_map_sam = os.path.join(tmp_dir, workflow_files.RRNA_MAP_SAM)
    log_file = os.path.join(logs_dir,
                            LOG_FORMAT.format(step, "hisat2_rrna.log"))
    workflow.map_to_r_rna(trim_fq, index_dir, r_rna_index,
                          r_rna_map_sam, non_r_rna_trim_fq, log_file,
                          run_config)
    step += 1

    orf_map_sam = os.path.join(tmp_dir, workflow_files.ORF_MAP_SAM)
    unaligned_fq = os.path.join(tmp_dir, workflow_files.UNALIGNED_FQ)
    log_file = os.path.join(logs_dir,
                            LOG_FORMAT.format(step, "hisat2_orf.log"))
    workflow.map_to_orf(non_r_rna_trim_fq, index_dir,
                        orf_index, orf_map_sam, unaligned_fq,
                        log_file, run_config)
    step += 1

    orf_map_sam_clean = os.path.join(tmp_dir, workflow_files.ORF_MAP_CLEAN_SAM)
    trim_5p_mismatch_tsv = os.path.join(
        tmp_dir, workflow_files.TRIM_5P_MISMATCH_TSV)
    log_file = os.path.join(logs_dir,
                            LOG_FORMAT.format(step, "trim_5p_mismatch.log"))
    workflow.trim_5p_mismatches(orf_map_sam,
                                orf_map_sam_clean,
                                trim_5p_mismatch_tsv, log_file,
                                run_config)
    step += 1

    log_file = os.path.join(logs_dir,
                            LOG_FORMAT.format(step, "samtools_view_sort.log"))
    sample_out_prefix = os.path.join(out_dir, sample)

    is_dedup_umis = value_in_dict(params.DEDUP_UMIS, config)
    if is_dedup_umis:
        # Create BAM file in temporary directory
        sample_bam = os.path.join(tmp_dir, workflow_files.PRE_DEDUP_BAM)
    else:
        sample_bam = sam_bam.BAM_FORMAT.format(sample_out_prefix)
        sample_out_bam = sample_bam

    workflow.sort_bam(orf_map_sam_clean, sample_bam, log_file,
                      run_config)
    step += 1
    log_file = os.path.join(logs_dir,
                            LOG_FORMAT.format(step, "samtools_index.log"))
    workflow.index_bam(sample_bam, log_file, run_config)
    step += 1

    if is_dedup_umis:
        if not is_extract_umis:
            LOGGER.warning(
                "WARNING: dedup_umis was TRUE but extract_umis was FALSE.")
        is_group_umis = value_in_dict(params.GROUP_UMIS, config)
        if is_group_umis:
            umi_groups = os.path.join(tmp_dir,
                                      workflow_files.PRE_DEDUP_GROUPS_TSV)
            log_file = os.path.join(
                logs_dir,
                LOG_FORMAT.format(step, "umi_tools_group.log"))
            workflow.group_umis(sample_bam, umi_groups,
                                log_file, run_config)
            step += 1

        sample_out_bam = sam_bam.BAM_FORMAT.format(sample_out_prefix)
        log_file = os.path.join(
            logs_dir,
            LOG_FORMAT.format(step, "umi_tools_dedup.log"))
        is_dedup_stats = True
        if params.DEDUP_STATS in config:
            is_dedup_stats = value_in_dict(params.DEDUP_STATS, config)
        dedup_stats_prefix = None
        if is_dedup_stats:
            dedup_stats_prefix = os.path.join(
                tmp_dir,
                workflow_files.DEDUP_STATS_PREFIX)
        workflow.deduplicate_umis(
            sample_bam, sample_out_bam, dedup_stats_prefix,
            log_file, run_config)
        step += 1

        log_file = os.path.join(
            logs_dir,
            LOG_FORMAT.format(step, "samtools_index.log"))
        workflow.index_bam(sample_out_bam, log_file, run_config)
        step += 1

        if is_group_umis:
            umi_groups = os.path.join(tmp_dir,
                                      workflow_files.POST_DEDUP_GROUPS_TSV)
            log_file = os.path.join(
                logs_dir,
                LOG_FORMAT.format(step, "umi_tools_group.log"))
            workflow.group_umis(sample_out_bam, umi_groups,
                                log_file, run_config)
            step += 1

    is_make_bedgraph = value_in_dict(params.MAKE_BEDGRAPH, config)
    if is_make_bedgraph:
        log_file = os.path.join(
            logs_dir,
            LOG_FORMAT.format(step, "bedtools_genome_cov_plus.log"))
        plus_bedgraph = os.path.join(out_dir, workflow_files.PLUS_BEDGRAPH)
        workflow.make_bedgraph(sample_out_bam, plus_bedgraph,
                               True, log_file, run_config)
        step += 1

        log_file = os.path.join(
            logs_dir,
            LOG_FORMAT.format(step, "bedtools_genome_cov_minus.log"))
        minus_bedgraph = os.path.join(out_dir, workflow_files.MINUS_BEDGRAPH)
        workflow.make_bedgraph(sample_out_bam, minus_bedgraph,
                               False, log_file, run_config)
        step += 1

    orf_gff_file = config[params.ORF_GFF_FILE]
    log_file = os.path.join(logs_dir,
                            LOG_FORMAT.format(step, "bam_to_h5.log"))
    sample_out_h5 = h5.H5_FORMAT.format(sample_out_prefix)
    workflow.bam_to_h5(sample_out_bam, sample_out_h5,
                       orf_gff_file, config, log_file, run_config)
    step += 1

    log_file = os.path.join(
        logs_dir,
        LOG_FORMAT.format(step, "generate_stats_figs.log"))
    workflow.generate_stats_figs(sample_out_h5, out_dir,
                                 config, log_file, run_config)

    LOGGER.info("Finished processing sample: %s", sample_fastq)


def process_samples(samples, in_dir, index_dir, r_rna_index,
                    orf_index, is_trimmed, config, tmp_dir, out_dir,
                    logs_dir, run_config, check_samples_exist=True):
    """
    Process FASTQ sample files. Any exceptions in the processing of
    any sample are logged but are not thrown from this function.

    (relevant configuration parameters are shown in brackets).

    * Processes a collection of FASTQ sample files. For each sample:
        - Checks that the sample exists (if doing a dry run with
          multiplexed data this check is skipped as the file won't
          exist).
        - Creates sample-specific temporary, output and logs
          directories (via
          :py:func:`riboviz.workflow.create_directory`).
        - Processes the sample using :py:func:`process_sample`.
        - If any errors arise when processing the sample, the error is
          logged but processing continues onto the other samples.
    * The number of successfully and unsuccessfully processed samples
      are counted and the number of successfully processed samples is
      returned.

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
    :param is_trimmed: Have adapters been cut and barcodes \
    and UMIs extracted?
    :type is_trimmed: bool
    :param config: Workflow configuration
    :type config: dict
    :param tmp_dir: Temporary directory
    :type tmp_dir: str or unicode
    :param out_dir: Output directory
    :type out_dir: str or unicode
    :param run_config: Run-related configuration
    :type run_config: RunConfigTuple
    :param check_samples_exist: If ``run_config.is_dry_run`` \
    is ``True``, should a check be made for the existence of sample \
    files?
    :type check_samples_exist: bool
    :return: Names of successfully-processed samples
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
            sample_logs_dir = os.path.join(logs_dir, sample)
            for directory in [sample_tmp_dir,
                              sample_out_dir,
                              sample_logs_dir]:
                workflow.create_directory(directory,
                                          run_config.cmd_file,
                                          run_config.is_dry_run)
            process_sample(sample, sample_fastq, index_dir,
                           r_rna_index, orf_index, is_trimmed,
                           config, sample_tmp_dir, sample_out_dir,
                           sample_logs_dir, run_config)
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


def run_workflow(config_file, is_dry_run=False):
    """
    Run the workflow.

    (relevant configuration parameters are shown in brackets).

    * Reads configuration information from configuration file.
    * Sets up command file (``cmd_file``) into which commands sent to
      the operating system are to be written (a bash script).
    * Creates index (``dir_in``), temporary (``dir_tmp``), output
      (``dir_output``) and logs (``dir_logs``) directories (via
      :py:func:`riboviz.workflow.create_directory`).
    * Checks existence of all non-sample-specific data files
      (``rrna_fasta_file``, ``orf_fasta_file``, ``orf_gff_file``,
      ``features_file``, ``t_rna_file``, ``codon_positions_file``,
      ``asite_disp_length_file``) to avoid having to repeatedly
      recheck these when processing each sample. This also allows the
      workflow to terminate before processing any samples, if any of
      these files cannot be found.
    * Sets up a run configuration \
      (:py:class:`riboviz.workflow.RunConfigTuple`) which provides \
      each function with common configuration, notably:
        - The command file (``cmd_file``) in which the commands sent
          to the operating system are to be written.
        - Whether the invocation is part of a dry run? If so the
          commands are recorded in the command file but are not
          submitted to the operating system.
        - The number of processes available (``num_processes``), for
          commands which allow the number of processeses to use to be
          specified.
        - R scripts directory.
    * Builds HISAT2 indices, if requested (``build_indices``), using
      ``hisat2 build``` and writes these into the index directory
      (``dir_index``) (via
      :py:func:`riboviz.workflow.build_indices`).
    * Checks if non-multiplexed FASTQ sample files (``fq_files``) or a
      multiplexed FASTQ sample file (``multiplex_fq_files``) have been
      specified.
    * If non-multiplexed sample files have been specified:
        - Processes samples using :py:func:`process_samples`.
        - Raises an error if no sample was processed successfully.
    * If a multiplexed sample file has been specified:
        - Checks for the existence of the multiplexed sample file and
          a sample sheet file (``sample_sheet``).
        - Cuts out sequencing library adapters (``adapters``) using
          ``cutadapt`` (via :py:func:`riboviz.workflow.cut_adapters`).
        - Extracts barcodes and UMIs using ``umi_tools extract``, if
          requested (``extract_umis``), using a UMI-tools-compliant
          regular expression pattern (``umi_regexp``) (via
          :py:func:`riboviz.workflow.extract_barcodes_umis`).
        - Demultiplexes the sample file with reference to the sample
          sheet, using :py:mod:`riboviz.tools.demultiplex_fastq` (via
          :py:func:`riboviz.workflow.demultiplex_fastq`) .
        - Parses the number of reads file produced when demultiplexing
          samples to get the names of all samples which had 1 or more
          reads.
        - Determines the names of the corresponding demultiplexed
          sample files.
        - Processes the demultiplexed samples using
          :py:func:`process_samples`.
        - Raises an error if no sample was processed successfully.
    * Collates TPMs across all processed samples using
      ``collate_tpms.R``  and writes these into the output directory
      (via :py:mod:`riboviz.workflow.collate_tpms`).
    * Counts the reads at each step using
      :py:mod:`riboviz.tools.count_reads`) and writes these into the
      output directory (via :py:mod:`riboviz.workflow.count_reads`).

    :param config_file: Configuration file path
    :type config_file: str or unicode
    :param is_dry_run: Is this a dry run? (if ``True`` workflow \
    commands will not be submitted to the operating system for \
    execution)
    :type is_dry_run: bool
    :raise FileNotFoundError: if an input file cannot be found
    :raise KeyError: if a configuration parameter is missing
    :raise ValueError: if a configuration parameter has an \
    invalid value
    :raise TypeError: if a configuration parameter has an invalid type
    :raise Exception: if any other error arises
    """
    LOGGER.info("Configuration file: %s", config_file)
    with open(config_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)

    if value_in_dict(params.CMD_FILE, config):
        cmd_file = config[params.CMD_FILE]
    else:
        cmd_file = workflow_files.DEFAULT_CMD_FILE
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

    # Check existence of all non-sample-specific files to avoid having
    # to recheck them when processing each sample.
    for input_key in [params.RRNA_FASTA_FILE,
                      params.ORF_FASTA_FILE,
                      params.ORF_GFF_FILE]:
        input_file = config[input_key]
        if not os.path.exists(input_file):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), input_file)
    # Optional files.
    for input_key in [params.FEATURES_FILE,
                      params.T_RNA_FILE,
                      params.CODON_POSITIONS_FILE,
                      params.ASITE_DISP_LENGTH_FILE]:
        if value_in_dict(input_key, config):
            input_file = config[input_key]
            if not os.path.exists(input_file):
                raise FileNotFoundError(errno.ENOENT,
                                        os.strerror(errno.ENOENT),
                                        input_file)

    run_config = workflow.RunConfigTuple(
        riboviz.R_SCRIPTS,
        cmd_file,
        is_dry_run,
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
            False, config, tmp_dir, out_dir, logs_dir, run_config)
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
            tmp_dir, workflow_files.ADAPTER_TRIM_FQ_FORMAT.format(multiplex_name))
        log_file = os.path.join(logs_dir, "cutadapt.log")
        workflow.cut_adapters(config[params.ADAPTERS],
                              multiplex_file, trim_fq, log_file,
                              run_config)

        extract_trim_fq = os.path.join(
            tmp_dir,
            workflow_files.UMI_EXTRACT_FQ_FORMAT.format(multiplex_name))
        log_file = os.path.join(logs_dir, "umi_tools_extract.log")
        workflow.extract_barcodes_umis(trim_fq, extract_trim_fq,
                                       config[params.UMI_REGEXP],
                                       log_file, run_config)

        deplex_dir = os.path.join(
            tmp_dir,
            workflow_files.DEPLEX_DIR_FORMAT.format(multiplex_name))
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
        # output by riboviz.tools.demultiplex_fastq which
        # uses this function too.
        file_format = fastq.FASTQ_FORMATS[utils.get_file_ext(extract_trim_fq)]
        sample_files = {sample: file_format.format(sample)
                        for sample in samples}
        processed_samples = process_samples(
            sample_files, deplex_dir, index_dir, r_rna_index,
            orf_index, True, config, tmp_dir, out_dir, logs_dir,
            run_config,
            False)
        if not processed_samples:
            raise Exception("No samples were processed successfully")

    log_file = os.path.join(logs_dir, "collate_tpms.log")
    workflow.collate_tpms(out_dir, processed_samples, log_file, run_config)

    if value_in_dict(params.COUNT_READS, config):
        log_file = os.path.join(logs_dir, "count_reads.log")
        read_counts_file = os.path.join(out_dir,
                                        workflow_files.READ_COUNTS_FILE)
        workflow.count_reads(config_file, in_dir, tmp_dir, out_dir,
                             read_counts_file, log_file, run_config)
    LOGGER.info("Completed")


def prep_riboviz(config_file, is_dry_run=False):
    """
    Run the workflow.

    * Checks whether Python 3 is being used and, if not, returns an
      exit code (:py:const:`EXIT_PYTHON_2_ERROR`).
    * Writes provenance (file name, Git commit hash and date) to
      workflow log file.
    * Invokes :py:func:`run_workflow` to run the workflow.
    * Catches and logs errors thrown by :py:func:`run_workflow` and
      maps these to the corresponding exit codes.
    * Returns :py:const:`EXIT_OK`, if no errors arose, or an exit code
      reflecting the error that arose.

    Exit codes are as follows:

    * :py:const:`EXIT_OK`
    * :py:const:`EXIT_FILE_NOT_FOUND_ERROR`
    * :py:const:`EXIT_CONFIG_ERROR`
    * :py:const:`EXIT_PROCESSING_ERROR`
    * :py:const:`EXIT_PYTHON_2_ERROR`

    :param config_file: Configuration file path
    :type config_file: str or unicode
    :param is_dry_run: Is this a dry run? (if ``True`` workflow \
    commands will not be submitted to the operating system for \
    execution)
    :return: exit code
    :rtype: int
    """
    LOGGER.info("Running under Python: %s", sys.version)
    if sys.version_info.major < 3:
        LOGGER.error("This script needs to be run under Python 3")
        return EXIT_PYTHON_2_ERROR
    LOGGER.info(provenance.write_provenance_to_str(__file__, " "))
    try:
        run_workflow(config_file, is_dry_run)
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
        LOGGER.error("Processing error: %s", config_file)
        exc_type, _, _ = sys.exc_info()
        LOGGER.exception(exc_type.__name__)
        return EXIT_PROCESSING_ERROR
    return EXIT_OK


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Run workflow")
    parser.add_argument("-c",
                        "--config-file",
                        dest="config_file",
                        required=True,
                        help="Configuration file")
    parser.add_argument("-d",
                        "--dry-run",
                        dest='is_dry_run',
                        action='store_true',
                        help="Dry run?")
    options = parser.parse_args()
    return options


def invoke_prep_riboviz():
    """
    Parses command-line options using
    :py:func:`parse_command_line_options` then invokes
    :py:func:`prep_riboviz` then exits using the exit code returned.

    Parse command-line options then invoke
    :py:func:`prep_riboviz`.
    """
    options = parse_command_line_options()
    config_file = options.config_file
    is_dry_run = options.is_dry_run
    exit_code = prep_riboviz(config_file, is_dry_run)
    sys.exit(exit_code)


if __name__ == "__main__":
    invoke_prep_riboviz()
