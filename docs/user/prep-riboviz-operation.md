# What the riboviz workflow does

This page describes what the riboviz workflow does.

Configuration parameters are shown in brackets and are described in [Configuring the riboviz workflow](./prep-riboviz-config.md).

---

## Coordinating local scripts and third-party components

The workflow prepares ribosome profiling data for by implementing a workflow that uses a combination of local Python and R scripts and third-party components. These are as follows:

* `hisat2-build`: build rRNA and ORF indices.
* `cutadapt`: cut adapters.
* `hisat2`: align reads.
* `riboviz.tools.trim_5p_mismatch`: trim 5' mismatches from reads and remove reads with more than a set number of mismatches (local script, in `riboviz/tools/`).
* `umi_tools` (`extract`, `dedup`, `group`): extract barcodes and UMIs, deduplicate reads and group reads.
* `riboviz.tools.demultiplex_fastq`: demultiplex multiplexed files (local script, in `riboviz/tools/`).
* `samtools` (`view`, `sort`, `index`): convert SAM files to BAM files and index.
* `bedtools` (`genomecov`): export transcriptome coverage as bedgraphs.
* `bam_to_h5.R`: convert BAM to compressed H5 format (local script, in `rscripts/`)
* `generate_stats_figs.R`: generate summary statistics, analyses plots and QC plots (local script, in `rscripts/`)
* `collate_tpms.R`: collate TPMs across samples (local script, in `rscripts/`)
* `riboviz.tools.count_reads`: count the number of reads (sequences) processed by specific stages of the workflow (local script, in `riboviz/tools/`).  
* `AnalysisOutputs.Rmd`: the `staticHTML` Nextflow process runs an R markdown document which generates an HTML output report for each sample processed, using the analysed data (local script, in `rmarkdown/`)

### Interactive Visualization

An additional script within the riboviz codebase allows the generation of interactive plots from the dataset after the riboviz workflow has been completed.  

A process within the main riboviz workflow (`createInteractiveVizParamsConfigFile`, within `prep-riboviz.nf`) generates a configuration yaml file specifically for this interactive visualization (`interactive_viz_config.yaml`, found in the `dir_out` directory), and users can then run `rscripts/run_shiny_server.R` using this additional YAML (as detailed within [How To Run the riboviz Interactive Data Visualization On Your Data](./run_shiny_server-operation.md)).  This generates a shiny app instance which users can view to explore their data.


---

## Process ribosome profiling sample data

If sample files (`fq_files`) are specified, then the workflow processes the sample files as follows:

1. Read configuration information from YAML configuration file.
2. Build hisat2 indices if requested (if `build_indices: TRUE`) using `hisat2 build` and save these into the index directory (`dir_index`).
3. Process each sample ID-sample file pair (`fq_files`) in turn:
   1. Cut out sequencing library adapters (`adapters`) using `cutadapt`.
   2. Extract UMIs using `umi_tools extract`, if requested (if `extract_umis: TRUE`), using a UMI-tools-compliant regular expression pattern (`umi_regexp`). The extracted UMIs are inserted into the read headers of the FASTQ records.
   3. Remove rRNA or other contaminating reads by alignment to rRNA index files (`rrna_index_prefix`) using `hisat2`.
   4. Align remaining reads to ORFs index files (`orf_index_prefix`). using `hisat2`.
   5. Trim 5' mismatches from reads and remove reads with more than 2 mismatches using `trim_5p_mismatch`, if requested (if `trim_5p_mismatches: TRUE`)
   6. Output UMI groups pre-deduplication using `umi_tools group` if requested (if `dedup_umis: TRUE` and `group_umis: TRUE`)
   7. Deduplicate reads using `umi_tools dedup`, if requested (if `dedup_umis: TRUE`), and output deduplication statistics, if requested (if `dedup_stats: TRUE`).  
   8. Output UMI groups post-deduplication using `umi_tools group` if requested (if `dedup_umis: TRUE` and `group_umis: TRUE`)
   9. Export bedgraph files for plus and minus strands, if requested (if `make_bedgraph: TRUE`) using `bedtools genomecov`.
   10. Write intermediate files produced above into a sample-specific directory, named using the sample ID, within the temporary directory (`dir_tmp`).
   11. Make length-sensitive alignments in compressed h5 format using `bam_to_h5.R`.
   12. Generate summary statistics, and analyses and QC plots for both RPF and mRNA datasets using `generate_stats_figs.R`. This includes estimated read counts, reads per base, and transcripts per million for each ORF in each sample.
   13. Write output files produced above into a sample-specific directory, named using the sample ID, within the output directory (`dir_out`).
4. Collate TPMs across results, using `collate_tpms.R` and write into output directory (`dir_out`). Only the results from successfully-processed samples are collated.
5. Count the number of reads (sequences) processed by specific stages if requested (if `count_reads: TRUE`).

[Workflow](../images/workflow.svg) (SVG) shows an images of the workflow with the key steps, inputs and outputs.

[Workflow with deduplication](../images/workflow-dedup.svg) (SVG) shows the workflow, if `dedup_umis: TRUE`.

---

## Process multiplexed ribosome profiling sample data

If a multiplexed file (`multiplex_fq_files`) is specified, then the workflow processes the multiplexed file as follows:

1. Read configuration information (as for 1. above).
2. Build hisat2 indices if requested (as for 2. above).
3. Read the multiplexed FASTQ file (`multiplex_fq_files`).
4. Cut out sequencing library adapters (`adapters`) using `cutadapt`.
5. Extract barcodes and UMIs using `umi_tools extract`, if requested (if `extract_umis: TRUE`), using a UMI-tools-compliant regular expression pattern (`umi_regexp`).
6. Demultiplex file with reference to the sample sheet (`sample_sheet`), using `demultiplex_fastq`. Sample IDs in the `SampleID` column in the sample sheet are used to name the demultiplexed files.
7. Process each demultiplexed FASTQ file which has one or more reads, in turn (as for 3.3 to 3.13 above)
8. Collate TPMs across results, using `collate_tpms.R` and write into output directory (`dir_out`) (as for 4. above.
9. Count the number of reads (sequences) processed by specific stages if requested (if `count_reads: TRUE`).

[Workflow with demultiplexing](../images/workflow-deplex.svg) (SVG) shows an images of the workflow with the key steps, inputs and outputs.

---

## Index files

Index files (HT2) are produced in the index directory (`dir_index`) (if `build_indices: TRUE`).

This directory provides easy access to the index files created and used in the workflow, and which can also be used for troubleshooting. By default, files in this index directory are **symbolic links** to files in the [Nextflow work/ directory](#nextflow-work-directory). This means that if the Nextflow `work/` directory is deleted, then the index files will no longer be accessible.

To request that Nextflow copy index files into the `dir_index` directory, set the `publish_index_tmp` parameter to `TRUE` in the workflow configuration file or provide the parameter `--publish_index_tmp` when running the workflow using Nextflow. Please note, however, that, as this creates copies of both the index files and temporary files (which can take up many gigabytes of space), setting `publish_index_tmp` is **not recommended** in general.

---

## Temporary files

Intermediate files are produced within the temporary directory (`dir_tmp`).

For each sample (`<SAMPLE_ID>`), intermediate files are produced in a sample-specific subdirectory (`<SAMPLE_ID>`):

* `trim.fq`: adapter trimmed reads (if `fq_files` and not `multiplex_fq_files` are specified).
* `nonrRNA.fq`: non-rRNA reads.
* `rRNA_map.sam`: rRNA-mapped reads.
* `unaligned.fq`: unaligned reads. These files can be used to find common contaminants or translated sequences not in the ORF annotation.
* `orf_map.sam`: ORF-mapped reads.
* `orf_map_clean.sam`: ORF-mapped reads with 5' mismatched nts trimmed (if `trim_5p_mismatches: TRUE`).
* `trim_5p_mismatch.tsv`: number of reads processed, discarded, trimmed and written when trimming 5' mismatches from reads and removing reads with more than a set number of mismatches (if `trim_5p_mismatches: TRUE`).
* `orf_map_clean.bam`: BAM file equivalent of `orf_map_clean.sam` if trimming is enabled (if `trim_5p_mismatches: TRUE`) OR `orf_map.sam` (if `trim_5p_mismatches: FALSE`). If deduplication is not enabled (if `dedup_umis: FALSE`) then this is copied to become the output file `<SAMPLE_ID>.bam` (see below).
* `orf_map_clean.bam.bai`: BAM index file for the above. If deduplication is not enabled (if `dedup_umis: FALSE`) then this is copied to become the output file `<SAMPLE_ID>.bam.bai` (see below).

If deduplication is enabled (if `dedup_umis: TRUE`) the following sample-specific files are also produced:

* `extract_trim.fq`: adapter trimmed reads with UMIs extracted (if `extract_umis: TRUE` and `fq_files` and not `multiplex_fq_files` are specified).
* `dedup.bam`: BAM file post deduplication. This is copied to become the output file `<SAMPLE_ID>.bam`.
* `dedup.bam.bai`: BAM index file for the above. This is copied to become the output file `<SAMPLE_ID>.bam.bai`.
* UMI groups pre- and post-deduplication (if `group_umis: TRUE`):
  - `pre_dedup_groups.tsv`: UMI groups before deduplication.
  - `post_dedup_groups.tsv`: UMI groups after deduplication.
* UMI deduplication statistics (if `dedup_stats: TRUE`):
  - `dedup_stats_edit_distance.tsv`: edit distance between UMIs at each position.
  - `dedup_stats_per_umi_per_position.tsv`: histogram of counts per position per UMI pre- and post-deduplication.
  - `dedup_stats_per_umi.tsv`: number of times each UMI was observed, total counts and median counts, pre- and post-deduplication
  - For more information on the `stats` files, see see UMI-tools [Dedup-specific options](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html) and [documentation on stats file #250](https://github.com/CGATOxford/UMI-tools/issues/250)

If a multiplexed file (`multiplex_fq_files`) is specified, then the following files and directories are also written into the temporary directory:

* `<FASTQ_FILE_NAME_PREFIX>_trim.fq`: FASTQ file post-adapter trimming, where `<FASTQ_FILE_NAME_PREFIX>` is the name of the file (without path or extension) in `multiplex_fq_files`.
* `<FASTQ_FILE_NAME_PREFIX>_extract_trim.fq`: `<FASTQ_FILE_NAME_PREFIX_trim.fq` post-barcode and UMI extraction (if `extract_umis: TRUE`).
* `<FASTQ_FILE_NAME_PREFIX>_deplex/`: demultiplexing results directory including:
   - `num_reads.tsv`: a tab-separated values file with columns:
     - `SampleID`, copied from the sample sheet.
     - `TagRead` (barcode), coped from the sample sheet.
     - `NumReads`, number of reads detected for each sample.
     - Row with `SampleID` with value `Unassigned` and `NumReads` value with the number of unassigned reads.
     - Row with `SampleID` with value `Total` and `NumReads` value with the total number of reads processed.
  - `<SAMPLE_ID>.fastq`: Files with demultiplexed reads, where `<SAMPLE_ID>` is a value in the `SampleID` column of the sample sheet. There will be one file per sample.
  - `Unassigned.fastq`: A FASTQ file with the reads that did not match any `TagRead` (barcode) in the sample sheet.

This directory provides easy access to the temporary, or intermediate, files created and used in the workflow which can also be used for troubleshooting. By default, files in this temporary directory are **symbolic links** to files in the [Nextflow work/ directory](#nextflow-work-directory). This means that if the Nextflow `work/` directory is deleted, then the temporary files will no longer be accessible.

To request that Nextflow copy temporary files into the `dir_tmp` directory, set the `publish_index_tmp` parameter to `TRUE` in the workflow configuration file or provide the parameter `--publish_index_tmp` when running the workflow using Nextflow. Please note, however, that, as this creates copies of both the index files and temporary files (which can take up many gigabytes of space), setting `publish_index_tmp` is **not recommended** in general.

---

## Output files

Output files are produced within the output directory (`dir_out`).

For each sample (`<SAMPLE_ID>`), intermediate files are produced in a sample-specific subdirectory (`<SAMPLE_ID>`):

* `<SAMPLE_ID>.bam`: BAM file of reads mapped to transcripts, which can be directly used in genome browsers.
* `<SAMPLE_ID>.bam.bai`: BAM index file for `<SAMPLE_ID>.bam`.
* `minus.bedgraph`: bedgraph of reads from minus strand (if `make_bedgraph: TRUE`).
* `plus.bedgraph`: bedgraph of reads from plus strand (if `make_bedgraph: TRUE`).
* `<SAMPLE_ID>.h5`, `<SAMPLE_ID>.h5.*`: length-sensitive alignments in compressed h5 format. The number of output files depends on the number of processes that `bam_to_h5.R` was run with (`num_processes`).
* `metagene_start_stop_read_counts.tsv`
* `metagene_start_stop_read_counts.pdf` (if `output_pdfs: TRUE`)
* `metagene_position_length_counts_5start.tsv`
* `read_counts_by_length.tsv`
* `read_counts_by_length.pdf` (if `output_pdfs: TRUE`)
* `nt_freq_per_read_position.tsv` (if `output_metagene_normalized_profile: TRUE`)
* `metagene_normalized_profile_start_stop.pdf` (if `output_pdfs: TRUE`)
* `metagene_normalized_profile_start_stop.tsv`
* `ORF_TPMs_vs_features.tsv` (if `features_file` was defined)
* `ORF_TPMs_vs_features.pdf` (if `features_file` was defined and `output_pdfs: TRUE`)
* `ORF_TPMs_and_counts.tsv`
* `normalized_density_APEsites_per_codon.tsv` (if `t_rna_file` and `codon_positions_file` were defined)
* `normalized_density_APEsites_per_codon.pdf` (if `t_rna_file` and `codon_positions-file` were defined and `output_pdfs: TRUE`)
 * `normalized_density_APEsites_per_codon_long.tsv`  (if `t_rna_file` and `codon_positions_file` were defined)
* `metagene_start_barplot_by_length.pdf` (if `output_pdfs: TRUE`)
* `metagene_start_ribogrid_by_length.pdf` (if `output_pdfs: TRUE`)
* `read_frame_per_ORF.tsv` (if `asite_disp_length_file` was defined)
* `read_frame_per_ORF_filtered.tsv` (if `asite_disp_length_file` was defined)
* `frame_proportions_per_ORF.pdf` (if `asite_disp_length_file` was defined and `output_pdfs: TRUE`)
 * `<SAMPLE_ID>_output_report.html` (if `run_static_html: TRUE`)

In addition, the following files are also put into the output directory:

* `TPMs_all_CDS_all_samples.tsv`: file with the transcripts per million (tpm) for all successfully processed samples.
* `read_counts_per_file.tsv`: a read counts file. (only if `count_reads: TRUE`).
* `interactive_viz_config.yaml`: this is a yaml file created by the workflow, for use with `rscripts/run_shiny_server.R` - an optional step which does not automatically run within the riboviz workflow and which allows users to generate interactive data visualization on the dataset.

More details on the output files can be found at [riboviz output files and figures](./riboviz-outputs.md).

---

## Nextflow `work/` directory

When Nextflow runs, it creates a unique step-specific directory for every step in the workflow. Each step-specific directory has symbolic links to the input files for the step and a bash script with the commands to be run by Nextflow for that step. Nextflow runs this bash script within this directory which creates the step's output files. Nextflow also creates files with output and error messages and the exit code. These step-specific directories are created within a Nextflow `work/` directory located, by default, within the directory within which Nextflow is run.

Every step-specific directory within Nextflow's `work/` directory has a name in common with the process identifier (e.g. `ad/1e7c54`) which is displayed when Nextflow runs. These identifiers can also be accessed using the `nextflow log` command. For example, to display the process identifier and `work/` subdirectory of a step called `collateTpms (WT3AT, WTnone)` in a workflow run named `big_majorana`:

```console
$ nextflow log big_majorana -f hash,workdir -filter "name == 'collateTpms (WT3AT, WTnone)'"

38/784d89	/home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73
```

This shows both the unique identifier, termed a "hash", of this step that was shown when the workflow was run, and also the corresponding subdirectory within `work/` for that step. Note that the identifier is a prefix of the subdirectory under `work/`.

These step-specific directories are where Nextflow runs each step. Each subdirectory has:

* Input files. These are symbolic links to the input files for the step which, depending on the step, can be:
  - Output files in other `work/` subdirectories. For example, the directory for an `hisat2rRNA` step will have input files which are symbolic links to the output files produced by a `cutAdapters` step.
  - Input files for the workflow. For example, the directory for a `cutAdapters` step will have an input file which is a symbolic link to a sample file in `vignettte/input`.
* Bash script (`.command.sh`) containing the specific commands invoked by Nextflow for that step.
* Output files, from the invocation of the step.
* Standard output (`.command.out`) and standard error (`.command.err`) files and combined standard output and standard error (`.command.log`) containing the output and error messages printed during invocation of the bash script and captured by Nextflow.
* Exit code (`.exitcode`) output from running the bash script.

For example, for a p `ad/1e7c54`, an invocation of `hisat2rRNA` for sample `WTnone`, the `work/` directory would include:

```console
$ ls -1a work/ad/1e7c54a889f21451cb07d29655e0be/ -printf '%P\t%l\n' | so
.command.begin
.command.err
.command.log
.command.out
.command.run
.command.sh
.command.trace
.exitcode
nonrRNA.fq
rRNA_map.sam
trim.fq	/home/ubuntu/riboviz/work/5b/ff17ba3e19c5d659e54b04b08dab85/trim.fq
yeast_rRNA.1.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.1.ht2
yeast_rRNA.2.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.2.ht2
yeast_rRNA.3.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.3.ht2
yeast_rRNA.4.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.4.ht2
yeast_rRNA.5.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.5.ht2
yeast_rRNA.6.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.6.ht2
yeast_rRNA.7.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.7.ht2
yeast_rRNA.8.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.8.ht2
```

The `.ht2` files are symbolic links to the outputs of the step `e5/ccf3e6`, an invocation of the process `buildIndicesrRNA`.

The riboviz workflow uses Nextflow's [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) directive which allows files to be published to specific directories outwith `work/`.

For index and temporary files, `publishDir` is configured using the value of the `publish_index_tmp` parameter. If `FALSE` then files in the index (`dir_index`) and temporary (`dir_tmp`) directories are symbolically linked to those in `work/`. If `TRUE` then they are copied. Output files are always copied from `work/` into the output (`dir_out`) directory specified in the workflow configuration file.

**Caution:**

* If `publish_index_tmp` is `FALSE` and the `work/` directory is deleted then the index and temporary files will no longer be accessible.
* If `publish_index_tmp` is `TRUE` then `dir_tmp` will be populated with copies of temporary files which can take up many gigabytes of space.  Setting `publish_index_tmp` is **not recommended** in general.
* If the `work/` folder is deleted, then certain information will no longer be accessible via `nextflow log`.
* If the `work/` folder is deleted, then `the `-resume` flag has no effect and the whole workflow will be rerun.

### `Missing` files

If an optional file for `generate_stats_figs.R` is not provided within the YAML configuration file then a `Missing_<PARAM>` file (for example `Missing_features_file`) is created within the `work/` directory for a `generateStatsFigs` steps. This symbolically links to a non-existent `Missing_<PARAM>` file in the current directory. This is not an issue since the files will not be passed onto `generate_stats_figs.R` and no attempt is made to use them. They are a side-effect of using the Nextflow pattern for optional inputs, [optional inputs](https://github.com/nextflow-io/patterns/blob/master/optional-input.nf).

---

## Nextflow log files

Information on the execution of the workflow is added to a file `.nextflow.log`. Log files from previous runs are in files named `.nextflow.log.1`, `.nextflow.log.2` etc. Every time Nextflow is run, the log file names are adjusted - on each successive run `.nextflow.log.1` becomes `.nextflow.log.2` and `.nextflow.log` becomes `.nextflow.log.1`).

Log files for the invocation of each step in the workflow are captured in the [Nextflow `work/` directory](#nextflow-work-directory) described above.
