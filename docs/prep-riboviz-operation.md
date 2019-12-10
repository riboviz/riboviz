# What the RiboViz workflow does

This page describes what `prep_riboviz.py` does. 

Configuration parameters are shown in brackets and are described in [Configuring the RiboViz workflow](./prep-riboviz-config.md).

---

## Coordinating local scripts and third-party components

`prep_riboviz.py` prepares ribosome profiling data for by implementing a workflow that uses a combination of local Python and R scripts and third-party components. These are as follows:

* `hisat2-build`: build rRNA and ORF indices.
* `cutadapt`: cut adapters.
* `hisat2`: align reads.
* `trim_5p_mismatch.py`: trim reads (local script, in `riboviz/tools/`).
* `umi_tools` (`extract`, `dedup`, `group`): extract barcodes and UMIs, deduplicate reads and group reads.
* `demultiplex_fastq.py`: demultiplex multiplexed files (local script, in `riboviz/tools/`).
* `samtools` (`view`, `sort`, `index`): convert SAM files to BAM files and index.
* `bedtools` (`genomecov`): export transcriptome coverage as bedgraphs.
* `bam_to_h5.R`: convert BAM to compressed H5 format (local script, in `rscripts/`)
* `generate_stats_figs.R`: generate summary statistics, analyses plots and QC plots (local script, in `rscripts/`)
* `collate_tpms.R`: collate TPMs across samples (local script, in `rscripts/`)

---

## Process ribosome profiling sample data

If sample files (`fq_files`) are specified, then `prep_riboviz.py` processes the sample files as follows:

1. Read configuration information from YAML configuration file.
2. Build hisat2 indices if requested (if `build_indices: TRUE`) using `hisat2 build` and save these into the index directory (`dir_index`).
3. Process each sample ID-sample file pair (`fq_files`) in turn:
   1. Cut out sequencing library adapters (`adapters`) using `cutadapt`.
   2. Extract UMIs using `umi_tools extract`, if requested (if `extract_umis: TRUE`), using a UMI-tools-compliant regular expression pattern (`umi_regexp`). The extracted UMIs are inserted into the read headers of the FASTQ records.
   3. Remove rRNA or other contaminating reads by alignment to rRNA index file (`rRNA_index`) using `hisat2`.
   4. Align remaining reads to ORFs index file (`orf_index`). using `hisat2`.
   5. Trim 5' mismatches from reads and remove reads with more than 2 mismatches using `trim_5p_mismatch.py`.
   6. Output UMI groups pre-deduplication using `umi_tools group` if requested (if `dedup_umis: TRUE` and `group_umis: TRUE`)
   7. Deduplicate reads using `umi_tools dedup`, if requested (if `dedup_umis: TRUE`)
   8. Output UMI groups post-deduplication using `umi_tools group` if requested (if `dedup_umis: TRUE` and `group_umis: TRUE`)
   9. Export bedgraph files for plus and minus strands, if requested (if `make_bedgraph: TRUE`) using `bedtools genomecov`.
   10. Write intermediate files produced above into a sample-specific directory, named using the sample ID, within the temporary directory (`dir_tmp`).
   11. Make length-sensitive alignments in compressed h5 format using `bam_to_h5.R`.
   12. Generate summary statistics, and analyses and QC plots for both RPF and mRNA datasets using `generate_stats_figs.R`. This includes estimated read counts, reads per base, and transcripts per million for each ORF in each sample.
   13. Write output files produced above into an sample-specific directory, named using the sample ID, within the output directory (`dir_out`). 
4. Collate TPMs across results, using `collate_tpms.R` and write into output directory (`dir_out`). Only the results from successfully-processed samples are collated.

Images of the workflow with the key steps, inputs and outputs are available:

* [Workflow](./images/workflow.svg) (SVG).
* [Workflow with deduplication](./images/workflow-dedup.svg) (SVG), if `dedup_umis: TRUE`.

---

## Process multiplexed ribosome profiling sample data

If a multiplexed file (`multiplex_fq_files`) is specified, then `prep_riboviz.py`, the `prep_riboviz.py` processes the multiplexed file as follows:

1. Read configuration information (as for 1. above).
2. Build hisat2 indices if requested (as for 2. above).
3. Read the multiplexed FASTQ file (`multiplex_fq_files`).
4. Cut out sequencing library adapters (`adapters`) using `cutadapt`.
5. Extract barcodes and UMIs using `umi_tools extract`, if requested (if `extract_umis: TRUE`), using a UMI-tools-compliant regular expression pattern (`umi_regexp`).
6. Demultiplex file with reference to the sample sheet (`sample_sheet`), using `demultiplex_fastq.py`. Sample IDs in the `SampleID` column in the sample sheet are used to name the demultiplexed files.
7. Process each demultiplexed FASTQ file which has one or more reads, in turn (as for 3.3 to 3.13 above)
8. Collate TPMs across results, using `collate_tpms.R` and write into output directory (`dir_out`) (as for 4. above.

---

## Index files

Index files (HT2) are produced in the index directory (`dir_index`).

---

## Temporary files

Intermediate files are produced within the temporary directory (`dir_tmp`).

For each sample (`<SAMPLE_ID>`), intermediate files are produced in a sample-specific subdirectory (`<SAMPLE_ID>`):

* `<SAMPLE_ID>_trim.fq`: trimmed reads.
* `<SAMPLE_ID>_nonrRNA.fq`: trimmed non-rRNA reads.
* `<SAMPLE_ID>_rRNA_map.sam`: rRNA-mapped reads.
* `<SAMPLE_ID>_orf_map.sam`: ORF-mapped reads.
* `<SAMPLE_ID>_orf_map_clean.sam`: ORF-mapped reads with mismatched nt trimmed.
* `<SAMPLE_ID>_unaligned.sam`: unaligned reads. These files can be used to find common contaminants or translated sequences not in your ORF annotation.

If deduplication is enabled (if `dedup_umis: TRUE`) the following sample-specific files are also produced:

* `<SAMPLE_ID>_extract_trim.fq`: trimmed reads with UMIs extracted.
* UMI groups pre- and post-deduplication (if `group_umis: TRUE`):
  - `<SAMPLE_ID>_pre_dedup_groups.tsv`: UMI groups before deduplication.
  - `<SAMPLE_ID>_post_dedup_groups.tsv`: UMI groups after deduplication.
* UMI deduplication statistics:
  - `<SAMPLE_ID>_dedup_stats_edit_distance.tsv`: edit distance between UMIs at each position.
  - `<SAMPLE_ID>_dedup_stats_per_umi_per_position.tsv`: histogram of counts per position per UMI pre- and post-deduplication.
  - `<SAMPLE_ID>_dedup_stats_per_umi.tsv`: number of times each UMI was observed, total counts and median counts, pre- and post-deduplication
* For more information on the `stats` files, see see UMI-tools [Dedup-specific options](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html) and [documentation on stats file #250](https://github.com/CGATOxford/UMI-tools/issues/250)

If a multiplexed file (`multiplex_fq_files`) is specified, then the following files and directories are also written into the temporary directory:

* `<FASTQ_FILE_NAME_PREFIX>_trim.fq`: FASTQ file post-adapter trimming, where `<FASTQ_FILE_NAME_PREFIX>` is the name of the file (without path or extension) in `multiplex_fq_files`.
* `<FASTQ_FILE_NAME_PREFIX>_extract_trim.fq`: `<FASTQ_FILE_NAME_PREFIX_trim.fq` post-barcode and UMI extraction.
* `<FASTQ_FILE_NAME_PREFIX>_deplex/`: demultiplexing results directory including:
   - `num_reads.tsv`: a tab-separated values file with columns:
     - `SampleID`, copied from the sample sheet.
     - `TagRead` (barcode), coped from the sample sheet.
     - `NumReads`, number of reads detected for each sample.
     - Row with `SampleID` with value `Unassigned` and `NumReads` value with the number of unassigned reads.
     - Row with `SampleID` with value `Total` and `NumReads` value with the total number of reads processed. 
  - `<SAMPLE_ID>.fastq`: Files with demultiplexed reads, where `<SAMPLE_ID>` is a value in the `SampleID` column of the sample sheet. There will be one file per sample.
  - `Unassigned.fastq`: A FASTQ file with the reads that did not match any `TagRead` (barcode) in the sample sheet.

---

## Output files

Output files are produced within the output directory (`dir_out`).

For each sample (`<SAMPLE_ID>`), intermediate files are produced in a sample-specific subdirectory (`<SAMPLE_ID>`):

* `<SAMPLE_ID>.bam`: BAM file of reads mapped to transcripts, which can be directly used in genome browsers.
* `<SAMPLE_ID>.bam.bai`: BAM index file for `<SAMPLE_ID>.bam`.
* `<SAMPLE_ID>_minus.bedgraph`: bedgraph of reads from minus strand (if `make_bedgraph: TRUE`).
* `<SAMPLE_ID>_plus.bedgraph`: bedgraph of reads from plus strand (if `make_bedgraph: TRUE`).
* `<SAMPLE_ID>.h5`: length-sensitive alignments in compressed h5 format.
* `<SAMPLE_ID>_3nt_periodicity.tsv`
* `<SAMPLE_ID>_3nt_periodicity.pdf`
* `<SAMPLE_ID>_read_lengths.tsv`
* `<SAMPLE_ID>_read_lengths.pdf`
* `<SAMPLE_ID>_pos_sp_nt_freq.tsv`
* `<SAMPLE_ID>_pos_sp_rpf_norm_reads.pdf`
* `<SAMPLE_ID>_pos_sp_rpf_norm_reads.tsv`
* `<SAMPLE_ID>_features.pdf`
* `<SAMPLE_ID>_tpms.tsv`
* `<SAMPLE_ID>_codon_ribodens.tsv`
* `<SAMPLE_ID>_codon_ribodens.pdf`
* `<SAMPLE_ID>_startcodon_ribogridbar.pdf`
* `<SAMPLE_ID>_startcodon_ribogrid.pdf`
* `<SAMPLE_ID>_3ntframe_bygene.tsv`
* `<SAMPLE_ID>_3ntframe_propbygene.pdf`

If deduplication is enabled (if `dedup_umis: TRUE`) the following sample-specific files are also produced:

* `<SAMPLE_ID>_dedup.bam`: deduplicated BAM file.
* `<SAMPLE_ID>_dedup.bam.bai`: BAM index file for `<SAMPLE_ID>_dedup.bam`.

A summary file is also put in the output directory (`dir_out`):

* `TPMs_collated.tsv`: file with the transcripts per million (tpm) for all successfully processed samples.

---

## Log files

Information on the execution of `prep_riboviz.py`, including the causes of any errors, is added to a timestamped log file in the current directory, named `riboviz-YYYYMMDD-HHMMSS.log` (for example, `riboviz.20190926-002455.log`).

Log files for each processing step are placed in a timestamped subdirectory (`YYYYMMDD-HHMMSS`) within the logs directory (`dir_logs`). 

For each sample (`<SAMPLE_ID>`), log files are produced in a sample-specific directory (`<SAMPLE_ID>`) within this timestamped subdirectory.

The following log files are produced:

```console
hisat2_build_r_rna.log
hisat2_build_orf.log
<SAMPLE_ID>/
  <SAMPLE_ID>_01_cutadapt.log
  <SAMPLE_ID>_02_hisat2_rrna.log
  <SAMPLE_ID>_03_hisat2_orf.log
  <SAMPLE_ID>_04_trim_5p_mismatch.log
  <SAMPLE_ID>_05_samtools_view_sort.log
  <SAMPLE_ID>_06_samtools_index.log
  <SAMPLE_ID>_07_bedtools_genome_cov_plus.log
  <SAMPLE_ID>_08_bedtools_genome_cov_minus.log
  <SAMPLE_ID>_09_bam_to_h5.log
  <SAMPLE_ID>_10_generate_stats_figs.log
collate_tpms.log
```

If deduplication is enabled (if `dedup_umis: TRUE`), then the following log files are produced:

```
hisat2_build_r_rna.log
hisat2_build_orf.log
<SAMPLE_ID>/
  <SAMPLE_ID>_01_cutadapt.log
  <SAMPLE_ID>_02_umi_tools_extract.log
  <SAMPLE_ID>_03_hisat2_rrna.log
  <SAMPLE_ID>_04_hisat2_orf.log
  <SAMPLE_ID>_05_trim_5p_mismatch.log
  <SAMPLE_ID>_06_samtools_view_sort.log
  <SAMPLE_ID>_07_samtools_index.log
  <SAMPLE_ID>_08_umi_tools_group.log
  <SAMPLE_ID>_09_umi_tools_dedup.log
  <SAMPLE_ID>_10_samtools_index.log
  <SAMPLE_ID>_11_umi_tools_group.log
  <SAMPLE_ID>_12_bedtools_genome_cov_plus.log
  <SAMPLE_ID>_13_bedtools_genome_cov_minus.log
  <SAMPLE_ID>_14_bam_to_h5.log
  <SAMPLE_ID>_15_generate_stats_figs.log
collate_tpms.log
```

If a multiplexed file (`multiplex_fq_files`) specified, then the following log files are produced:

```
hisat2_build_r_rna.log
hisat2_build_orf.log
cutadapt.log
umi_tools_extract.log
demultiplex_fastq.log
<SAMPLE_ID>/
  <SAMPLE_ID>_01_hisat2_rrna.log
  <SAMPLE_ID>_02_hisat2_orf.log
  <SAMPLE_ID>_03_trim_5p_mismatch.log
  <SAMPLE_ID>_04_samtools_view_sort.log
  <SAMPLE_ID>_05_samtools_index.log
  <SAMPLE_ID>_06_umi_tools_group.log
  <SAMPLE_ID>_07_umi_tools_dedup.log
  <SAMPLE_ID>_08_samtools_index.log
  <SAMPLE_ID>_09_umi_tools_group.log
  <SAMPLE_ID>_10_bedtools_genome_cov_plus.log
  <SAMPLE_ID>_11_bedtools_genome_cov_minus.log
  <SAMPLE_ID>_12_bam_to_h5.log
  <SAMPLE_ID>_13_generate_stats_figs.log
collate_tpms.log
```

