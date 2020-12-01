# Configuring the RiboViz workflow

This page describes the inputs that the RiboViz workflow requires and how it is configured.

These inputs apply to both the Python workflow and the Nextflow workflow.

---

## Input files

The workflow requires the following inputs (file formats are in brackets).

### Configuration

A configuration file containing all the information required for the workflow, including the locations of the other input files ([YAML](http://www.yaml.org/)). For details, see [Configuration parameters](#configuration-parameters) below.

### Organism data

Transcript sequences containing both coding regions and flanking regions (which could be fixed, or coincident with measured UTRs) (FASTA).

Matched genome feature file  specifying coding sequences locations (start and stop coordinates) within the transcripts (GTF/GFF3).

Ribosomal rRNA and other contaminant sequences to avoid aligning to (FASTA).

Additional organism-specific data, for creating statistics and figures:

* Position of codons within each gene (RData).
* Features to correlate with ORFs (tab-separated values, with `ORF`, `Length_log10`, `uATGs`, `FE_atg`, `FE_cap`, `utr`, `utr_gc`, `polyA` columns).
* tRNA estimates. (tab-separated values, with `AA`, `Codon`, `tRNA`, `tAI`, `Microarray`, `RNA.seq` columns).
* Summary of read frame displacement from 5' end to A-site for each read length based on "standard" yeast data from early ribosome profiling papers (tab-separated values, with `read_length`, `asite_disp` columns).

Ribosome profiling data - one or more files (FASTQ).

### Sample sheet

For processing multiplexed FASTQ files, a sample sheet (tab-separated values with, at least, `SampleID` and `TagRead` (barcode) columns).

---

## Configuration parameters

The workflow also supports the following configuration parameters. All directory and file paths can be relative or absolute. If relative then they are relative to the directory the workflow is invoked from.

| Parameter | Description | Mandatory | Default |
| --------- | ----------- | --------- | ------- |
| `adapters` | Illumina sequencing adapter(s) to remove | Yes | |
| `aligner` | Short read aligner to use (currently ignored, hisat2 is used) | No | |
| `asite_disp_length_file` | Summary of read frame displacement from 5' end to A-site for each read length based on "standard" yeast data from early ribosome profiling papers (tab-separated values file with `read_length`, `asite_disp` columns) | No | |
| `buffer` | Length of flanking region around the CDS | No | `250` |
| `build_indices` | Rebuild indices from FASTA files? If `false` then `dir_index` is expected to contain the index files | No | `true` |
| `cmd_file` | Bash commands file, to log bash commands executed by the workflow (Python workflow only) | No | `run_riboviz_vignette.sh` |
| `codon_positions_file` | Position of codons within each gene (RData file) | Only if `t_rna_file` is also provided | |
| `count_reads` | Scan input, temporary and output files and produce counts of reads in each FASTQ, SAM, and BAM file processed? | No | `true` |
| `count_threshold` | Remove genes with a read count below this threshold, when generating statistics and figures | No | `1` |
| `dataset` | Human-readable name of the dataset | No | `dataset` |
| `dedup_stats` | Output UMI deduplication statistics? | No | `true` |
| `dedup_umis` | Deduplicate reads using UMI-tools? | No | `false` |
| `dir_in` | Input directory | Yes | |
| `dir_index` | Built indices directory | No | `index` |
| `dir_logs` | Log files directory (Python workflow only) | Yes | |
| `dir_out` | Output directory | No | `output` |
| `dir_tmp` | Intermediate files directory | No | `tmp` |
| `do_pos_sp_nt_freq` | Calculate position-specific nucleotide freqeuency? | No | `true` |
| `extract_umis` | Extract UMIs after adapter trimming? | No | `false` |
| `features_file` | Features to correlate with ORFs (tab-separated values file) | No | |
| `fq_files` |  List of FASTQ files to be processed, relative to `<dir_in>`. Each list member consists of identifier key with a file name value (e.g. `WT3AT: SRR1042864_s1mi.fastq.gz`). | Only if `multiplex_fq_files` is not provided | |
| `group_umis` | Summarise UMI groups both pre- and post-deduplication, using UMI-tools? Useful for debugging. | No | `false` |
| `is_riboviz_gff` | Does the GFF file contain 3 elements per gene - UTR5, CDS, and UTR3? | No | `true` |
| `is_test_run` | Is this a test run? (unused) | No | |
| `job_email_events` | Events triggering emails about batch job. Any combination of `b`(begin), `e` (end), `a` (abort), `s` (suspend). (see [Create job submission script from template](./docs/user/create-job-script.md)) | No | `beas` |
| `job_email` | E-mail address for batch job events (see [Create job submission script from template](./docs/user/create-job-script.md)) | No | `null` |
| `job_memory` | Requested memory for batch job (see [Create job submission script from template](./docs/user/create-job-script.md)) | No | `8GB` |
| `job_name` | Name of batch job (see [Create job submission script from template](./docs/user/create-job-script.md)) | No | `riboviz` |
| `job_num_cpus` | Requested number of CPUs for batch job (see [Create job submission script from template](./docs/user/create-job-script.md)) | No | `4` |
| `job_runtime` | Maximum runtime for batch job (see [Create job submission script from template](./docs/user/create-job-script.md)) | No | `48:00:00` |
| `make_bedgraph` | Output bedgraph data files in addition to H5 files? | No | `true` |
| `max_read_length` | Maximum read length in H5 output | No | `50` |
| `min_read_length` | Minimum read length in H5 output | No | `10` |
| `multiplex_fq_files` | List with a single multiplexed FASTQ file, relative to `<dir_in>`. If this is provided then the `fq_files` parameter must not be present in the configuration and the `sample_sheet` parameter must be present. | Only if `fq_files` is not provided | |
| `nextflow_report_file` | Nextflow report file (see [Create job submission script from template](./docs/user/create-job-script.md)) (Nextflow workflow only) | No | `nextflow-report.html` |
| `nextflow_work_dir` | Nextflow work directory (see [Create job submission script from template](./docs/user/create-job-script.md)) (Nextflow workflow only) | No | `work` |
| `num_processes` | Number of processes to parallelize over, used by specific steps in the workflow | No | `1` |
| `orf_fasta_file` | Transcript sequences file containing both coding regions and flanking regions (FASTA file) | Yes | |
| `orf_gff_file` | Matched genome feature file, specifying coding sequences locations (start and stop coordinates) within the transcripts (GTF/GFF3 file) | Yes | |
| `orf_index_prefix` | Prefix for ORF index files, relative to `<dir_index>` | Yes | |
| `primary_id` | Primary gene IDs to access the data (YAL001C, YAL003W, etc.) | No | `Name` |
| `publish_index_tmp` | Publish index and temporary files to `<dir_index>` and `<dir_tmp>`? If `true` copy index and temporary files from Nextflow's `work/` directory, else use symbolic links only. (Nextflow workflow only - see [Nextflow `work/` directory](../user/prep-riboviz-operation.md#nextflow-work-directory) | No | `false` |
| `rpf` | Is the dataset an RPF or mRNA dataset? | No | `true` |
| `rrna_fasta_file` | Ribosomal rRNA and other contaminant sequences to avoid aligning to (FASTA file) | Yes | |
| `rrna_index_prefix` | Prefix for rRNA index files, relative to `<dir_index>` | Yes | |
| `samsort_memory` | Memory to give to `samtools sort` (Nextflow workflow only) | No | `null` (`samtools sort` uses built-in default `768M`, [samtools sort](http://www.htslib.org/doc/samtools-sort.html) |
 | `sample_sheet` | A sample sheet, relative to `<dir_in>` (tab-separated values file) | Only if `multiplex_fq_files` is provided | |
| `secondary_id` | Secondary gene IDs to access the data (COX1, EFB1, etc. or `NULL`) | No | `NULL` |
| `skip_inputs` | When validating configuration (see `validate_only` below) skip checks for existence of ribosome profiling data files (`fq_files`, `multiplexed_fq_files`, `sample_sheet`)? (Nextflow workflow only) | No | `false` |
| `stop_in_cds` | Are stop codons part of the CDS annotations in GFF? | No | `false` |
| `trim_5p_mismatches` | Trim mismatched 5' base? (Nextflow workflow only) | No | `true` |
| `t_rna_file` | tRNA estimates file (tab-separated values file) | Only if `codon_positions_file` is also provided | |
| `umi_regexp` | UMI-tools-compliant regular expression to extract barcodes and UMIs. For details on the regular expression format, see UMI-tools documentation on [Barcode extraction](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction) | Only if `extract_umis` is `true` | |
| `validate_only ` | Validate configuration, check that mandatory parameters have been provided and that input files exist, then exit without running the workflow? (Nextflow workflow only) | No | `false` |

### Examples

An example `fq_files` list is:

```yaml
fq_files:
  WTnone: SRR1042855_s1mi.fastq.gz
  WT3AT: SRR1042864_s1mi.fastq.gz
```

An example `multiplex_fq_files` list is:

```
multiplex_fq_files:
- multiplex_umi_barcode_adaptor.fastq
```

Example `umi_regexp` are:

* `^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$` extracts a 4nt UMI from the 5' end of a read and a 4nt UMI from the 3' end.
* `^(?P<umi_1>.{4}).+(?P<umi_2>.{4})(?P<cell_1>.{3})$` extracts a 3nt barcode from the 3' end of a read then extracts a 4nt UMI from the 5' end and a 4nt UMI from the 3' end.

### Constraints

If `dedup_umis` is `TRUE` but `extract_umis` is `FALSE` then a warning will be displayed, but processing will continue.

If both `fq_files` and `multiplex_fq_files` parameters are provided then the workflow will exit. Only a group of non-multiplexed files or a single multiplexed file can be provided.

If `fq_files` are provided then `umi_regexp` should extract only UMIs (i.e. it should contain `<umi>` elements only).

If `multiplex_fq_files` is provided then `umi_regexp` should extract both barcodes and UMIs (i.e. it should contain both `<cell>` and `<umi>` elements).

While both `codon_positions_file` and `t_rna_file` are optional, either both must be specified or neither must be specified.
