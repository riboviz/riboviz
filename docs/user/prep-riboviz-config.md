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

| Parameter | Description |
| --------- | ----------- |
| `adapters` | Illumina sequencing adapter(s) to remove |
| `aligner` | Short read aligner to use (currently ignored - hisat2 is used) |
| `asite_disp_length_file` | Summary of read frame displacement from 5' end to A-site for each read length based on "standard" yeast data from early ribosome profiling papers (tab-separated values file) (optional) |
| `buffer` | Length of flanking region around the CDS |
| `build_indices` | `TRUE` or `FALSE`, rebuild indices from FASTA files? If `FALSE` then `dir_index` is expected to contain the index files |
| `cmd_file` | Bash commands file, to log bash commands executed by the workflow (if omitted, `run_riboviz_vignette.sh` is used) (Python workflow only) |
| `codon_positions_file` | Position of codons within each gene (RData file)  (optional) |
| `count_reads` | `TRUE` or `FALSE`, scan input, temporary and output files and produce counts of reads in each FASTQ, SAM, and BAM file processed? |
| `count_threshold` | Remove genes with a read count below this threshold, when generating statistics and figures |
| `dataset` | Human-readable name of the dataset |
| `dedup_umis` | `TRUE` or `FALSE`, deduplicate reads using UMI-tools? |
| `dir_in` | Input directory |
| `dir_index` | Built indices directory |
| `dir_logs` | Log files directory (Python workflow only) |
| `dir_out` | Output directory |
| `dir_tmp` | Intermediate files directory |
| `do_pos_sp_nt_freq` | `TRUE` or `FALSE`, calculate position-specific nucleotide freqeuency? |
| `extract_umis` | `TRUE` or `FALSE`, extract UMIs after adapter trimming? |
| `features_file` | Features to correlate with ORFs (tab-separated values file)  (optional) |
| `fq_files` |  List of FASTQ files to be processed, relative to `<dir_in>`. Each list member consists of identifier key with a file name value (e.g. `WT3AT: SRR1042864_s1mi.fastq.gz`). | 
| `group_umis` | `TRUE` or `FALSE`, summarise UMI groups both pre- and post-deduplication, using UMI-tools? Useful for debugging. |
| `is_riboviz_gff` | `TRUE` or `FALSE`, does the GFF file contain 3 elements per gene - UTR5, CDS, and UTR3? |
| `is_test_run` | `TRUE` or `FALSE`, is this a test run? (unused) |
| `make_bedgraph` | `TRUE` or `FALSE`, output bedgraph data files in addition to H5 files? |
| `max_read_length` | Maximum read length in H5 output |
| `min_read_length` | Minimum read length in H5 output |
| `multiplex_fq_files` | List with a single multiplexed FASTQ file, relative to `<dir_in>`. If this is provided then the `fq_files` parameter must not be present in the configuration and the `sample_sheet` parameter must be present. |
| `num_processes` | Number of processes to parallelize over, used by specific steps in the workflow |
| `orf_fasta_file` | Transcript sequences file containing both coding regions and flanking regions (FASTA file) |
| `orf_gff_file` | Matched genome feature file, specifying coding sequences locations (start and stop coordinates) within the transcripts (GTF/GFF3 file) |
| `orf_index_prefix` | Prefix for ORF index files, relative to `<dir_index>` |
| `primary_id` | Primary gene IDs to access the data (YAL001C, YAL003W, etc.) |
| `rpf` | `TRUE` or `FALSE`, is the dataset an RPF or mRNA dataset? |
| `rrna_fasta_file` | Ribosomal rRNA and other contaminant sequences to avoid aligning to (FASTA file) |
| `rrna_index_prefix` | Prefix for rRNA index files, relative to `<dir_index>` |
| `sample_sheet` | A sample sheet, relative to `<dir_in>`, mandatory if `multiplex_fq_files` is used (tab-separated values file) |
| `secondary_id` | Secondary gene IDs to access the data (COX1, EFB1, etc. or `NULL`) |
| `stop_in_cds` | `TRUE` or `FALSE`, are stop codons part of the CDS annotations in GFF? |
| `t_rna_file` | tRNA estimates file (tab-separated values file)  (optional) |
| `umi_regexp` | UMI-tools-compliant regular expression to extract barcodes and UMIs. For details on the regular expression format, see UMI-tools documentation on [Barcode extraction](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction). Only required if `extract_umis` is `TRUE`. |

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
