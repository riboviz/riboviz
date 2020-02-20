# Configuring the RiboViz workflow

This page describes the inputs that `riboviz.tools.prep_riboviz` requires and how it is configured.

---

## Input files

`prep_riboviz` requires the following inputs (file formats are in brackets).

### Configuration

A configuration file containing all the information required for the analysis workflow, including the locations of the other input files ([YAML](http://www.yaml.org/)). For details, see [Configuration parameters](#configuration-parameters) below.

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

`prep_riboviz` supports the following configuration parameters. All directory and file paths can be relative or absolute. If relative then they are relative to the directory `prep_riboviz` is invoked from.

| Parameter | Description |
| --------- | ----------- |
| `dir_in` | Input directory |
| `dir_out` | Output directory |
| `dir_tmp` | Intermediate files directory |
| `dir_logs` | Log files directory |
| `cmd_file` | Bash commands file, to log bash commands executed by the workflow (if omitted, `run_riboviz_vignette.sh` is used) |
| `rrna_fasta_file` | Ribosomal rRNA and other contaminant sequences to avoid aligning to (FASTA file) |
| `orf_fasta_file` | Transcript sequences file containing both coding regions and flanking regions (FASTA file) |
| `orf_gff_file` | Matched genome feature file, specifying coding sequences locations (start and stop coordinates) within the transcripts (GTF/GFF3 file) |
| `build_indices` | `TRUE|FALSE`, rebuild indices from FASTA files? |
| `dir_index` | Built indices directory |
| `rrna_index_prefix` | Prefix for rRNA index files, relative to `<dir_index>` |
| `orf_index_prefix` | Prefix for ORF index files, relative to `<dir_index>` |
| `fq_files` |  List of FASTQ files to be processed, relative to `<dir_in>`. Each list member consists of identifier key with a file name value (e.g. `WT3AT: SRR1042864_s1mi.fastq.gz`). | 
| `multiplex_fq_files` | List with a single multiplexed FASTQ file, relative to `<dir_in>`. If this is provided then the `fq_files` parameter must not be present in the configuration and the `sample_sheet` parameter must be present. |
| `adapters` | Illumina sequencing adapter(s) to remove |
| `extract_umis` | `TRUE|FALSE`, extract UMIs after adapter trimming? |
| `umi_regexp` | UMI-tools-compliant regular expression to extract barcodes and UMIs. For details on the regular expression format, see UMI-tools documentation on [Barcode extraction](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction). Only required if `extract_umis` is `TRUE`. |
| `dedup_umis` | `TRUE|FALSE`, deduplicate reads using UMI-tools? |
| `group_umis` | `TRUE|FALSE`, summarise UMI groups both pre- and post-deduplication, using UMI-tools? Useful for debugging. |
| `sample_sheet` | A sample sheet, relative to `<dir_in>`, mandatory if `multiplex_fq_files` is used (tab-separated values file) |
| `num_processes` | Number of processes to parallelize over, used by specific steps in the workflow |
| `aligner` | Short read aligner to use (currently ignored - hisat2 is used) |
| `make_bedgraph` | `TRUE|FALSE`, output bedgraph data files in addition to H5 files? |
| `min_read_length` | Minimum read length in H5 output |
| `max_read_length` | Maximum read length in H5 output |
| `buffer` | Length of flanking region around the CDS |
| `primary_id` | Primary gene IDs to access the data (YAL001C, YAL003W, etc.) |
| `secondary_id` | Secondary gene IDs to access the data (COX1, EFB1, etc. or `NULL`) |
| `dataset` | Human-readable name of the dataset |
| `stop_in_cds` | `TRUE|FALSE`, are stop codons part of the CDS annotations in GFF? |
| `is_test_run` | `TRUE|FALSE`, is this a test run? (unused) |
| `rpf` | `TRUE|FALSE`, is the dataset an RPF or mRNA dataset? |
| `is_riboviz_gff` | `TRUE|FALSE`, does the GFF file contain 3 elements per gene - UTR5, CDS, and UTR3? |
| `features_file` | Features to correlate with ORFs (tab-separated values file) |
| `do_pos_sp_nt_freq` | `TRUE|FALSE`, calculate position-specific nucleotide freqeuency? |
| `t_rna_file` | tRNA estimates file (tab-separated values file) |
| `codon_positions_file` | Position of codons within each gene (RData file) |
| `count_threshold` | Remove genes with a read count below this threshold, when generating statistics and figures |
| `asite_disp_length_file` | Summary of read frame displacement from 5' end to A-site for each read length based on "standard" yeast data from early ribosome profiling papers (tab-separated values file) |
| `count_reads` | `TRUE|FALSE`, scan input, temporary and output files and produce counts of reads in each FASTQ, SAM, and BAM file processed? |

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

If both `fq_files` and `multiplex_fq_files` parameters are provided then `prep_riboviz` will exit. Only a group of non-multiplexed files or a single multiplexed file can be provided.

If `fq_files` are provided then `umi_regexp` should extract only UMIs (i.e. it should contain `<umi>` elements only).

If `multiplex_fq_files` is provided then `umi_regexp` should extract both barcodes and UMIs (i.e. it should contain both `<cell>` and `<umi>` elements).
