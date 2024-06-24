# riboviz command-line tools

* [Python command-line tools](#python-command-line-tools)
* [R command-line tools](#r-command-line-tools)

---

## Python command-line tools

The following command-line tools are invoked in the riboviz workflow:

| Tool | Description |
| ---- | ----------- |
| `count_reads` | Scan input, temporary and output directories and count the number of reads (sequences) processed by specific stages of a workflow |
| `demultiplex_fastq` | Demultiplex FASTQ files using UMI-tools-compliant barcodes present within the FASTQ headers and a sample sheet file |
| `trim_5p_mismatch` | Remove a single 5' mismatched nt and filter reads with more than a specified mismatches from a SAM file and save the trimming summary to a file |

The following additional command-line tools are also available:

| Tool | Description |
| ---- | ----------- |
| `check_fasta_gff` | [Check FASTA and GFF files for coding sequence (CDS) features](./check-fasta-gff.md) |
| `create_barcode_pairs` | Create barcode pairs and write each pair plus the Hamming distance between then to a file of tab-separated values |
| `create_fastq_simdata` | Create simulated FASTQ files to test UMI/deduplication, adaptor trimming, and demultiplexing. Files in `data/simdata/` were created using this tool |
| `create_job_script` | [Create job submission script from template](./create-job-script.md) |
| `get_cds_codons` | Extract coding sequence codons and export as a tab-separated values file |
| `subsample_bioseqfile` | Subsample an input FASTQ (or other sequencing) file, to produce a smaller file whose reads are randomly sampled from of the input with a fixed probability |
| `upgrade_config_file]` | [Upgrade configuration files to current version](./upgrade-config.md) |

For usage, run one of:

```console
$ python -m riboviz.tools.<SCRIPT> -h
```
```console
$ <SCRIPT> -h
```

For example:

```console
$ python -m riboviz.tools.check_fasta_gff -h
```
```console
$ check_fasta_gff -h
```

---

## R command-line tools

The following command-line tools are invoked in the riboviz workflow:

| Tool | Description |
| ---- | ----------- |
| `rscripts/bam_to_h5.R` | Given a GFF file and a BAM file, create a riboviz H5 file with information about a feature |
| `rscripts/collate_tpms.R` | Collate sample-specific TPMs files into a single file with data for all the samples |
| `rscripts/generate_stats_figs.R` | Create riboviz data files and image files |

For usage, run:

```console
$ Rscript --vanilla <SCRIPT> -h
```

For example:

```console
$ Rscript --vanilla rscripts/bam_to_h5.R -h
```
