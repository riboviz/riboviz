# Deduplication and demultiplexing

As described in [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md), `prep_riboviz.py` can be configured to do barcode and UMI extraction, using [UMI-tools](https://umi-tools.readthedocs.io/), demultiplexing, using the custom [demultiplex-fastq.py fastq demultiplexer](./demultiplex-fastq.md), and deduplication, using [UMI-tools](https://umi-tools.readthedocs.io/). This page describes deduplication and demultiplexing in more detail.

In the following, `<dir_in>`, `<dir_tmp>`, `<dir_out>` and `<dir_logs>` denote the values of the `dir_in`, `dir_tmp`, `dir_out` and `dir_logs` configuration parameters.

---

## UMI extraction and deduplication

Configuration parameters pertinent to UMI extraction and deduplication are as follows:

* `fq_files`: A dictionary with pairs of sample IDs and FASTQ (`.fastq[.gz]`) files (assumed to be in `<dir_in>`).
* `extract_umis`: If `TRUE` then UMIs will be extracted. The UMIs are extracted using `umi_tools extract` on the trimmed FASTQ file, using the regular expression specified in `umi_regexp` below. The extracted UMIs are inserted into the read headers of the FASTQ records.
* `umi_regexp`: A UMI-tools-compliant regular expression to extract UMIs. For example `^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$` extracts a 4nt UMI from the 5' end of a read and a 4nt UMI from the 3' end. For details on the regular expression format, see UMI-tools documentation on [Barcode extraction](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction).
* `dedup_umis`: If `TRUE` then the reads will be deduplicated using `umi_tools dedup`.
* `group_umis`: If `TRUE` then the UMI groups are summarised both pre- and post-deduplication, using `umi_tools group`. This is provided for debugging purposes.

Note: if `dedup_umis` is `TRUE` but `extract_umis` is `FALSE` then a warning will be displayed, but processing will continue.

After a run, the following UMI extraction and deduplication-specific files will exist, in addition to the other files output during a run of `prep_riboviz.py`:

* UMI groups pre- and post-deduplication (in `<dir_tmp>`):
  - `<SAMPLE>_pre_dedup_groups.tsv`: UMI groups before deduplication
  - `<SAMPLE>_post_dedup_groups.tsv`: UMI groups after deduplication
* UMI deduplication statistics (in `<dir_tmp>`):
 - `<SAMPLE>_dedup_stats_edit_distance.tsv`: edit distance between UMIs at each position.
  - `<SAMPLE>_dedup_stats_per_umi_per_position.tsv`: histogram of counts per position per UMI pre- and post-deduplication.
 - `<SAMPLE>_dedup_stats_per_umi.tsv`: number of times each UMI was observed, total counts and median counts, pre- and post-deduplication
 - For more information see UMI-tools [Dedup-specific options](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html) and [documentation on stats file #250](https://github.com/CGATOxford/UMI-tools/issues/250)
* Log files (in `<dir_logs>/YYYYMMNN-HHMMSS/`):

```
<SAMPLE>_02_umi_tools_extract.log
<SAMPLE>_08_umi_tools_group.log
<SAMPLE>_09_umi_tools_dedup.log
<SAMPLE>_11_umi_tools_group.log
```

### Troubleshooting: `WARNING: dedup_umis was TRUE but extract_umis was FALSE`

This error in the log file means that in your YAML configuration file you have defined:

```yaml
extract_umis: FALSE
dedup_umis: TRUE
```

Unless you explicitly want this you should:

* Either, set `extract_umis` to `TRUE`, if you want UMI deduplication to occur.
* Or, set `dedup_umis` to `FALSE`, if you do not want UMI deduplication to occur.

### Run example

[vignette/example_config.yaml](../vignette/example_config.yaml) has a sample configuration file which runs an analysis of `data/example/umi5_umi3_umi_adaptor.fastq` with UMI extraction and deduplication enabled. For more details on the contents on `data/example/`, see [Simulated FASTQ test files](./data.md#simulated-fastq-test-files).

The example can be run as follows:

```console
$ python -m riboviz.tools.prep_riboviz rscripts/ vignette/example_config.yaml 
```

---

## Barcode and UMI extraction, demultiplexing and deduplication

Configuration parameters pertinent to barcode and UMI extraction, demultiplexing and deduplication are as follows:

* `multiplex_fq_files`: A list with a single multiplexed FASTQ (`fastq[.gz]`) file (assumed to be in `<dir_in>`).
* `sample_sheet`: A sample sheet (assumed to be in `<dir_in>)`. This must be a tab-separated values file with, at least, `SampleID` and `TagRead` (barcode) columns.
* `extract_umis`: as described above.
* `umi_regexp`: as described above, but specifying a regular expression that extracts barcodes too (via the use of `<cell_<N>>` patterns. For example `^(?P<umi_1>.{4}).+(?P<umi_2>.{4})(?P<cell_1>.{3})$` extracts a 3nt barcode from the 3' end of a read then extracts a 4nt UMI from the 5' end and a 4nt UMI from the 3' end.
* `dedup_umis`: as described above.
* `group_umis`: as described above.

Note: if  both `fq_files` and `multiplex_fq_files` parameters are provided then `prep_riboviz.py` will exit. Only a group of non-multiplexed files or a single multiplexed file can be provided.

After a run, the following barcode and UMI extraction, demultiplexing, and deduplication-specific files will exist, in addition to the other files output during a run of `prep_riboviz.py`:

* FASTQ file post-barcode and UMI extraction (in `<dir_tmp>`):
  - `<FASTQ_FILE_NAME_PREFIX>_extract_trim.fq` where `<FASTQ_FILE_NAME_PREFIX>` is the name of the file (without path or extension) in `multiplex_fq_files`.
* Demultiplexing results (in `<dir_tmp>/<FASTQ_FILE_NAME_PREFIX>_deplex`):
  - `num_reads.tsv`: a tab-separated values file with:
     - `SampleID` column, from the sample sheet.
     - `TagRead` (barcode) column, from the sample sheet.
     - `NumReads` column, with the number of reads detected for each sample.
     - Row with `SampleID` with value `Unassigned` and `NumReads` value with the number of unassigned reads.
     - Row with `SampleID` with value `Total` and `NumReads` value with the total number of reads processed. 
  - `<SAMPLE_ID>.fastq`: Files with demultiplexed reads. Each file name has format `<SAMPLE_ID>.fastq` where `<SAMPLE_ID>` is a value in the `SampleID` column of the sample sheet. There will be one file per sample.
  - `Unassigned.fastq`: A FASTQ file with the reads that did not match any `TagRead` (barcode) in the sample sheet.
* Sample-specific temporary files (in `<dir_tmp>/<SAMPLE_ID>`):
  - There will be one directory per file that had one or more reads.
  - This contains the temporary files specific to processing this sample, including those arising from deduplication described above.
* Sample-specific output files (in `<dir_out>/<SAMPLE_ID>`):
  - There will be one directory per file that had one or more reads.
  - This contains the output files specific to processing this sample, including those arising from deduplication described above.
  - This will also contain a `<SAMPLE_ID>_TPMs_collated.tsv` file for the sample.
* Log files (in `<dir_logs>/YYYYMMNN-HHMMSS/`):

```
umi_tools_extract.log
demultiplex_fastq.log

<SAMPLE_ID>_06_umi_tools_group.log
<SAMPLE_ID>_07_umi_tools_dedup.log
<SAMPLE_ID>_09_umi_tools_group.log
<SAMPLE_ID>_14_collate_tpms.log
```

### Run example

[vignette/example_multiplex_config.yaml](../vignette/example_multiplex_config.yaml) has a sample configuration file which runs an analysis of `data/example/multiplex.fastq` with barcode and UMI extraction, demultiplexing and deduplication enabled. For more details on the contents on `data/example/`, see [Simulated FASTQ test files](./data.md#simulated-fastq-test-files).

The example can be run as follows:

```console
$ python -m riboviz.tools.prep_riboviz rscripts/ vignette/example_multiplex_config.yaml 
```

### Demultiplexing and dry runs

When run with demultiplexing enabled, `prep_riboviz.py` determines the demultiplexed files to process based on the samples in the `num_reads.tsv` file output by `demultipled_fastq.py` which have one or more reads.

If running `prep_riboviz.py` in dry run mode (with the `--dry-run` flag) then this file does not exist. In this case, `prep_riboviz.py` will determine the samples based on the samples in the sample sheet.

This means that the bash script output by `prep_riboviz.py` with `--dry-run` may not match that output during a live run, if, during the live run, some samples had 0 associated reads.
