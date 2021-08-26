# Run UMI extraction, deduplication and demultiplexing examples

Examples of configurations that request that the RiboViz workflow do UMI extraction and deduplication and also demultiplexing of multiplexed data are provided. These use [Simulated FASTQ test files](../reference/data.md#simulated-fastq-test-files) located within `data/simdata/`.

**Note:** These currently do not output `codon_ribodens.pdf|tsv` during the invocation of `generate_stats_figs.R`. This is a tiny mock data with very few reads and failures arise if these files are requested.

## UMI extraction and deduplication

`vignette/simdata_umi_config.yaml` has a sample configuration file which runs an analysis of `data/simdata/umi5_umi3_umi_adaptor.fastq` with UMI extraction and deduplication enabled.

To run this example using the Nextflow workflow:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/simdata_umi_config.yaml -ansi-log false
```

## Barcode and UMI extraction, demultiplexing and deduplication

`vignette/simdata_multiplex_config.yaml` has a sample configuration file which runs an analysis of `data/simdata/multiplex.fastq` with barcode and UMI extraction, demultiplexing and deduplication enabled.

To run this example using the Nextflow workflow:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/simdata_multiplex_config.yaml -ansi-log false
```
