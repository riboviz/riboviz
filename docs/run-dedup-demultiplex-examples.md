# Run UMI extraction, deduplication and demultiplexing examples

Examples of using `prep_riboviz.py` to do UMI extraction and deduplication and also demultiplexing of multiplexed data are provided. These use [Simulated FASTQ test files](./data.md#simulated-fastq-test-files) located within `data/simdata/`.

## UMI extraction and deduplication

`vignette/simdata_umi_config.yaml` has a sample configuration file which runs an analysis of `data/simdata/umi5_umi3_umi_adaptor.fastq` with UMI extraction and deduplication enabled.

The example can be run as follows:

```console
$ python -m riboviz.tools.prep_riboviz rscripts/ vignette/simdata_umi_config.yaml 
```

## Barcode and UMI extraction, demultiplexing and deduplication

`vignette/simdata_multiplex_config.yaml` has a sample configuration file which runs an analysis of `data/simdata/multiplex.fastq` with barcode and UMI extraction, demultiplexing and deduplication enabled.

The example can be run as follows:

```console
$ python -m riboviz.tools.prep_riboviz rscripts/ vignette/simdata_multiplex_config.yaml 
```
