# Run an example of a configuration using environment variable tokens

[Environment variables and configuration tokens](./prep-riboviz-config.md#environment-variables-and-configuration-tokens) describes how the use of configuration tokens representing environment variables with their paths, and the use of symbolic links are supported. This page walks you through a runnable example, using the "vignette" dataset from [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md).

Create and populate data, organisms, and samples directories:

```console
$ cp -r riboviz/data/ ribo-data
$ mkdir ribo-organisms
$ cp riboviz/vignette/input/yeast_* ribo-organisms/
$ mkdir -p ribo-samples/input
$ cp riboviz/vignette/input/SRR10428* ribo-samples/input/
```

Create a directory with a configuration file:

```console
$ mkdir ribo-run/
$ cp riboviz/vignette/riboviz_env.yaml ribo-run
$ cd ribo-run
```

`riboviz_env.yaml` is a sample configuration, based on `vignette_config.yaml` but with tokens for paths to the input files expected by riboviz.

Run riboviz, specifying values for the configuration tokens via environment variables:

```console
$ RIBOVIZ_DATA=$HOME/ribo-data \
  RIBOVIZ_ORGANISMS=$HOME/ribo-organisms \
  RIBOVIZ_SAMPLES=$HOME/ribo-samples \
  nextflow run $HOME/riboviz/prep_riboviz.nf  -ansi-log false -params-file riboviz_env.yaml
```

Optionally, run integration tests:

```console
$ RIBOVIZ_DATA=$HOME/ribo-data \
  RIBOVIZ_ORGANISMS=$HOME/ribo-organisms \
  RIBOVIZ_SAMPLES=$HOME/ribo-samples \
  pytest -vs $HOME/riboviz/riboviz/test/integration/test_integration.py \
  --expected=$HOME/test-data-2.2 --skip-workflow --config-file riboviz_env.yaml
```
