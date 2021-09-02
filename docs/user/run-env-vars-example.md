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
$ cp riboviz/vignette/vignette_config.yaml ribo-run/ribo_config.yaml
$ cd ribo-run
```

Edit `ribo_config.yaml` and update the following parameters to the values shown:

```yaml
asite_disp_length_file: ${RIBOVIZ_DATA}/yeast_standard_asite_disp_length.txt
codon_positions_file: ${RIBOVIZ_DATA}/yeast_codon_pos_i200.RData
dir_index: ${RIBOVIZ_SAMPLES}/index
dir_in: ${RIBOVIZ_SAMPLES}/input
dir_out: ${RIBOVIZ_SAMPLES}/output
dir_tmp: ${RIBOVIZ_SAMPLES}/tmp
features_file: ${RIBOVIZ_DATA}/yeast_features.tsv
orf_gff_file: ${RIBOVIZ_ORGANISMS}/yeast_YAL_CDS_w_250utrs.gff3
orf_fasta_file: ${RIBOVIZ_ORGANISMS}/yeast_YAL_CDS_w_250utrs.fa	
rrna_fasta_file: ${RIBOVIZ_ORGANISMS}/yeast_rRNA_R64-1-1.fa
t_rna_file: ${RIBOVIZ_DATA}/yeast_tRNAs.tsv
```

Run RiboViz, specifying values for the configuration tokens via environment variables:

```console
$ RIBOVIZ_DATA=$HOME/ribo-data \
  RIBOVIZ_ORGANISMS=$HOME/ribo-organisms \
  RIBOVIZ_SAMPLES=$HOME/ribo-samples \
  nextflow run $HOME/riboviz/prep_riboviz.nf  -ansi-log false -params-file ribo_config.yaml
```

Optionally, run integration tests:

```console
$ RIBOVIZ_DATA=$HOME/ribo-data \
  RIBOVIZ_ORGANISMS=$HOME/ribo-organisms \
  RIBOVIZ_SAMPLES=$HOME/ribo-samples \
  pytest -vs $HOME/riboviz/riboviz/test/integration/test_integration.py \
  --expected=$HOME/test-data-2.1 --skip-workflow --config-file ribo_config.yaml
```
