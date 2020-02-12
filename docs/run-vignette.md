# Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file

`prep_riboviz.py` is a Python implementation of the RiboViz analysis workflow. This page describes how you can run a "vignette" of the workflow on a sample data set - *Saccharomyces cerevisiae* reads - to output of HDF5 files and generate summary statistics.

---

## Input files

The vignette uses the following input files.

Configuration:

* `vignette/vignette_config.yaml`: all the information required for the analysis workflow, including the locations of the other input files.

Downsampled genome and annotation data for *Saccharomyces cerevisiae* (yeast):

* `vignette/input/yeast_YAL_CDS_w_250utrs.fa`: transcript sequences containing both coding regions and flanking regions from just the left arm of chromosome 1.
* `vignette/input/yeast_YAL_CDS_w_250utrs.gff3`: matched genome feature file specifying coding sequences locations (start and stop coordinates) within the transcripts.
* For information on the provenance of these files see [Saccharomyces cerevisiae (yeast) genome and annotation data](./data.md#saccharomyces-cerevisiae-yeast-genome-and-annotation-data) and [Downsampled Saccharomyces cerevisiae (yeast) genome and annotation data](./data.md#downsampled-saccharomyces-cerevisiae-yeast-genome-and-annotation-data).

Ribosomal rRNA and other contaminant sequences to avoid aligning to:

* `vignette/input/yeast_rRNA_R64-1-1.fa`. For information on the provenance of this file see [Ribosomal RNA (rRNA) contaminants to remove](./data.md#ribosomal-rna-rrna-contaminants-to-remove)

Additional organism-specific data, for creating statistics and figures:

* `data/yeast_codon_pos_i200.RData`: position of codons within each gene (the numbering ignores the first 200 codons).
* `data/yeast_features.tsv`: features to correlate with ORFs.
* `data/yeast_tRNAs.tsv`: tRNA estimates.
* `data/yeast_standard_asite_disp_length.txt`: summary of read frame displacement from 5' end to A-site for each read length based on "standard" yeast data from early ribosome profiling papers.
* For information on the provenance of these files see [Additional yeast-specific data](./data.md#additional-yeast-specific-data).

Downsampled ribosome profiling data for *Saccharomyces cerevisiae* (yeast):

* `vignette/input/SRR1042855_s1mi.fastq.gz`: ~1mi-sampled RPFs wild-type no additive.
* `vignette/input/SRR1042864_s1mi.fastq.gz`: ~1mi-sampled RPFs wild-type + 3-AT.
* For information on the provenance of these files see [Downsampled ribosome profiling data from Saccharomyces cerevisiae](./data.md#downsampled-ribosome-profiling-data-from-saccharomyces-cerevisiae)

For full information on how to configure `prep_riboviz.py` and its inputs, see [Configuring the RiboViz workflow](./prep-riboviz-config.md).

---

## Set up your environment

If you have not already done so, activate your Python environment:

* Miniconda Python 3.6+:

```console
$ source $HOME/miniconda3/bin/activate
```

If you have not already done so, set the paths to Hisat2 and Bowtie:

* If you followed [Create `setenv.sh` to configure Hisat2 and Bowtie paths](./install.md#create-setenvsh-to-configure-hisat2-and-bowtie-paths), then run:

```console
$ source $HOME/setenv.sh
```

* Otherwise, run the following (your directory names may be different, depending on the versions of Hisat2 and Bowtie you have):

```console
$ export PATH=~/hisat2-2.1.0:$PATH
$ export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
```

---

## Configure number of processes (optional)

`prep_riboviz.py` can instruct the tools it invokes to use a specific number of processes. By default this is 1.

To configure `prep_riboviz.py` to use additional processes, in `vignette/vignette_config.yaml` change:

```yaml
num_processes: 1
```

* to the desired number of processes:

```yaml
num_processes: 4
```

---

## Dry run `prep_riboviz.py`

`prep_riboviz.py` supports a `--dry-run` command-line parameter which can be used to validate the configuration. 

**Tip:** we strongly recommend doing a dry run before doing a live run on data you have not processed before.

Run `prep_riboviz.py` with `--dry-run` enabled:

* Either:

```console
$ python -m riboviz.tools.prep_riboviz --dry-run rscripts/ vignette/vignette_config.yaml
```

* Or:

```console
$ PYTHONPATH=. python riboviz/tools/prep_riboviz.py --dry-run rscripts/ \
    vignette/vignette_config.yaml
```

where:

* `riboviz/tools/prep_riboviz.py`: path to `prep_riboviz.py`, relative to the RiboViz home directory.
* `--dry-run`: flag to enable the dry run.
* `rscripts/`: path to the directory with RiboViz's R scripts, relative to the RiboViz home directory.
* `vignette/vignette_config.yaml`: path to the vignette configuration file.

### Troubleshooting: `This script needs to be run under Python 3`

This warning arises if you try and run `prep_riboviz.py` under Python 2. You can only run `prep_riboviz.py` with Python 3.

### Troubleshooting: `File not found: vignette/input/example_missing_file.fastq.gz`

This warning is expected and can be ignored. The vignette configuration file intentionally includes a reference to a missing input file. This demonstrates that if processing of one sample fails it does not block the processing of the other samples.

---

## Run `prep_riboviz.py`

Run `prep_riboviz.py`:

* Either:

```console
$ python -m riboviz.tools.prep_riboviz rscripts/ vignette/vignette_config.yaml
```

* Or:

```console
$ PYTHONPATH=. python riboviz/tools/prep_riboviz.py rscripts/ \
    vignette/vignette_config.yaml
```

For full information on how to run `prep_riboviz.py` and the options available, see [Running the RiboViz workflow](./prep-riboviz-running.md)

---

## Check the expected files

You can expect to see the following files produced.

Index files in `vignette/index`:

```
YAL_CDS_w_250.1.ht2
YAL_CDS_w_250.2.ht2
YAL_CDS_w_250.3.ht2
YAL_CDS_w_250.4.ht2
YAL_CDS_w_250.5.ht2
YAL_CDS_w_250.6.ht2
YAL_CDS_w_250.7.ht2
YAL_CDS_w_250.8.ht2
yeast_rRNA.1.ht2
yeast_rRNA.2.ht2
yeast_rRNA.3.ht2
yeast_rRNA.4.ht2
yeast_rRNA.5.ht2
yeast_rRNA.6.ht2
yeast_rRNA.7.ht2
yeast_rRNA.8.ht2
```

Intermediate outputs in `vignette/tmp`:

```
WT3AT/
  nonrRNA.fq
  orf_map_clean.sam
  orf_map.sam
  rRNA_map.sam
  trim.fq
  unaligned.sam
WTnone/
  nonrRNA.fq
  orf_map_clean.sam
  orf_map.sam
  rRNA_map.sam
  trim.fq
  unaligned.sam
```

Outputs in `vignette/output`:

```
WT3AT/
  3ntframe_bygene.tsv
  3ntframe_propbygene.pdf
  3nt_periodicity.pdf
  3nt_periodicity.tsv
  WT3AT.bam
  WT3AT.bam.bai
  codon_ribodens.pdf
  codon_ribodens.tsv
  features.pdf
  WT3AT.h5
  minus.bedgraph
  plus.bedgraph
  pos_sp_nt_freq.tsv
  pos_sp_rpf_norm_reads.pdf
  pos_sp_rpf_norm_reads.tsv
  read_lengths.pdf
  read_lengths.tsv
  startcodon_ribogridbar.pdf
  startcodon_ribogrid.pdf
  tpms.tsv
WTnone/
  3ntframe_bygene.tsv
  3ntframe_propbygene.pdf
  3nt_periodicity.pdf
  3nt_periodicity.tsv
  WTnone.bam
  WTnone.bam.bai
  codon_ribodens.pdf
  codon_ribodens.tsv
  features.pdf
  WTnone.h5
  minus.bedgraph
  plus.bedgraph
  pos_sp_nt_freq.tsv
  pos_sp_rpf_norm_reads.pdf
  pos_sp_rpf_norm_reads.tsv
  read_lengths.pdf
  read_lengths.tsv
  startcodon_ribogridbar.pdf
  startcodon_ribogrid.pdf
  tpms.tsv
TPMs_collated.tsv
```

For full information on what `prep_riboviz.py` does and the files it outputs, see [What the RiboViz workflow does](./prep-riboviz-operation.md).

---

## Cleaning up to run again

Before rerunning the vignette, delete the auto-generated index, temporary, logs and output directories:

```console
$ rm -rf vignette/index
$ rm -rf vignette/logs
$ rm -rf vignette/tmp
$ rm -rf vignette/output
$ rm -rf vignette/logs
```

You might also want to do this if you have run the vignette with a missing R package, and then want to run it again from scratch.

---

## Customising the vignette

We suggest copying `vignette/vignette_config.yaml` and the rest of the `vignette/` directory, and then customising it for use with your own datasets.
