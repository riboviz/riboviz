# Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file

This page describes how you can run a "vignette" of the workflow on a sample data set - *Saccharomyces cerevisiae* reads - to output of HDF5 files and generate summary statistics.

---

## Input files

The vignette uses the following input files.

Configuration:

* `vignette/vignette_config.yaml`: all the information required for the analysis workflow, including the locations of the other input files.

Downsampled genome and annotation data for *Saccharomyces cerevisiae* (yeast):

* `vignette/input/yeast_YAL_CDS_w_250utrs.fa`: transcript sequences containing both coding regions and flanking regions from just the left arm of chromosome 1.
* `vignette/input/yeast_YAL_CDS_w_250utrs.gff3`: matched genome feature file specifying coding sequences locations (start and stop coordinates) within the transcripts.
* For information on the provenance of these files see [Saccharomyces cerevisiae (yeast) genome and annotation data](../reference/data.md#saccharomyces-cerevisiae-yeast-genome-and-annotation-data) and [Downsampled Saccharomyces cerevisiae (yeast) genome and annotation data](../reference/data.md#downsampled-saccharomyces-cerevisiae-yeast-genome-and-annotation-data).

Ribosomal rRNA and other contaminant sequences to avoid aligning to:

* `vignette/input/yeast_rRNA_R64-1-1.fa`. For information on the provenance of this file see [Ribosomal RNA (rRNA) contaminants to remove](../reference/data.md#ribosomal-rna-rrna-contaminants-to-remove)

Additional organism-specific data, for creating statistics and figures:

* `data/yeast_codon_pos_i200.RData`: position of codons within each gene (the numbering ignores the first 200 codons).
* `data/yeast_features.tsv`: features to correlate with ORFs.
* `data/yeast_tRNAs.tsv`: tRNA estimates.
* `data/yeast_standard_asite_disp_length.txt`: summary of read frame displacement from 5' end to A-site for each read length based on "standard" yeast data from early ribosome profiling papers.
* For information on the provenance of these files see [Additional yeast-specific data](../reference/data.md#additional-yeast-specific-data).

Downsampled ribosome profiling data for *Saccharomyces cerevisiae* (yeast):

* `vignette/input/SRR1042855_s1mi.fastq.gz`: ~1mi-sampled RPFs wild-type no additive.
* `vignette/input/SRR1042864_s1mi.fastq.gz`: ~1mi-sampled RPFs wild-type + 3-AT.
* For information on the provenance of these files see [Downsampled ribosome profiling data from Saccharomyces cerevisiae](../reference/data.md#downsampled-ribosome-profiling-data-from-saccharomyces-cerevisiae)

For full information on how to configure the workflow and its inputs, see [Configuring the RiboViz workflow](./prep-riboviz-config.md).

---

## Set up your environment

If you have not already done so, activate your Python environment:

```console
$ source $HOME/miniconda3/bin/activate
$ conda activate riboviz
```

If you have not already done so, set the paths to Hisat2 and Bowtie:

* If you followed [Create `set-riboviz-env.sh` to configure paths](./install.md#create-set-riboviz-envsh-to-configure-paths), then run:

```console
$ source $HOME/set-riboviz-env.sh
```

* Otherwise, run the following (your directory names may be different, depending on the versions of Hisat2 and Bowtie you have):

```console
$ export PATH=~/hisat2-2.1.0:$PATH
$ export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
```

---

## Configure number of processes (optional)

The workflow can instruct the tools it invokes to use a specific number of processes. By default this is 1.

To configure the workflow to use additional processes, in `vignette/vignette_config.yaml` change:

```yaml
num_processes: 1
```

* to the desired number of processes:

```yaml
num_processes: 4
```

---

## Validate configuration

The Nextflow workflow supports a `--validate_only` command-line parameter which allows for the workflow configuration to be validated without running the workflow.

**Tip:** we strongly recommend validating the configuration before doing a live run on data you have not processed before.

Validate configuration:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/vignette_config.yaml --validate_only
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [compassionate_hoover] - revision: caa1a2bf21
Validating configuration only
No such sample file (NotHere): example_missing_file.fastq.gz
Validated configuration
```

### Troubleshooting: `No such file (NotHere): example_missing_file.fastq.gz`

This warning is expected and can be ignored. The vignette configuration file intentionally includes a reference to a missing input file. This demonstrates that if processing of one sample fails it does not block the processing of the other samples.

---

## Run the workflow

To run the Nextflow workflow:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [goofy_jones] - revision: 0e94b1fa62
No such file (NotHere): example_missing_file.fastq.gz
[93/e96307] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[f1/38397f] Submitted process > cutAdapters (WT3AT)
[1a/a01841] Submitted process > buildIndicesrRNA (yeast_rRNA)
[f0/5779d2] Submitted process > cutAdapters (WTnone)
[c2/497da3] Submitted process > hisat2rRNA (WTnone)
[14/479d3b] Submitted process > hisat2rRNA (WT3AT)
[cd/2e10ec] Submitted process > hisat2ORF (WTnone)
[02/1fef3f] Submitted process > trim5pMismatches (WTnone)
[9e/90df9d] Submitted process > samViewSort (WTnone)
[1e/1a3f4d] Submitted process > outputBams (WTnone)
[d8/3d85de] Submitted process > makeBedgraphs (WTnone)
[26/055ef0] Submitted process > bamToH5 (WTnone)
[d6/355a50] Submitted process > hisat2ORF (WT3AT)
[90/78ca31] Submitted process > trim5pMismatches (WT3AT)
[03/dfb679] Submitted process > samViewSort (WT3AT)
[85/ccbb35] Submitted process > outputBams (WT3AT)
[75/7e865e] Submitted process > bamToH5 (WT3AT)
[32/a55251] Submitted process > makeBedgraphs (WT3AT)
[d5/7a7e8b] Submitted process > generateStatsFigs (WTnone)
[f4/5188e3] Submitted process > generateStatsFigs (WT3AT)
Finished processing sample: WTnone
[c1/184176] Submitted process > renameTpms (WTnone)
Finished processing sample: WT3AT
[6f/716c9e] Submitted process > renameTpms (WT3AT)
[c6/bd5296] Submitted process > collateTpms (WTnone, WT3AT)
[a6/3cb086] Submitted process > countReads
Workflow finished! (OK)
```

* For full information on how to run the Nextflow workflow, the options available, and troubleshooting advice, see [Running the RiboViz Nextflow workflow](./prep-riboviz-run-nextflow.md)

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
  orf_map_clean.bam
  orf_map_clean.bam.bai
  orf_map_clean.sam
  orf_map.sam
  rRNA_map.sam
  trim_5p_mismatch.tsv
  trim.fq
  unaligned.fq
WTnone/
  nonrRNA.fq
  orf_map_clean.bam
  orf_map_clean.bam.bai
  orf_map_clean.sam
  orf_map.sam
  rRNA_map.sam
  trim_5p_mismatch.tsv
  trim.fq
  unaligned.fq
```

Outputs in `vignette/output`:

```
WT3AT/
  3ntframe_bygene.tsv
  3ntframe_propbygene.pdf
  metagene_start_stop_read_counts.pdf
  metagene_start_stop_read_counts.tsv
  WT3AT.bam
  WT3AT.bam.bai
  codon_ribodens.pdf
  codon_ribodens.tsv
  features.pdf
  WT3AT.h5
  minus.bedgraph
  plus.bedgraph
  pos_sp_nt_freq.tsv
  metagene_normalized_profile_start_stop.pdf
  metagene_normalized_profile_start_stop.tsv
  read_lengths.pdf
  read_lengths.tsv
  startcodon_ribogridbar.pdf
  startcodon_ribogrid.pdf
  tpms.tsv
WTnone/
  3ntframe_bygene.tsv
  3ntframe_propbygene.pdf
  metagene_start_stop_read_counts.pdf
  metagene_start_stop_read_counts.tsv
  WTnone.bam
  WTnone.bam.bai
  codon_ribodens.pdf
  codon_ribodens.tsv
  features.pdf
  WTnone.h5
  minus.bedgraph
  plus.bedgraph
  pos_sp_nt_freq.tsv
  metagene_normalized_profile_start_stop.pdf
  metagene_normalized_profile_start_stop.tsv
  read_lengths.pdf
  read_lengths.tsv
  startcodon_ribogridbar.pdf
  startcodon_ribogrid.pdf
  tpms.tsv
read_counts_per_file_tsv
TPMs_all_CDS_all_samples.tsv
```

For full information on what the workflow does and the files it outputs, see [What the RiboViz workflow does](./prep-riboviz-operation.md).

---

## Cleaning up to run again

Before rerunning the vignette, delete the auto-generated index, temporary, logs and output directories:

```console
$ rm -rf vignette/index
$ rm -rf vignette/tmp
$ rm -rf vignette/output
```

Delete the working directories:

```console
$ rm -rf work
```

You might also want to do this if you have run the vignette with an error (for example, a missing R package), and then want to run it again from scratch.

---

## Customising the vignette

We suggest copying `vignette/vignette_config.yaml` and the rest of the `vignette/` directory, and then customising it for use with your own datasets.
