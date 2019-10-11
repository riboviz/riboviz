# Content and provenance of repository data files

This page documents the content and provenance of the data files within the repository.

## *Saccharomyces cerevisiae* (yeast) genome and annotation data

```
data/yeast_CDS_w_250utrs.fa
data/yeast_CDS_w_250utrs.gff3
data/yeast_codon_pos_i200.RData
```

These files hold S288c annotations and ORF sequences:

* `data/yeast_CDS_w_250utrs.fa`: transcript sequences.
* `data/yeast_CDS_w_250utrs.gff3`: coding sequences locations.
* `data/yeast_codon_pos_i200.RData`: position of codons within each gene (the numbering ignores the first 200 codons).

These files were created as follows:

* The file [genome release R64-2-1](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz) (file name `S288C_reference_genome_R64-2-1_20150113.tgz`) was downloaded from the [Saccharomyces Genome Database](https://www.yeastgenome.org/).
* The files `saccharomyces_cerevisiae_R64-2-1_20150113.gff` and `S288C_reference_sequence_R64-2-1_20150113.fsa` were extracted from the `.tgz` file.
* The sequence and annotation files for the whole approximate *Saccharomyces cerevisiae* transcriptome were prepared using [script_for_transcript_annotation.Rmd](../rmarkdown/script_for_transcript_annotation.Rmd). 

`yeast_codon_pos_i200.RData` is read by [generate_stats_figs.R](../rscripts/generate_stats_figs.R) to help with generating plots and tables of results data.

The files can be used as inputs to RiboViz. However, `yeast_CDS_w_250utrs.fa` and `yeast_CDS_w_250utrs.gff3` were downsampled to provide a manageable data set for demonstration purposes, as described in the next section.

---

## Downsampled *Saccharomyces cerevisiae* (yeast) genome and annotation data

```
vignette/input/yeast_YAL_CDS_w_250utrs.fa
vignette/input/yeast_YAL_CDS_w_250utrs.gff3
```

As the yeast data files described in the previous section are very large, these files were downsampled for demonstration processes. The data files `yeast_CDS_w_250utrs.fa` and `yeast_CDS_w_250utrs.gff3` were processed by filtering only ORFs in the left arm of chromosome 1, for which the ORF names start with `YALnnnx`. This produced the yeast genome and annotation data files:

* `vignette/input/yeast_YAL_CDS_w_250utrs.fa`: transcript sequences to align to, from just the left arm of chromosome 1.
* `vignette/input/yeast_YAL_CDS_w_250utrs.gff3`: matched genome feature file, specifying coding sequences locations (start and stop coordinates).

The document [Appendix A1: Yeast Nomenclature Systematic Open Reading Frame (ORF) and Other Genetic Designations](https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527636778.app1) describes the ORF naming convention.

The files can be used as inputs to RiboViz.

---

## Ribosomal RNA (rRNA) contaminants to remove

```
vignette/input/yeast_rRNA_R64-1-1.fa
```

This file holds rRNA sequences to avoid aligning to.

This files was created as follows:

* The file [genome release R64-1-1](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz) (file name `S288C_reference_genome_R64-1-1_20110203.tgz`) was downloaded from the [Saccharomyces Genome Database](https://www.yeastgenome.org/).
* The file `rna_coding_R64-1-1_20110203.fasta` was extracted from the `.tgz` file.
* Selected `RDN-n-n` sequences were copied and pasted from this file.

The file can be used as an input to RiboViz.

---

## Downsampled ribosome profiling data from *Saccharomyces cerevisiae*

```
vignette/input/SRR1042855_s1mi.fastq.gz
vignette/input/SRR1042864_s1mi.fastq.gz
```

These read data files are downsampled ribosome profiling data from *Saccharomyces cerevisiae*. The data has been downsampled to provide a dataset that is realistic, but small enough to run quickly.

The data is from the paper Guydosh N.R. and Green R. "[Dom34 rescues ribosomes in 3' untranslated regions](https://www.ncbi.nlm.nih.gov/pubmed/24581494)", Cell. 2014 Feb 27;156(5):950-62. doi: [10.1016/j.cell.2014.02.006](https://doi.org/10.1016/j.cell.2014.02.006). The NCBI accession for the whole dataset is #GSE52968:

* SRX386986: GSM1279570: wild-type no additive, [SRR1042855](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1042855)
* SRX386995: GSM1279579: wild-type plus 3-AT, [SRR1042864](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1042864)

In July 2017, these files were imported using NCBI's [fastq-dump](https://ncbi.github.io/sra-tools/fastq-dump.html) and gzipped to produce:

```
SRR1042855.fastq.gz
SRR1042864.fastq.gz
```

(these files are not in the repository)

These files can alternatively be accessed via [SRA Explorer](https://ewels.github.io/sra-explorer/#):

* Search for: SRR1042855
* Select GSM1279570: wild-type no additive; Saccharomyces cerevisiae; OTHER
* Click Add 1 to collection
* Search for: SRR1042864
* Select GSM1279579: wild-type 3-AT; Saccharomyces cerevisiae; OTHER	
* Click Add 1 to collection
* Click 2 saved datasets
* Click Bash script for downloading FastQ files
* Click Download
* Run the bash script, e.g.

```console
$ source sra_explorer_fastq_download.sh
```

* **Warning:** the total download time may take ~40 minutes or more. The files are  1.5GB and 2.2GB respectively.

The data was sampled uniformly at random 1/50 reads from each file, producing about 1 million reads total, to produce the downsampled read data files:

* `vignette/input/SRR1042855_s1mi.fastq.gz`: ~1mi-sampled RPFs wild-type no additive.
* `vignette/input/SRR1042864_s1mi.fastq.gz`: ~1mi-sampled RPFs wild-type + 3-AT.

---

## Simulated FASTQ test files

```
data/simdata_UMI5and3_4nt_adaptor.fastq
data/simdata_UMI5and3_4nt.fastq
data/simdata_extracted_UMI5and3_4nt.fastq
data/simdata_extracted_UMI3_4nt.fastq
data/simdata_UMI3_4nt.fastq
```

These files are simple simulated FASTQ files to test adaptor trimming, UMI extraction and deduplication using UMI-tools when invoked from within the RiboViz workflow.

These files were created by running `riboviz/tools/create_fastq_examples.py`.

See [Create simulated FASTQ files](developer-guide.md#create-simulated-fastq-files).

The files can be used as inputs to RiboViz.

---

## `trim_5p_mismatch.py` test files

```
data/testdata_trim_5p_mismatch.sam
data/testdata_trim_5pos5neg.sam
```

These files are used for testing [trim_5p_mismatch.py](../riboviz/tools/trim_5p_mismatch.py).

These files were created by running RiboViz using `vignette/vignette_config.yaml` and the data in `vignette/input`. Lines were copied and pasted from the SAM files output then these lines were manually edited to produce a desired range of outcomes.
