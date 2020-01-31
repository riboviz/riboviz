# Content and provenance of repository data files

This page documents the content and provenance of the data files within the repository.

---

## *Saccharomyces cerevisiae* (yeast) genome and annotation data

Transcript sequences:

```
data/yeast_CDS_w_250utrs.fa
```

Coding sequence locations:

```
data/yeast_CDS_w_250utrs.gff3
```

These files hold S288c annotations and ORF sequences.

These files were created as follows:

* The file [genome release R64-2-1](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz) (file name `S288C_reference_genome_R64-2-1_20150113.tgz`) was downloaded from the [Saccharomyces Genome Database](https://www.yeastgenome.org/).
* The files `saccharomyces_cerevisiae_R64-2-1_20150113.gff` and `S288C_reference_sequence_R64-2-1_20150113.fsa` were extracted from the `.tgz` file.
* The sequence and annotation files for the whole approximate *Saccharomyces cerevisiae* transcriptome were prepared using [script_for_transcript_annotation.Rmd](../rmarkdown/script_for_transcript_annotation.Rmd).

The files can be used as inputs to RiboViz. However, `yeast_CDS_w_250utrs.fa` and `yeast_CDS_w_250utrs.gff3` were downsampled to provide a manageable data set for demonstration purposes, as described in the next section.

---

## Downsampled *Saccharomyces cerevisiae* (yeast) genome and annotation data

Transcript sequences to align to, from just the left arm of chromosome 1:

```
vignette/input/yeast_YAL_CDS_w_250utrs.fa
```

Matched genome feature file, specifying coding sequences locations (start and stop coordinates):

```
vignette/input/yeast_YAL_CDS_w_250utrs.gff3
```

As the yeast data files described in the previous section are very large, these files were downsampled for demonstration processes. The data files `yeast_CDS_w_250utrs.fa` and `yeast_CDS_w_250utrs.gff3` were processed by filtering only ORFs in the left arm of chromosome 1, for which the ORF names start with `YALnnnx`. This produced the above yeast genome and annotation data files.

The document [Appendix A1: Yeast Nomenclature Systematic Open Reading Frame (ORF) and Other Genetic Designations](https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527636778.app1) describes the ORF naming convention.

The files can be used as inputs to RiboViz.

---

## Ribosomal RNA (rRNA) contaminants to remove

rRNA sequences to avoid aligning to:

```
vignette/input/yeast_rRNA_R64-1-1.fa
```

This files was created as follows:

* The file [genome release R64-1-1](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz) (file name `S288C_reference_genome_R64-1-1_20110203.tgz`) was downloaded from the [Saccharomyces Genome Database](https://www.yeastgenome.org/).
* The file `rna_coding_R64-1-1_20110203.fasta` was extracted from the `.tgz` file.
* Selected `RDN-n-n` sequences were copied and pasted from this file.

The file can be used as an input to RiboViz.

---

## Additional yeast-specific data

Position of codons within each gene (the numbering ignores the first 200 codons):

```
data/yeast_codon_pos_i200.RData
```

This file was produced using [script_for_transcript_annotation.Rmd](../rmarkdown/script_for_transcript_annotation.Rmd) as part of the preparation described in [Saccharomyces cerevisiae (yeast) genome and annotation data](#saccharomyces-cerevisiae-yeast-genome-and-annotation-data) above.

Features to correlate with ORFs:

```
data/yeast_features.tsv
```

Data within this file was derived as follows:

1. `Length_log10`: [genome release R64-2-1](https://downloads.yeastgenome.org/sequence/S288C_reference/) (GFF) from the [Saccharomyces Genome Database](https://www.yeastgenome.org/).
2. `utr`, 5' UTR length: Arribere, J.A. and Wendy V. Gilbert, W.V. "Roles for transcript leaders in translation and mRNA decay revealed by transcript leader sequencing", Genome Res. 2013. 23: 977-987 doi:[10.1101/gr.150342.112](http://doi.org/10.1101/gr.150342.112)
3. `polyA` length: Subtelny, A.O. et al. "Poly(A)-tail profiling reveals an embryonic switch in translational control", Nature, 508(66), 29/01/2019 doi:[10.1038/nature13007](http://doi.org/10.1038/nature13007)
4. `uATGs`: Estimated from 2 by counting the upstream ATGs in the annotated 5'UTR
5. `utr_gc`: Estimated from 2 by calculating proportion of G/C in the annotated 5' UTR.
6. `FE_cap`: Estimated from 2. using sequences of length 70 nts from the 5' end of the mRNA transcript with folding energies calculated at 37 degress Centigrade following [Supplementary Methods](https://www.cell.com/cms/10.1016/j.celrep.2016.01.043/attachment/257faf34-ff8f-4071-a642-bfdb531c75b8/mmc1) for Weinberg et al. 2016 "Improved Ribosome-Footprint and mRNA Measurements Provide Insights into Dynamics and Regulation of Yeast Translation", Cell Reports, 14(7), 23 February 2016, 1787-1799 doi: [10.1016/j.celrep.2016.01.043](https://doi.org/10.1016/j.celrep.2016.01.043). Calculations were done using [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html) in the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) package.
7. `FE_atg`: Estimated from 30 nt upstream from ATG.

tRNA estimates:

```
data/yeast_tRNAs.tsv
```

A-site displacement values (based on standard yeast data from Ingolia 2009, Weinberg & Shah 2016):

```
data/yeast_standard_asite_disp_length.txt
```

These files are all read by [generate_stats_figs.R](../rscripts/generate_stats_figs.R) to help with generating plots and tables of results data.

---

## Downsampled ribosome profiling data from *Saccharomyces cerevisiae*

~1mi-sampled RPFs wild-type no additive:

```
vignette/input/SRR1042855_s1mi.fastq.gz
```

~1mi-sampled RPFs wild-type + 3-AT:

```
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

The data was sampled uniformly at random 1/50 reads from each file, producing about 1 million reads total, to produce the downsampled read data files.

---

## Simulated FASTQ test files

`data/simdata/` folder:

```
multiplex_barcodes.tsv
multiplex.fastq
multiplex_umi_barcode_adaptor.fastq
multiplex_umi_barcode.fastq
umi3.fastq
umi3_umi_adaptor.fastq
umi3_umi.fastq
umi5_umi3.fastq
umi5_umi3_umi_adaptor.fastq
umi5_umi3_umi.fastq
deplex/num_reads.tsv
deplex/Tag0.fastq
deplex/Tag1.fastq
deplex/Tag2.fastq
deplex/Unassigned.fastq
```

These files are simple simulated FASTQ files to test adaptor trimming, UMI extraction and deduplication using UMI-tools when invoked from within the RiboViz workflow.

These files were created by running `riboviz/tools/create_fastq_simdata.py`.

See [Create simulated FASTQ files](developer-guide.md#create-simulated-fastq-files).

The files can be used as inputs to RiboViz.

---

## Demultiplexing test files

```
data/demultiplex/
```

Test data for `riboviz/tools/demultiplex_fastq.py`.

Data was imported from https://github.com/ewallace/pyRNATagSeq, commit 6ffd465fb0a80d2134bad9d2147c877c3b363720 (Thu May 11 23:44:13 2017).

* `Sample_4reads_R1.fastq.gz`: artificial sample with 4 read 1s.
* `Sample_4reads_R2.fastq.gz`: 4 read 2s corresponding to `Sample_4reads_R1.fastq.gz`.
* `Sample_init10000_R1.fastq.gz`: Initial 10000 read 1s from a paired-end S. cerevisiae dataset.
* `Sample_init10000_R2.fastq.gz`: 10000 read 2s corresponding to `Sample_init10000_R1.fastq.gz`.
* `TagSeqBarcodedOligos2015.txt`: TagSeq barcoded oligos used in Shishkin, et al. (2015). "Simultaneous generation of many RNA-seq libraries in a single reaction", Nature Methods, 12(4), 323–325. doi: [10.1038/nmeth.3313](http://doi.org/10.1038/nmeth.3313)

---

## `riboviz/test/` test data files

### `riboviz.test.test_trim_5p_mismatch` test data files

```
riboviz/test/data/trim_5p_mismatch.sam
riboviz/test/data/trim_5pos5neg.sam
```

These files are used by `riboviz.test.test_trim_5p_mismatch` for testing `riboviz.trim_5p_mismatch`.

These files were created by running RiboViz using `vignette/vignette_config.yaml` and the data in `vignette/input/`. Lines were copied and pasted from the SAM files output then these lines were manually edited to produce a desired range of outcomes.

### `riboviz.test.test_sam_bam` test data files

```
WTnone_rRNA_map_20.sam
WTnone_rRNA_map_20.bam
WTnone_rRNA_map_20.bam.bai
WTnone_rRNA_map_6_primary.sam
WTnone_rRNA_map_6_primary.bam
WTnone_rRNA_map_6_primary.bam.bai
WTnone_rRNA_map_14_secondary.sam
WTnone_rRNA_map_14_secondary.bam
WTnone_rRNA_map_14_secondary.bam.bai
```

The SAM files were created from the file `tmp/WTnone/rRNA_map.sam` from a run of the vignette (using RiboViz version commit 9efaf93, 08/10/2020):

* `WTnone_rRNA_map_20.sam`: the first 20 sequences from `rRNA_map.sam`.
* `WTnone_rRNA_map_6_primary.sam`: the 6 mapped (primary) sequences from `WTnone_rRNA_map_20.sam`.
* `WTnone_rRNA_map_14_secondary.sam`: the 14 remaining unmapped, non-primary, sequences from `WTnone_rRNA_map_20.sam`.

The BAM and BAI files were created as follows:

```console
$ samtools view -b <FILE>.sam  | samtools sort -@1 -O bam -o <FILE>.bam
$ samtools index <FILE>.bam 
```
