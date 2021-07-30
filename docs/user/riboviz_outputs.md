# riboviz output files and figures

This document lists the output files from a typical riboviz run, along with short descriptions.

After a riboviz run, many output files are produced within the output directory.
The output directory is specified by the parameter `dir_out` in `config.yaml`.

There are a few output files that collect information for an entire run.
There are many output files that are specific to each sample, which are organized into a separate subdirectory for each sample. 


# Output files for an entire run


## `TPMs_collated.tsv` 

File with the transcripts per million (tpm) for all successfully processed samples.

## `read_counts.tsv` 

A [read counts file](#read-counts-file) (only if `count_reads: TRUE`).


# Output files for each sample

For each sample (`<SAMPLE_ID>`), intermediate files are produced in a sample-specific subdirectory (`<SAMPLE_ID>`).


## `<SAMPLE_ID>_output_report.html` 

This output report in .html format contains multiple figures output by riboviz:

* first figure (link to description below?)
* ... (fill this in!)

Only output if `run_static_html: TRUE`.


## `<SAMPLE_ID>.bam` 

BAM file of reads mapped to transcripts, which can be directly used in genome browsers.


## `<SAMPLE_ID>.bam.bai` 

BAM index file for `<SAMPLE_ID>.bam`.


## `minus.bedgraph` 

Bedgraph of reads from minus strand (if `make_bedgraph: TRUE`).

Because riboviz aligns to the transcriptome, which represents single-stranded positive-sense RNA, there should be very few reads counted in `minus.bedgraph`.


## `plus.bedgraph` 

Bedgraph of reads from plus strand (if `make_bedgraph: TRUE`).

Almost all translated reads should be counted in open reading frames within `plus.bedgraph`, again because riboviz aligns to the transcriptome, which represents single-stranded positive-sense RNA.


## `<SAMPLE_ID>.h5` 

Length-sensitive alignments in compressed h5 format.

## `3nt_periodicity.tsv`


## `3nt_periodicity.pdf`


## `read_lengths.tsv`


## `read_lengths.pdf`


## `pos_sp_nt_freq.tsv`


## `pos_sp_rpf_norm_reads.pdf`


## `pos_sp_rpf_norm_reads.tsv`


## `features.pdf` 

Only output if `--features-file` was defined.


## `tpms.tsv`


## `codon_ribodens.tsv` 

Only output if `--t-rna-file` and `--codon-positions-file` were defined.


## `codon_ribodens.pdf` 

Only output if `--t-rna-file` and `--codon-positions-file` were defined.


## `startcodon_ribogridbar.pdf`


## `startcodon_ribogrid.pdf`


## `3ntframe_bygene.tsv` 

Only output if `--asite-disp-length-file` was defined.


## `3ntframe_propbygene.pdf` 

Only output if `--asite-disp-length-file` was defined.
