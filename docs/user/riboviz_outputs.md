# riboviz output files and figures

This document lists the output files from a typical riboviz run, along with short descriptions.

After a riboviz run, many output files are produced within the output directory.
The output directory is specified by the parameter `dir_out` in `config.yaml`.

# Output files for an entire run

* `TPMs_collated.tsv`: file with the transcripts per million (tpm) for all successfully processed samples.
* `read_counts.tsv`: a [read counts file](#read-counts-file) (only if `count_reads: TRUE`).

# Output files for each sample

For each sample (`<SAMPLE_ID>`), intermediate files are produced in a sample-specific subdirectory (`<SAMPLE_ID>`):

* `<SAMPLE_ID>.bam`: BAM file of reads mapped to transcripts, which can be directly used in genome browsers.
* `<SAMPLE_ID>.bam.bai`: BAM index file for `<SAMPLE_ID>.bam`.
* `minus.bedgraph`: bedgraph of reads from minus strand (if `make_bedgraph: TRUE`).
* `plus.bedgraph`: bedgraph of reads from plus strand (if `make_bedgraph: TRUE`).
* `<SAMPLE_ID>.h5`: length-sensitive alignments in compressed h5 format.
* `3nt_periodicity.tsv`
* `3nt_periodicity.pdf`
* `read_lengths.tsv`
* `read_lengths.pdf`
* `pos_sp_nt_freq.tsv`
* `pos_sp_rpf_norm_reads.pdf`
* `pos_sp_rpf_norm_reads.tsv`
* `features.pdf`: only output if `--features-file` was defined.
* `tpms.tsv`
* `codon_ribodens.tsv`: only output if `--t-rna-file` and `--codon-positions-file` were defined.
* `codon_ribodens.pdf`: only output if `--t-rna-file` and `--codon-positions-file` were defined.
* `startcodon_ribogridbar.pdf`
* `startcodon_ribogrid.pdf`
* `3ntframe_bygene.tsv`: only output if `--asite-disp-length-file` was defined.
* `3ntframe_propbygene.pdf`: only output if `--asite-disp-length-file` was defined.
* `<SAMPLE_ID>_output_report.html`: only output if `run_static_html: TRUE`.
