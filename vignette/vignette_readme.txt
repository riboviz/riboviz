vignette_readme.txt

Folder vignette contains the files to run a vignette of the Riboviz package up to the h5 output.

To run the vignette, from the Riboviz root directory:
bash scripts/prepRiboviz.sh vignette/vignette_config.yaml


Contents:

Config etc in root of vignette:
  - vignette_readme.txt (this file)
  - vignette_config.yaml, contains all config information
  
Input data, in vignette/input:
  - yeast_YAL_CDS_w_250utrs.fa sequence file to align to from just left arm of chromosome 1
  - yeast_YAL_CDS_w_250utrs.gff3 matched genome feature file to specify start and stop co-ordinates
  - yeast_rRNA_R64-1-1.fa rRNA sequence file to avoid aligning to
  - SRR1042855_s1mi.fastq.gz about 1mi sampled RPFs wild-type no additive from Guydosh & Green 2014
  - SRR1042864_s1mi.fastq.gz about 1mi sampled RPFs wild-type + 3-AT from Guydosh & Green 2014

Running the vignette produces the following outputs in folders vignette/index, vignette/tmp, and vignette/output.
To rerun the vignette, first delete these folders.

Index files, in vignette/index:
  - YAL_CDS_w_250.*.ht2, hisat2 index from yeast_YAL_CDS_w_250utrs
  - yeast_rRNA.*.ht2, hisat2 index from yeast_rRNA_R64-1-1

Intermediate outputs, in vignette/tmp:
  - *_trim.fq, trimmed reads
  - *_nonrRNA.fq, trimmed non-rRNA reads
  - *_rRNA_map.sam, rRNA-mapped reads
  - *_orf_map.sam, orf-mapped reads
  - *_orf_map_clean.sam, orf-mapped reads with mismatched nt trimmed
  - *_unaligned.sam, unaligned reads
(note these are uncompressed and large)

Outputs, in vignette/output
  - WTnone.bam
  - WTnone.bam.bai
  - WTnone_minus.bedgraph
  - WTnone_plus.bedgraph
  - WTnone.h5
  - WTnone_3nt_periodicity.tsv
  - WTnone_3nt_periodicity.pdf
  - WTnone_read_lengths.tsv
  - WTnone_read_lengths.pdf
  - WTnone_pos_sp_nt_freq.tsv
  - WTnone_pos_sp_rpf_norm_reads.pdf
  - WTnone_pos_sp_rpf_norm_reads.tsv
  - WTnone_features.pdf
  - WTnone_features.tsv
  - WTnone_codon_ribodens.tsv
  - WTnone_codon_ribodens.pdf

  - WT3AT.bam
  - WT3AT.bam.bai
  - WT3AT_minus.bedgraph
  - WT3AT_plus.bedgraph
  - WT3AT.h5
  - WT3AT_3nt_periodicity.tsv
  - WT3AT_3nt_periodicity.pdf
  - WT3AT_read_lengths.tsv
  - WT3AT_read_lengths.pdf
  - WT3AT_pos_sp_nt_freq.tsv
  - WT3AT_pos_sp_rpf_norm_reads.pdf
  - WT3AT_pos_sp_rpf_norm_reads.tsv
  - WT3AT_features.pdf
  - WT3AT_features.tsv
  - WT3AT_codon_ribodens.tsv
  - WT3AT_codon_ribodens.pdf

