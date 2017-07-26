vignette_readme.txt

Folder vignette contains the files necessary to run a vignette of the Riboviz package.

To run the vignette, from the riboviz root directory:
prepRiboviz vignette/vignette_config.yaml
## Currently this is ambition not reality.

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

Running the vignette produces the following outputs:

Index files, in vignette/index:


Intermediate outputs, in vignette/tmp:
  - vignette_trim.fq, trimmed reads
  - vignette_nonrRNA.fq, trimmed non-rRNA reads
  - vignette_rRNA_map.sam, rRNA-mapped reads
  - vignette_orf_map.sam, orf-mapped reads
  - vignette_orf_map_clean.sam, orf-mapped reads with mismatched nt trimmed
  - vignette_unaligned.sam, unaligned reads

Outputs, in vignette/output
  - vignette.bam
  - vignette.h5