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
  - vignette_2009_Ingolia_rpf-rich-1.fastq.gz of example reads
  - scer_10genes.fasta sequence file to align to
  - scer_rRNA.fasta rRNA sequence file to avoid aligning to
  - scer_10genes.gff genome feature file to specify start and stop co-ordinates

Running the vignette produces the following outputs:

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