rm(list=ls())

library(rhdf5)
library(ggplot2)
library(dplyr)

# load riboviz functions
riboviz_rscripts_dirname <- "."
source(file.path(riboviz_rscripts_dirname, "stats_figs_block_functions.R"))
source(file.path(riboviz_rscripts_dirname, "read_count_functions.R"))

# data from YAML read with optparse in generate_stats_figs.R
# asite_disp_length_file <- "../data/yeast_standard_asite_disp_length.txt"
hd_file <- "../green/output/sample_1/sample_1.h5"
orf_gff_file <- "../green/input/scer.transcripts.20cds20.gff3"
dataset <- "green"
min_read_length <- 10
max_read_length <- 50
nnt_buffer <- 20
nnt_gene <- 0

# process files
# asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(asite_disp_length_file)
gff_df <- readGFFAsDf(orf_gff_file)
gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
read_range <- min_read_length:max_read_length

# determine range of read lengths
read_lengths_counts <- CalculateReadLengths(gene_names, dataset, hd_file)
min_obs_length <- sapply(seq(nrow(read_lengths_counts)),
                         function(x) {
                           all(read_lengths_counts$Counts[1:x]==0)
                         })
min_obs_length <- read_lengths_counts$Length[which.min(min_obs_length)]
max_obs_length <- sapply(seq(nrow(read_lengths_counts)),
                         function(x) {
                           all(read_lengths_counts$Counts[(x+1):nrow(read_lengths_counts)]==0)
                         })
max_obs_length <- read_lengths_counts$Length[which.max(max_obs_length)]
if(max_obs_length < min_obs_length) { max_obs_length <- max_read_length }

### start codon
# subset to reads around start codon
start_codon <- AllGenes5StartPositionLengthCountsTibble(gene_names, dataset, hd_file, gff_df)
# determine range of valid distances between 5' end and A site
start_codon_max_dist5 <- sapply(min_obs_length:max_obs_length,
                                function(x) {
                                  rpf_length_cts <- subset(start_codon, ReadLen==x)
                                  rpf_length_cts$dist5 <- -rpf_length_cts$Pos + 3
                                  max_dist5 <- sapply(seq(nrow(rpf_length_cts)),
                                                      function(x) {
                                                        any(rpf_length_cts$Counts[1:x]!=0)
                                                      })
                                  max_dist5 <- rpf_length_cts$dist5[which.max(max_dist5)]
                                  return(max_dist5)
                                })


start_codon_plot <- ggplot(start_codon, aes(x=-Pos, y=ReadLen, fill=Counts)) +
  geom_tile(linetype=1, col=1) + theme_classic() +
  scale_fill_gradient(low="white", high="darkblue") +
  # scale_x_continuous(breaks=seq(min(start_codon$start_dist), max(start_codon$start_dist))) +
  # scale_y_continuous(breaks=fp_min_length:fp_max_length) +
  xlab("position relative to start codon") + ylab("footprint length") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  ggtitle("Riboviz vignette data",
          subtitle="start codon footprint density") +
  scale_x_reverse()
