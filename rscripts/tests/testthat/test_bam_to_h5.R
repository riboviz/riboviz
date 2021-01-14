# testthat tests for bam_to_h5.R
#
# This assumes the following files are in the path:
#
# rscripts/bam_to_h5.R
# vignette/input/yeast_YAL_CDS_w_250utrs.gff3
# vignette/output/WTnone/WTnone.bam
#
# The test runs bam_to_h5.R using the GFF and BAM file then validates
# the .h5 file created based upon its expected qualities.
#
# Note: for rapid development this does not currently run bam_to_h5.R
# but validates vignette/output/WTnone/WTnone.h5. TODO remove this
# comment when development is complete.
#
# To run interactively:
#
# test_file("rscripts/tests/testthat/test_bam_to_h5.R")
#
# To run from console:
#
# Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_bam_to_h5.R")'

library(glue)
library(here)
library(testthat)
library(withr)

source(here::here("rscripts", "read_count_functions.R"))

print(here())
bam_to_h5 <- here::here("rscripts/bam_to_h5.R")
print(bam_to_h5)
gff_file = here::here("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
print(gff_file)
bam_file = here::here("vignette/output/WTnone/WTnone.bam")
print(bam_file)
h5_file = here::here("test.h5")
h5_file = here::here("vignette/output/WTnone/WTnone.h5") # TODO remove
print(h5_file)

context("test_bam_to_h5.R")

delete_file <- function(file_name) {
  print("Invoking delete_file fixture...")
  if (file.exists(file_name)) {
    # file.remove(file_name) # TODO uncomment
  }
}

# TODO fix read_count_functions.R
# GetGeneReadLength <- function(gene, hd_file)
# =>
# GetGeneReadLength <- function(gene, dataset, hd_file)

# TODO add to read_count_functions.R
GetBufferLeft <- function(gene, dataset, hd_file){
  rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["buffer_left"]]
}
GetBufferRight <- function(gene, dataset, hd_file){
  rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["buffer_right"]]
}
GetStartCodonPos <- function(gene, dataset, hd_file){
  rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["start_codon_pos"]]
}
GetStopCodonPos <- function(gene, dataset, hd_file){
  rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["stop_codon_pos"]]
}
GetMappedReadLengths <- function(gene, dataset, hd_file){
  rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["lengths"]]
}

test_that("Run bam_to_h5.R and validate H5 file", {
  withr::defer(delete_file(h5_file))

  expect_equal(0, 0, info = "Example assertion") # TODO remove

  min_read_length <- 10
  max_read_length <- 50
  buffer <- 250
  primary_id <- "Name"
  secondary_id <- "NULL"
  dataset <- "vignette"
  is_riboviz_gff <- TRUE
  stop_in_cds <- FALSE
  
  if (FALSE) # TODO uncomment
  {
  bam_to_h5_cmd_template <- "Rscript --vanilla {bam_to_h5} --num-processes=1 --min-read-length={min_read_length} --max-read-length={max_read_length} --buffer={buffer} --primary-id={primary_id} --secondary-id={secondary_id} --dataset={dataset} --bam-file={bam_file} --hd-file={h5_file} --orf-gff-file={gff_file} --is-riboviz-gff={is_riboviz_gff} --stop-in-cds={stop_in_cds}"
  print(bam_to_h5_cmd_template)
  cmd <- glue(bam_to_h5_cmd_template)
  print(cmd)
  exit_code <- system(cmd)
  print(glue("bam_to_h5.R exit code: {exit_code}"))
  expect_equal(exit_code, 0, info = "Unexpected exit code from bam_to_h5.R")
  }

  ##### EXTRACT GFF #####

  print("========== GFF - bam_to_h5.R-style ==========")
  gff <- rtracklayer::readGFFAsGRanges(gff_file)
  print(gff)
  gff_names <- mcols(gff)
  print(gff_names)
  gff_names <- gff_names["Name"] # --primary_id  
  print(gff_names)
  gff_names <- gff_names[,1] 
  print(gff_names)
  gff_names <- unique(gff_names)
  print(gff_names)
  print(length(gff_names)) # 68
  gff_pid <- mcols(gff)["Name"][,1]
  print(gff_pid)

  gene <- "YAL062W"
  gene_location <- gff[gff_pid == gene]
  print(gene_location)
  gene_location <- gff["Name" == gene]
  print(gene_location)

  print("========== GFF - generate_stats_figs.R/read_count_functions.R-style ==========")
  # From generate_stats_figs.R calling read_count_functions.R
  gff_df <- readGFFAsDf(gff_file)
  print(gff_df)
  gff_names_2 <- unique(gff_df$seqnames) # Tibble data frame
  print(gff_names_2)
  print(length(gff_names_2)) # 68

  # TODO

  ##### EXTRACT BAM #####

  # TODO
  
  ##### EXTRACT AND VALIDATE H5 #####

  print("========== H5 ==========")

  # From generate_stats_figs.R
  h5_data <- rhdf5::h5ls(h5_file, recursive = 1)
  print(class(h5_data)) # data.frame
  print(typeof(h5_data)) # list
  print(h5_data)
  #  group      name     otype dclass dim
  # 0      /   YAL001C H5I_GROUP           
  # 58     /   YAL062W H5I_GROUP
  h5_names <- h5_data$name
  print(class(h5_names)) # character
  print(typeof(h5_names)) # character
  print(h5_names)
  #  [1] "YAL001C"   "YAL002W"   ...
  #  ...
  # [67] "YAL067W-A" "YAL068C"  
  print(length(h5_names)) # 68

  expect_equal(length(gff_names_2), length(h5_names), info = "Mismatch in sequence names between GFF and H5")

  # From read_count_functions.R
  gene <- "YAL062W"

  print("--- H5 'data' ---")
  h5_data <- GetGeneDatamatrix(gene, dataset, h5_file)
  # print(h5_data) # Verbose
  # TODO print rows which were expected to have 1s
  print(dim(h5_data)) # 41 1874

  print("--- H5 'reads_by_len' ---")
  h5_reads_by_len <- GetGeneReadLength(gene, h5_file)
  print(h5_reads_by_len)
  # [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  # [39] 0 0 0
  print(length(h5_reads_by_len)) # 41

  print("--- H5 length ('stop_codon_pos' - 'start_codon_pos') ---")
  h5_gene_length <- GetGeneLength(gene, dataset, h5_file)
  print(h5_gene_length) # 1371

  print("--- H5 'reads_total' ---")
  h5_reads_total <- GetGeneReadsTotal(gene, dataset, h5_file)
  print(h5_reads_total) # 2

  print("--- H5 density ('reads_total' / length) ---")
  h5_density <- GetGeneReadDensity(gene, dataset, h5_file, buffer = buffer)
  print(h5_density) # 0.00140746

  print("--- H5 'buffer_left' ---")
  h5_buffer_left <- GetBufferLeft(gene, dataset, h5_file)
  print(h5_buffer_left) # 250
  print("--- H5 'buffer_right' ---")
  h5_buffer_right <- GetBufferRight(gene, dataset, h5_file)
  print(h5_buffer_right) # 250
  print("--- H5 'start_codon_pos' ---")
  h5_start_codon_pos <- GetStartCodonPos(gene, dataset, h5_file)
  print(h5_start_codon_pos) # 251 252 253
  print(length(h5_start_codon_pos)) # 3
  print("--- H5 'stop_codon_pos' ---")
  h5_stop_codon_pos <- GetStopCodonPos(gene, dataset, h5_file)
  print(h5_stop_codon_pos) # 1622 1623 1624
  print(length(h5_stop_codon_pos)) # 3
  print("--- H5 'lengths' ---")
  h5_lengths <- GetMappedReadLengths(gene, dataset, h5_file)
  print(h5_lengths)
  # [1] 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
  # [26] 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
  print(length(h5_lengths)) # 41
})

# How H5 template is populated - basis for tests
#
# * 'buffer_left': number of nucleotides upstream of the start codon (ATG) (UTR5 length) (from bam_to_h5.R command-line)
# * 'buffer_right': number of nucleotides downstream of the stop codon (TAA/TAG/TGA) (UTR3 length) (from bam_to_h5.R command-line)
# * 'start_codon_pos': Positions corresponding to start codon of CDS in organism sequence (from GFF)
# * 'stop_codon_pos': Positions corresponding to stop codon of CDS in organism sequence (from GFF)
# * 'lengths' : Lengths of mapped reads.
# * 'reads_by_len': Counts of number of ribosome sequences of each length (from BAM).
# * 'reads_total': Total number of ribosome sequences (from BAM, equal to number of non-zero reads in 'reads_by_len').
# * 'data': Positions and lengths of ribosome sequences within the organism data (from BAM).
# 
# * Sequence names in H5 equal sequence names in GFF
# * Every sequence in H5 equal or are superset of sequence names in BAM
# * Every sequence in BAM has "non-zero" DATA sequence in H5
# * Every sequence not in BAM has "zero" DATA sequence in H5
# 
# * H5 has '<sequence ID from FASTA|GFF|reference sequence name from BAM>' GROUP for each sequence
# * Each sequence GROUP has '<dataset> sub-GROUP # vignette default
# * '<dataset> sub-GROUP has 'reads' sub-GROUP
# * 'reads sub-GROUP has attributes:
#   - 'buffer_left': '<buffer|GFF UTR5 length>' # 250 default
#   - 'buffer_right': '<buffer|GFF UTR5 length>' # 250 default
#   - 'lengths':
#      - Length: <max-read-length>-<min-read-length>+1) / (<max-read-length>-<min-read-length>+1
#      - Values: <min-read-length>, <min-read-length>+1, ... ,<min-read-length>+<m>-1, <min-read-length>+<m>, <min-read-length>+<m>+1, ... ,<min-read-length>+<n>-1, <min-read-length>+<n>, <min-read-length>+<n>+1, ... ,<min-read-length>+<max-read-length>
#   - 'reads_by_len':
#      - Length: '<max-read-length>-<min-read-length>+1' # Default, 50-10+1=41
#      - Values: 'reads_by_len[i]' = sum of sequences in BAM which have length equal to 'lengths[i]'. Also equals sum of 'DATA[*, i]' i.e. sum across all positions for a specific length.
#   - 'reads_total': <number of non-zero values in reads_by_len>
#   - 'start_codon_pos>': <position of 1st nt of CDS start codon from GFF>,...,<position of 3rd nt of CDS start codon from GFF>
#   - 'stop_codon_pos': <position of 1st nt of CDS stop codon from GFF>,...,<position of 3rd nt of CDS stop codon from GFF>
# * 'reads sub-GROUP has DATA 'data':
#    - Dimensions: '<position of final codon of UTR3 from GFF|Length of sequence from BAM header LN value>' x '<max-read-length>-<min-read-length>+1' ) /
#    - Values:
#      - 'DATA[p, i]' = 1 if there is a sequence from BAM at position 'p'+1 which has length equal to 'lengths[i]', else 0.
#      - 0 <= 'p' <= position of final codon of UTR3 from GFF|Length of sequence from BAM header 'LN' value.
