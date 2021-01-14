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
  
  # GFF - bam_to_h5.R
  print("========== GFF - bam_to_h5.R ==========")
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

  # GFF - read_count_functions.R
  print("========== GFF - read_count_functions.R ==========")

  gff_df <- readGFFAsDf(gff_file)
  print(gff_df)
  # TODO how to extract content from gff_df?

  gff_names_2 <- unique(gff_df$seqnames)
  print(gff_names_2)
  print(length(gff_names_2)) # 68


  # BAM

#  bam_what <- c("strand", "pos", "qwidth")
#  bam_param <- ScanBamParam(which = gene_location, what = bam_what)
#  bam_data <- scanBam(bam_file, param=bam_param)
#  bam_data <- scanBam(bam_file)
#  print(bam_data)

  # H5

  print("========== H5 - read_count_functions.R ==========")

  from_h5 <- rhdf5::h5ls(h5_file, recursive = 1)
  print(from_h5)
# TODO how to extract content from from_h5?

  h5_names <- from_h5$name
  print(h5_names)
  print(length(h5_names)) # 68

  expect_equal(length(gff_names_2), length(h5_names), info = "Mismatch in sequence names")



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
