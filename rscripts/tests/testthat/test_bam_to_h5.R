# testthat tests for bam_to_h5.R
#
# This assumes the following files are in the path:
#
# rscripts/bam_to_h5.R
# vignette/input/yeast_YAL_CDS_w_250utrs.gff3
# vignette/output/WTnone/WTnone.bam
#
# The test runs bam_to_h5.R using the GFF and BAM file then validates
# the .h5 file created based upon its expected qualities given those
# of the input GFF and BAM files and the bam_to_h5.R command-line
# parameters.
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
gff_file <- here::here("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
print(gff_file)
bam_file <- here::here("vignette/output/WTnone/WTnone.bam")
print(bam_file)
h5_file <- here::here("test.h5")
h5_file <- here::here("vignette/output/WTnone/WTnone.h5") # TODO remove
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

  num_read_counts <- max_read_length - min_read_length + 1

  bam_to_h5_cmd_template <- "Rscript --vanilla {bam_to_h5} --num-processes=1 --min-read-length={min_read_length} --max-read-length={max_read_length} --buffer={buffer} --primary-id={primary_id} --secondary-id={secondary_id} --dataset={dataset} --bam-file={bam_file} --hd-file={h5_file} --orf-gff-file={gff_file} --is-riboviz-gff={is_riboviz_gff} --stop-in-cds={stop_in_cds}"
  print(bam_to_h5_cmd_template)
  cmd <- glue(bam_to_h5_cmd_template)
  print(cmd)
  if (FALSE) # TODO uncomment
  {
  exit_code <- system(cmd)
  print(glue("bam_to_h5.R exit code: {exit_code}"))
  expect_equal(exit_code, 0, info = "Unexpected exit code from bam_to_h5.R")
  }

  ##### EXTRACT GFF #####

  ## bam_to_h5.R-style
  
  print("========== GFF - bam_to_h5.R-style ==========")
  gff <- rtracklayer::readGFFAsGRanges(gff_file)
  print(gff)
  print(length(gff)) # 204
  print(class(gff)) # GRanges attr(,"package") GenomicRanges
  print(typeof(gff)) # S4
  print(levels(gff)) # NULL
  print(dim(gff)) # NULL
  gff_names <- mcols(gff)
  print(gff_names)
  print(length(gff_names)) # 5 
  print(class(gff_names)) # DFrame attr(,"package") S4Vectors
  print(typeof(gff_names)) # S4
  print(levels(gff_names)) # NULL
  print(dim(gff_names)) # 204 5
  gff_names <- gff_names["Name"] # --primary_id  
  print(gff_names)
  gff_names <- gff_names[,1] 
  print(gff_names)
  gff_names <- unique(gff_names)
  print(gff_names)
  print(length(gff_names)) # 68
  print(class(gff_names)) # character
  print(typeof(gff_names)) # character
  print(levels(gff_names)) # NULL
  print(dim(gff_names)) # NULL
  gff_pid <- mcols(gff)["Name"][,1]
  print(gff_pid)
  print(length(gff_pid)) # 204
  print(class(gff_pid)) # character
  print(typeof(gff_pid)) # character
  print(levels(gff_names)) # NULL
  print(dim(gff_pid)) # NULL

  gene <- "YAL062W"
  gene_location <- gff[gff_pid == gene]
  print(gene_location)
  gene_location <- gff["Name" == gene]
  print(gene_location)

  ## read_count_functions.R-style

  print("========== GFF ==========")
  gff_df <- readGFFAsDf(gff_file)
  print(gff_df)
  gff_names <- unique(gff_df$seqnames) # Tibble data frame
  print(gff_names)
  # [1] YAL068C YAL067W-A YAL067C YAL065C YAL064W-B YAL064C-A YAL064W ...
  # [64] YAL007C YAL005C   YAL003W YAL002W YAL001C  
  # 68 Levels: YAL001C YAL002W YAL003W YAL005C YAL007C YAL008W YAL009W ... YAL068C
  print(length(gff_names)) # 68
  print(class(gff_names)) # factor
  print(typeof(gff_names)) # integer
  print(levels(gff_names))
  # [1] "YAL068C" "YAL067W-A" "YAL067C" "YAL065C" "YAL064W-B" "YAL064C-A" ...
  # [67] "YAL067W-A" "YAL068C"  
  print(dim(gff_names)) # NULL

  # TODO
 
  ##### EXTRACT BAM #####
  
  # TODO

  ##### EXTRACT AND VALIDATE H5 #####

  print("========== H5 ==========")

  h5_data <- rhdf5::h5ls(h5_file, recursive = 1)
  print(h5_data)
  #  group      name     otype dclass dim
  # 0      /   YAL001C H5I_GROUP           
  # 58     /   YAL062W H5I_GROUP
  print(length(h5_data)) # 5
  print(class(h5_data)) # data.frame
  print(typeof(h5_data)) # list
  print(levels(h5_data)) # NULL
  print(dim(h5_data)) # 68 5
  h5_names <- h5_data$name
  print(h5_names)
  print(length(h5_names)) # 68
  print(class(h5_names)) # character
  print(typeof(h5_names)) # character
  print(levels(h5_names)) # NULL
  print(dim(h5_names)) # NULL
  #  [1] "YAL001C"   "YAL002W"   ...
  #  ...
  # [67] "YAL067W-A" "YAL068C"  

  expect_equal(length(h5_names), length(gff_names),
    info = "Mismatch in number of sequence names between GFF and H5")
  expect_equal(as.factor(sort(h5_names)), sort(gff_names),
    info = "Mismatch in sequence names between GFF and H5")

  # TODO check h5_names superset of names in BAM.
  
  #
  # Validate values for YAL062W
  #

  gene <- "YAL062W"

  # 'buffer_left': number of nucleotides upstream of the start codon (ATG) (UTR5 length) (from bam_to_h5.R command-line)
  print("buffer_left:")
  h5_buffer_left <- GetGeneBufferLeft(gene, dataset, h5_file) # double
  print(h5_buffer_left) # 250
  expect_equal(h5_buffer_left, buffer, info = "Unexpected buffer_left")
  # TODO Alternative, check against GFF UTR5 length

  # 'buffer_right': number of nucleotides downstream of the stop codon (TAA/TAG/TGA) (UTR3 length) (from bam_to_h5.R command-line)
  print("buffer_right:")
  h5_buffer_right <- GetGeneBufferRight(gene, dataset, h5_file) # integer
  print(h5_buffer_right) # 250
  expect_equal(h5_buffer_right, buffer, info = "Unexpected buffer_right")
  # TODO Alternative, check against GFF UTR3 length

  # TODO generalise
  # 'start_codon_pos': Positions corresponding to start codon of CDS in organism sequence (from GFF)
  # TODO Get position of 1st nt to 3rd nt of CDS start codon from GFF
  expected_start_codons = array(c(251, 252, 253))
  print("start_codon_pos:")
  h5_start_codon_pos <- GetGeneStartCodonPos(gene, dataset, h5_file) # 1D array of 3 integer
  print(h5_start_codon_pos) # 251 252 253
  expect_equal(length(h5_start_codon_pos), 3,
    info = "Unexpected number of start_codon_pos")
  expect_equal(h5_start_codon_pos, expected_start_codons,
    info = "Unexpected start_codon_pos")

  # TODO generalise
  # 'stop_codon_pos': Positions corresponding to stop codon of CDS in organism sequence (from GFF)
  # TODO Get position of 1st nt to 3rd nt of CDS stop codon from GFF
  expected_stop_codons <- array(c(1622, 1623, 1624))
  print("stop_codon_pos:")
  h5_stop_codon_pos <- GetGeneStopCodonPos(gene, dataset, h5_file) # 1D array of 3 integer
  print(h5_stop_codon_pos) # 1622 1623 1624
  expect_equal(length(h5_stop_codon_pos), 3,
    info = "Unexpected number of stop_codon_pos")
  expect_equal(h5_stop_codon_pos, expected_stop_codons,
    info = "Unexpected stop_codon_pos")

  # 'lengths' : Lengths of mapped reads.
  print("lengths:")
  expected_lengths <- as.array(seq(min_read_length, max_read_length))
  h5_lengths <- GetGeneMappedReadLengths(gene, dataset, h5_file) # 1D array of <max_read_length - min_read_length + 1> integer
  print(h5_lengths)
  # [1] 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
  # [26] 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
  expect_equal(length(h5_lengths), num_read_counts, info = "Unexpected number of lengths")
  expect_equal(h5_lengths, expected_lengths, info = "Unexpected lengths")

  # TODO generalise
  # 'reads_by_len': Counts of number of ribosome sequences of each length (from BAM).
  print("reads_by_len:")
  h5_reads_by_len <- GetGeneReadLength(gene, dataset, h5_file) # 1D array of <max_read_length - min_read_length + 1> double
  print(h5_reads_by_len)
  # [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  # [39] 0 0 0
  expect_equal(length(h5_reads_by_len), num_read_counts,
    info = "Unexpected number of reads_by_len")
  # TODO Deduce positions of non-zero values from BAM (reads_by_len[i] = sum of sequences in BAM which have length equal to lengths[i])

  # TODO generalise
  expect_equal(h5_reads_by_len[18], 1,
    info = "Unexpected value for reads_by_len[17]") # 1-indexed R vs 0-indexed H5
  expect_equal(h5_reads_by_len[19], 1,
    info = "Unexpected value for reads_by_len[18]") # 1-indexed R vs 0-indexed H5

  # 'reads_total': Total number of ribosome sequences (from BAM, equal to number of non-zero reads in 'reads_by_len').
  h5_reads_total <- GetGeneReadsTotal(gene, dataset, h5_file) # 1D array of 1 double
  print("reads_total:")
  print(h5_reads_total) # 2
  h5_reads_len_total <- Reduce("+", h5_reads_by_len)
  print("SUM(reads_by_len):")
  print(h5_reads_len_total)
  expect_equal(h5_reads_total[1], h5_reads_len_total,
    info = "reads_total does not equal sum of totals in reads_by_len ")

  # TODO generalise
  # TODO Cross-check against count of BAM sequences which have Flag = 0
  expected_reads_total = 2
  expect_equal(length(h5_reads_total), 1, info = "Unexpected number of reads_total")
  expect_equal(h5_reads_total[1], expected_reads_total, info = "Unexpected reads_total")

  # 'data': Positions and lengths of ribosome sequences within the organism data (from BAM).
  print("data:")
  h5_data <- GetGeneDatamatrix(gene, dataset, h5_file)
  # print(h5_data) # Verbose
  print(length(h5_data)) # 76384
  print(class(h5_data)) # matrix
  print(typeof(h5_data)) # integer
  print(levels(h5_data)) # NULL
  print(dim(h5_data)) # 41 1874
  print(nrow(h5_data)) # 41
  print(ncol(h5_data)) # 1874
  expect_equal(nrow(h5_data), num_read_counts,
    info = "Unexpected number of data rows")
  expected_num_data_cols = h5_stop_codon_pos[3] + buffer
  print(expected_num_data_cols) # 1874
  expect_equal(ncol(h5_data), expected_num_data_cols,
    info = "Unexpected number of data columns")
  # TODO Alternative, get position of final codon of UTR3 from GFF3
  # TODO Alternative, get number of columns from length of sequence from BAM header LN value
  # TODO Check DATA[p, i] = 1 if there is a sequence from BAM at position p+1 which has length equal to lengths[i], else 0.
  
  # TODO check sequence with "non-zeros" is in BAM.
  # TODO check sequence with "zeros" only is not in BAM.

  # TODO generalise
  expected_row <- integer(num_read_counts)
  expected_row[19] = 1 #  1-indexed R vs 0-indexed H5 (18)
  row = h5_data[,752] # 1-indexed R vs 0-indexed H5
  print(row)
  expect_equal(row, expected_row, info = "Unexpected data[751]")

  expected_row <- integer(num_read_counts)
  expected_row[18] = 1 #  1-indexed R vs 0-indexed H5 (17)
  row = h5_data[,753] # 1-indexed R vs 0-indexed H5
  print(row)
  expect_equal(row, expected_row, info = "Unexpected data[752]")

  # TODO Cross-check reads_by_len[i] = sum of DATA[*, i] i.e. sum across all positions for a specific length.
})
