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

suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(withr, quietly = T))
suppressMessages(library(GenomicAlignments, quietly = T))
suppressMessages(library(Rsamtools, quietly = T))

source(here::here("rscripts", "read_count_functions.R"))
print(paste0("here: ", here()))
bam_to_h5 <- here::here("rscripts/bam_to_h5.R")
print(paste0("bam_to_h5.R: ", bam_to_h5))
gff_file <- here::here("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
print(paste0("GFF: ", gff_file))
bam_file <- here::here("vignette/output/WTnone/WTnone.bam") # TODO remove
bam_file <- here::here("WTnone.bam") # TODO remove
print(paste0("BAM: ", bam_file))
h5_file <- here::here("vignette/output/WTnone/WTnone.h5") # TODO remove
h5_file <- here::here("WTnone.h5") # TODO remove
print(paste0("HDF5: ", h5_file))

context("test_bam_to_h5.R")

#' delete_file(): Delete a file.
#' 
#' Delete a file if it exists.
#' @param file_name File name.
#' @export
delete_file <- function(file_name)
{
  # print(paste("Deleting ", file_name)) # TODO uncomment
  if (file.exists(file_name))
  {
    # file.remove(file_name) # TODO uncomment
  }
}

#' validate_h5_sequence(): Validate H5 data for a specific sequence.
#' 
#' @param sequence: Sequence name (character, character)
#' @param h5_file: H5 file with data on sequence to be validated
#' (character, character)
#' @param gff: GFF data (tbl_df tbl data.frame, list)
#' @param bam_hdr_seq_info: Data on sequences from BAM file header
#' (GenomeInfoDb::Seqinfo, S4) 
#' @param bam: Data on alignments from BAM file
#' (GenomicAlignments::GAlignments, S4) 
#' @param dataset: Dataset name (character, character)
#' @param buffer: Buffer size (numeric, double)
#' @param min_read_length: Minimum read length (numeric, double)
#' @param max_read_length: Maximum read length (numeric, double)
#' @export
validate_h5_sequence <- function(sequence, h5_file, gff,
  bam_hdr_seq_info, bam, dataset, buffer, min_read_length,
  max_read_length) 
{

  num_read_counts <- max_read_length - min_read_length + 1

  # Get sequence positions from GFF
  gff_utr5_start <- GetCDS5start(sequence, gff, ftype="UTR5")
  gff_utr5_end <- GetCDS3end(sequence, gff, ftype="UTR5")
  gff_utr5_length <- gff_utr5_end - gff_utr5_start + 1
  gff_cds_start <- GetCDS5start(sequence, gff, ftype="CDS")
  gff_cds_end <- GetCDS3end(sequence, gff, ftype="CDS")
  gff_cds_length <- gff_cds_end - gff_cds_start + 1
  gff_utr3_start <- GetCDS5start(sequence, gff, ftype="UTR3")
  gff_utr3_end <- GetCDS3end(sequence, gff, ftype="UTR3")
  gff_utr3_length <- gff_utr3_end - gff_utr3_start + 1
  print(paste0("UTR5 start/length/end: ", gff_utr5_start, " ",
    gff_utr5_length, " ", gff_utr5_end))
  print(paste0("CDS start/length/end: ", gff_cds_start, " ",
    gff_cds_length, " ", gff_cds_end))
  print(paste0("UTR3 start/length/end: ", gff_utr3_start, " ",
    gff_utr3_length, " ", gff_utr3_end))

  # Get sequence length from BAM header
  bam_hdr_sequence <- bam_hdr_seq_info[sequence] # GenomeInfoDb::Seqinfo, S4
  bam_hdr_sequence_seq_length <- bam_hdr_sequence@seqlengths
  print(paste0("Sequence length: ", bam_hdr_sequence_seq_length))

  # Get sequence entries from BAM
  bam_sequence = bam[(GenomicAlignments::seqnames(bam) == sequence)] # GenomicAlignments::GAlignments, S4
  print(paste0("Number of alignments: ", length(bam_sequence)))
  bam_sequence_flag_zero = bam[(GenomicAlignments::seqnames(bam) == sequence)
    & (mcols(bam)$flag == 0)]
  print(paste0("Number of alignments (Flag = 0): ", length(bam_sequence_flag_zero)))
  bam_sequence_flag_non_zero =
    bam[(GenomicAlignments::seqnames(bam) == sequence)
    & (mcols(bam)$flag != 0)]
  print(paste0("Number of alignment (Flag != 0): ", length(bam_sequence_flag_non_zero)))
 
  # Validate buffer_left: number of nucleotides upstream of the start
  # codon (ATG) (UTR5 length)
  h5_buffer_left <- GetGeneBufferLeft(sequence, dataset, h5_file) # double
  print(paste0("buffer_left: ", h5_buffer_left))
  expect_equal(h5_buffer_left, buffer,
    info = "Unexpected buffer_left, compared to bam_to_h5.R parameter")
  expect_equal(h5_buffer_left, gff_utr5_length,
    info = "Unexpected buffer_left, compared to GFF UTR5 length")

  # Validate buffer_right: number of nucleotides downstream of the
  # stop codon (TAA/TAG/TGA) (UTR3 length) 
  h5_buffer_right <- GetGeneBufferRight(sequence, dataset, h5_file) # integer
  print(paste0("buffer_right: ", h5_buffer_right))
  expect_equal(h5_buffer_right, buffer,
    info = "Unexpected buffer_right, compared to bam_to_h5.R parameter")
  expect_equal(h5_buffer_left, gff_utr3_length,
    info = "Unexpected buffer_left, compared to GFF UTR3 length")

  # Validate start_codon_pos: Positions corresponding to start codon
  # of CDS in organism sequence
  gff_start_codon_pos <- as.array(seq(gff_cds_start, gff_cds_start + 2))
  h5_start_codon_pos <- GetGeneStartCodonPos(sequence, dataset, h5_file) # 1D array of 3 integer
  print(paste0("start_codon_pos: ", toString(h5_start_codon_pos)))
  expect_equal(length(h5_start_codon_pos), 3,
    info = "Unexpected start_codon_pos length")
  expect_equal(h5_start_codon_pos, gff_start_codon_pos,
    info = "Unexpected start_codon_pos, compared to GFF CDS start codon positions")

  # Validate stop_codon_pos: Positions corresponding to stop codon of
  # CDS in organism sequence
  gff_stop_codon_pos <- as.array(seq(gff_cds_end - 2, gff_cds_end))
  h5_stop_codon_pos <- GetGeneStopCodonPos(sequence, dataset, h5_file) # 1D array of 3 integer
  print(paste0("stop_codon_pos: ", toString(h5_stop_codon_pos)))
  expect_equal(length(h5_stop_codon_pos), 3,
    info = "Unexpected stop_codon_pos length")
  expect_equal(h5_stop_codon_pos, gff_stop_codon_pos,
    info = "Unexpected stop_codon_pos, compared to GFF CDS stop codon positions")

  # Validate lengths: Lengths of mapped reads.
  lengths <- as.array(seq(min_read_length, max_read_length))
  h5_lengths <- GetGeneMappedReadLengths(sequence, dataset, h5_file) # 1D array of <max_read_length - min_read_length + 1> integer
  print(paste0("lengths: ", toString(h5_lengths)))
  expect_equal(length(h5_lengths), num_read_counts,
    info = "lengths length does not equal max_read_length - min_read_length + 1")
  expect_equal(h5_lengths, lengths, info = "Unexpected lengths")

  # Validate reads_by_len: Counts of number of ribosome sequences of
  # each length
  h5_reads_by_len <- GetGeneReadLength(sequence, dataset, h5_file) # 1D array of <max_read_length - min_read_length + 1> double
  print(paste0("reads_by_len: ", toString(h5_reads_by_len)))
  expect_equal(length(h5_reads_by_len), num_read_counts,
    info = "reads_by_len length does not equal max_read_length - min_read_length + 1")

  # Calculate expected reads_by_len based on information from BAM
  reads_by_len_bam <- as.array(integer(num_read_counts))
  for (width in sort(GenomicAlignments::qwidth(bam_sequence_flag_zero)))
  {
    index <- width - min_read_length + 1
    reads_by_len_bam[index] <- reads_by_len_bam[index] + 1
  }
  expect_equal(h5_reads_by_len, reads_by_len_bam,
    info = "Unexpected reads_by_len, compared to those computed from BAM")

  # Validate reads_total: Total number of ribosome sequences
  h5_reads_total <- GetGeneReadsTotal(sequence, dataset, h5_file) # 1D array of 1 double
  print(paste0("reads_total: ", h5_reads_total))
  expect_equal(length(h5_reads_total), 1,
    info = "Unexpected reads_total length")
  h5_reads_len_total <- Reduce("+", h5_reads_total)
  expect_equal(h5_reads_total[1], h5_reads_len_total,
    info = "reads_total does not equal sum of totals in reads_by_len")
  expect_equal(h5_reads_total[1], length(bam_sequence_flag_zero),
    info = "reads_total does not equal number of BAM alignments with Flag = 0")

  # Validate data: Positions and lengths of ribosome sequences within
  # the organism data
  h5_data <- GetGeneDatamatrix(sequence, dataset, h5_file) # matrix, integer
  print(paste0("data rows/columns: ", toString(dim(h5_data))))
  num_data_cols <- h5_stop_codon_pos[3] + buffer
  expect_equal(nrow(h5_data), num_read_counts,
    info = "Number of data rows does not equal max_read_length - min_read_length + 1")
  expect_equal(ncol(h5_data), num_data_cols,
    info = "Number of data columns does not equal stop_codon_pos[3] + buffer")
  expect_equal(ncol(h5_data), gff_utr3_end,
    info = "Number of data columns does not equal GFF UTR3 final nt position")
  expect_equal(ncol(h5_data), bam_hdr_sequence_seq_length,
    info = "Number of data columns does not equal BAM sequence length")
  h5_reads_by_len_data <- rowSums(h5_data)
  expect_equal(h5_reads_by_len, as.array(h5_reads_by_len_data),
    info = "reads_by_len is not consistent with data")

  # Calculate expected data based on information from BAM
  data <- matrix(0, nrow = num_read_counts, ncol = num_data_cols)
  if (sequence %in% GenomicAlignments::seqnames(bam))
  {
    print("Sequence has alignments in BAM.")
    for (align in names(bam_sequence_flag_zero))
    {
      start <- GenomicAlignments::start(bam_sequence_flag_zero[align])
      width <- GenomicAlignments::qwidth(bam_sequence_flag_zero[align])
      width <- width - min_read_length + 1
      data[width, start] = data[width, start] + 1
    }
    expect_equal(h5_data, data,
      info = "Unexpected data, compared to that computed from BAM")
  }
  else
  {
      print("Sequence has no alignments in BAM.")
      expect_equal(h5_data, data,
        info = "Unexpected data, expected 0s as no alignments in BAM")
  }

  reads_by_len_data = rowSums(data)
  expect_equal(h5_reads_by_len_data, reads_by_len_data,
    info = "Unexpected reads_by_len length, compared to those computed from data computed from BAM")
  expect_equal(h5_reads_by_len, as.array(reads_by_len_data),
    info = "Unexpected reads_by_len, compared to those computed from data computed from BAM")

}

testthat::test_that("Run bam_to_h5.R and validate H5 file", {
  withr::defer(delete_file(h5_file))

  min_read_length <- 10
  max_read_length <- 50
  buffer <- 250
  primary_id <- "Name"
  secondary_id <- "NULL"
  dataset <- "vignette"
  is_riboviz_gff <- TRUE
  stop_in_cds <- FALSE

  bam_to_h5_cmd_template <- "Rscript --vanilla {bam_to_h5} --num-processes=1 --min-read-length={min_read_length} --max-read-length={max_read_length} --buffer={buffer} --primary-id={primary_id} --secondary-id={secondary_id} --dataset={dataset} --bam-file={bam_file} --hd-file={h5_file} --orf-gff-file={gff_file} --is-riboviz-gff={is_riboviz_gff} --stop-in-cds={stop_in_cds}"
  print(bam_to_h5_cmd_template)
  cmd <- glue(bam_to_h5_cmd_template)
  print(cmd)
  if (FALSE) # TODO uncomment
  {
  exit_code <- system(cmd)
  print(paste0("bam_to_h5.R exit code: ", exit_code))
  expect_equal(exit_code, 0, info = "Unexpected exit code from bam_to_h5.R")
  }

  ##
  ## READ GFF
  ##
  
  print("========== GFF ==========")
  gff <- readGFFAsDf(gff_file) # tbl_df tbl data.frame, list
  gff_names <- unique(gff$seqnames) # factor, integer
  print(paste0("GFF sequence names (", length(gff_names), "):"))
  print(gff_names)

  ##
  ## READ BAM
  ##
  
  print("========== BAM ==========")
  bam_file_f <- Rsamtools::BamFile(bam_file)
  bam_hdr_seq_info <- Rsamtools::seqinfo(bam_file_f) # GenomeInfoDb::Seqinfo, S4
  bam_hdr_names <- bam_hdr_seq_info@seqnames # character, character
  print(paste0("BAM header sequence names (", length(bam_hdr_names), "):"))
  print(bam_hdr_names)

  # By default readGAlignments returns: seqnames, strand, cigar,
  # qwidth, start, end, width. Also want "flag" so specify
  # explicitly.
  bam_params <- Rsamtools::ScanBamParam(what = c("flag"))
  bam <- GenomicAlignments::readGAlignments(bam_file,
    param = bam_params, use.names = T) # GenomicAlignments::GAlignments, S4
  print(paste0("Number of BAM alignments: ", length(bam)))
  bam_names <- unique(sort(GenomicAlignments::seqnames(bam))) # factor, integer
  print(paste0("BAM sequence names (", length(bam_names), "):"))
  print(bam_names)

  ##
  ## READ H5
  ##

  print("========== H5 ==========")
  h5_data <- rhdf5::h5ls(h5_file, recursive = 1) # data.frame, list
  h5_names <- h5_data$name # character, character
  print(paste0("H5 sequence names (", length(h5_names), "):"))
  print(h5_names)

  ##
  ## VALIDATE H5
  ##

  expect_equal(length(h5_names), length(gff_names),
    info = "Unexpected number of sequence names, compared to GFF")
  expect_equal(as.factor(sort(h5_names)), sort(gff_names),
    info = "Unexpected sequence names, compared to GFF")

  expect_equal(length(h5_names), length(bam_hdr_names),
    info = "Unexpected number of sequence names, compared to BAM header")
  expect_equal(sort(h5_names), sort(bam_hdr_names),
    info = "Unexpected sequence names, compared to BAM header")

  expect_true(length(bam_hdr_names) <= length(h5_names),
    info = "Number of sequence names should be >= to those in BAM")
  expect_true(all(sort(bam_names) %in% as.factor(sort(h5_names))),
    info = "Sequence names should be superset of those in BAM")

  ##
  ## VALIDATE H5 (sequence-specific)
  ##

  for (sequence in h5_names)
  {
  # print(paste0("Sequence: ", sequence))
  # validate_h5_sequence(sequence, h5_file, gff, bam_hdr_seq_info, bam, 
  # dataset, buffer, min_read_length, max_read_length)
  }
  sequence <- "YAL062W"
  validate_h5_sequence(sequence, h5_file, gff, bam_hdr_seq_info, bam,
    dataset, buffer, min_read_length, max_read_length)
  sequence <- "YAL001C" # 6 BAM (4 Flag = 0, 2 Flag != 0), 4 HDF5
  validate_h5_sequence(sequence, h5_file, gff, bam_hdr_seq_info, bam,
    dataset, buffer, min_read_length, max_read_length)
  sequence <- "YAL018C" # GFF, BAM header, H5, no BAM body.
  validate_h5_sequence(sequence, h5_file, gff, bam_hdr_seq_info, bam,
    dataset, buffer, min_read_length, max_read_length)
})
