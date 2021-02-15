# testthat tests for bam_to_h5.R
#
# This assumes the following files are in the path:
#
# rscripts/bam_to_h5.R
#
# The test runs bam_to_h5.R using the GFF and BAM file then validates
# the .h5 file created based upon its expected qualities given those
# of the input GFF and BAM files and the bam_to_h5.R command-line
# parameters.
#
# To run interactively:
#
# test_file("rscripts/tests/testthat/test_bam_to_h5.R")
#
# To run from console:
#
# Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_bam_to_h5.R")'
#
# This test assumes the following test data files exist:
#
# data/Mok-tinysim-gffsam/A.sam
# data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3
#
# To use this test with other data files see comments on
# test_that("Run bam_to_h5.R and validate H5 file"...) below.

suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(withr, quietly = T))
suppressMessages(library(GenomicAlignments, quietly = T))
suppressMessages(library(Rsamtools, quietly = T))

source(here::here("rscripts", "read_count_functions.R"))
print(paste0("here: ", here()))
bam_to_h5 <- here::here("rscripts/bam_to_h5.R")

context("test_bam_to_h5.R")

#' delete_file(): Delete a file.
#' 
#' Delete a file if it exists.
#' @param file_name File name.
#' @export
delete_file <- function(file_name)
{
  print(paste0("Deleting ", file_name))
  if (file.exists(file_name))
  {
      file.remove(file_name)
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
  bam_sequence <- bam[(GenomicAlignments::seqnames(bam) == sequence)] # GenomicAlignments::GAlignments, S4
  print(paste0("Number of alignments: ", length(bam_sequence)))
  bam_sequence_kept <- bam_sequence[
    (mcols(bam_sequence)$flag %in% c(0, 256))]
  print(paste0("Number of alignments (Flag = 0|256): ",
    length(bam_sequence_kept)))
  bam_sequence_discard <- bam_sequence[
    (!(mcols(bam_sequence)$flag %in% c(0, 256)))]
  print(paste0("Number of alignments (Flag != 0|256): ",
    length(bam_sequence_discard)))

  # Validate buffer_left: number of nucleotides upstream of the start
  # codon (ATG) (UTR5 length)
  h5_buffer_left <- GetGeneBufferLeft(sequence, dataset, h5_file) # double
  print(paste0("buffer_left: ", h5_buffer_left))
  expect_equal(h5_buffer_left, buffer,
    info = paste0(sequence,
      ": Unexpected buffer_left, compared to bam_to_h5.R parameter"))
  expect_equal(h5_buffer_left, gff_utr5_length,
    info = "Unexpected buffer_left, compared to GFF UTR5 length")

  # Validate buffer_right: number of nucleotides downstream of the
  # stop codon (TAA/TAG/TGA) (UTR3 length) 
  h5_buffer_right <- GetGeneBufferRight(sequence, dataset, h5_file) # integer
  print(paste0("buffer_right: ", h5_buffer_right))
  expect_equal(h5_buffer_right, buffer,
    info = paste0(sequence,
      ": Unexpected buffer_right, compared to bam_to_h5.R parameter"))
  expect_equal(h5_buffer_right, gff_utr3_length,
    info = paste0(sequence,
      ": Unexpected buffer_right, compared to GFF UTR3 length"))

  # Validate start_codon_pos: Positions corresponding to start codon
  # of CDS in organism sequence
  gff_start_codon_pos <- as.array(seq(gff_cds_start, gff_cds_start + 2))
  h5_start_codon_pos <- GetGeneStartCodonPos(sequence, dataset, h5_file) # 1D array of 3 integer
  print(paste0("start_codon_pos: ", toString(h5_start_codon_pos)))
  expect_equal(length(h5_start_codon_pos), 3,
    info = paste0(sequence, ": Unexpected start_codon_pos length"))
  expect_equal(h5_start_codon_pos, gff_start_codon_pos,
    info = paste0(sequence,
      ": Unexpected start_codon_pos, compared to GFF CDS start codon positions"))

  # Validate stop_codon_pos: Positions corresponding to stop codon of
  # CDS in organism sequence
  gff_stop_codon_pos <- as.array(seq(gff_cds_end - 2, gff_cds_end))
  h5_stop_codon_pos <- GetGeneStopCodonPos(sequence, dataset, h5_file) # 1D array of 3 integer
  print(paste0("stop_codon_pos: ", toString(h5_stop_codon_pos)))
  expect_equal(length(h5_stop_codon_pos), 3,
    info = paste0(sequence, ": Unexpected stop_codon_pos length"))
  expect_equal(h5_stop_codon_pos, gff_stop_codon_pos,
    info = paste0(sequence,
      ": Unexpected stop_codon_pos, compared to GFF CDS stop codon positions"))

  # Validate lengths: Lengths of mapped reads.
  lengths <- as.array(seq(min_read_length, max_read_length))
  h5_lengths <- GetGeneMappedReadLengths(sequence, dataset, h5_file) # 1D array of <max_read_length - min_read_length + 1> integer
  print(paste0("lengths: ", toString(h5_lengths)))
  expect_equal(length(h5_lengths), num_read_counts,
    info = paste0(sequence,
      ": lengths length does not equal max_read_length - min_read_length + 1"))
  expect_equal(h5_lengths, lengths,
    info = paste0(sequence, ": Unexpected lengths"))

  # Validate reads_by_len: Counts of number of ribosome sequences of
  # each length
  h5_reads_by_len <- GetGeneReadLength(sequence, dataset, h5_file) # 1D array of <max_read_length - min_read_length + 1> double
  print(paste0("reads_by_len: ", toString(h5_reads_by_len)))
  expect_equal(length(h5_reads_by_len), num_read_counts,
    info = paste0(sequence,
      ": reads_by_len length does not equal max_read_length - min_read_length + 1"))

  # Calculate expected reads_by_len based on information from BAM
  reads_by_len_bam <- as.array(integer(num_read_counts))
  for (width in sort(GenomicAlignments::qwidth(bam_sequence_kept)))
  {
    index <- width - min_read_length + 1
    reads_by_len_bam[index] <- reads_by_len_bam[index] + 1
  }
  print(paste0("reads_by_len_bam: ", toString(reads_by_len_bam)))
  expect_equal(h5_reads_by_len, reads_by_len_bam,
    info = paste0(sequence,
      ": Unexpected reads_by_len, compared to those computed from BAM"))
  
  # Validate reads_total: Total number of ribosome sequences
  h5_reads_total <- GetGeneReadsTotal(sequence, dataset, h5_file) # 1D array of 1 double
  print(paste0("reads_total: ", h5_reads_total))
  expect_equal(length(h5_reads_total), 1,
    info = paste0(sequence, ": Unexpected reads_total length"))
  h5_reads_len_total <- Reduce("+", h5_reads_total)
  expect_equal(h5_reads_total[1], h5_reads_len_total,
    info = paste0(sequence,
      ": reads_total does not equal sum of totals in reads_by_len"))
  expect_equal(h5_reads_total[1], length(bam_sequence_kept),
    info = paste0(sequence,
      ": reads_total does not equal number of BAM alignments with Flag = 0"))

  # Validate data: Positions and lengths of ribosome sequences within
  # the organism data
  h5_data <- GetGeneDatamatrix(sequence, dataset, h5_file) # matrix, integer
  print(paste0("data rows/columns: ", toString(dim(h5_data))))
  num_data_cols <- h5_stop_codon_pos[3] + buffer
  expect_equal(nrow(h5_data), num_read_counts,
    info = paste0(sequence,
      ": Number of data rows does not equal max_read_length - min_read_length + 1"))
  expect_equal(ncol(h5_data), num_data_cols,
    info = paste0(sequence,
      ": Number of data columns does not equal stop_codon_pos[3] + buffer"))
  expect_equal(ncol(h5_data), gff_utr3_end,
    info = paste0(sequence,
      ": Number of data columns does not equal GFF UTR3 final nt position"))
  expect_equal(ncol(h5_data), bam_hdr_sequence_seq_length,
    info = paste0(sequence,
      ": Number of data columns does not equal BAM sequence length"))
  h5_reads_by_len_data <- rowSums(h5_data)
  expect_equal(h5_reads_by_len, as.array(h5_reads_by_len_data),
    info = paste0(sequence, ": reads_by_len is not consistent with data"))

  # Calculate expected data based on information from BAM
  data <- matrix(0, nrow = num_read_counts, ncol = num_data_cols)
  if (sequence %in% GenomicAlignments::seqnames(bam))
  {
    print("Sequence has alignments in BAM.")
    for (align in names(bam_sequence_kept))
    {
      start <- GenomicAlignments::start(bam_sequence_kept[align])
      width <- GenomicAlignments::qwidth(bam_sequence_kept[align])
      width <- width - min_read_length + 1
      data[width, start] <- data[width, start] + 1
    }
    expect_equal(h5_data, data,
      info = paste0(sequence,
        ": Unexpected data, compared to that computed from BAM"))
  }
  else
  {
      print("Sequence has no alignments in BAM.")
      expect_equal(h5_data, data,
        info = paste0(sequence,
	  ": Unexpected data, expected 0s as no alignments in BAM"))
  }
  reads_by_len_data <- rowSums(data)
  expect_equal(h5_reads_by_len_data, reads_by_len_data,
    info = paste0(sequence,
      ": Unexpected reads_by_len length, compared to those computed from data computed from BAM"))
  expect_equal(h5_reads_by_len, as.array(reads_by_len_data),
    info = paste0(sequence,
      ": Unexpected reads_by_len, compared to those computed from data computed from BAM"))
}

#' validate_h5(): Validate H5 data.
#' 
#' @param h5_file: H5 file with data on sequence to be validated
#' (character, character)
#' @param gff_file: GFF file (character, character)
#' @param bam_file: BAM file (character, character)
#' @param dataset: Dataset name (character, character)
#' @param buffer: Buffer size (numeric, double)
#' @param min_read_length: Minimum read length (numeric, double)
#' @param max_read_length: Maximum read length (numeric, double)
#' @export
validate_h5 <- function(h5_file, gff_file, bam_file, dataset, buffer,
  min_read_length, max_read_length)
{

  # READ GFF
  
  print("========== GFF ==========")
  gff <- readGFFAsDf(gff_file) # tbl_df tbl data.frame, list
  gff_names <- unique(gff$seqnames) # factor, integer
  print(paste0("GFF sequence names (", length(gff_names), "):"))
  print(gff_names)
  
  # READ BAM
  
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

  # READ H5

  print("========== H5 ==========")
  h5_data <- rhdf5::h5ls(h5_file, recursive = 1) # data.frame, list
  h5_names <- h5_data$name # character, character
  print(paste0("H5 sequence names (", length(h5_names), "):"))
  print(h5_names)

  ## VALIDATE H5

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

  ## VALIDATE H5 (sequence-specific)

  for (sequence in h5_names)
  {
    print(paste0("Sequence: ", sequence))
    validate_h5_sequence(sequence, h5_file, gff, bam_hdr_seq_info, bam,
      dataset, buffer, min_read_length, max_read_length)
  }
}


testthat::test_that("Run bam_to_h5.R and validate H5 file", {
  withr::defer(delete_file(h5_file)) # Delete H5 when test completes.
  withr::defer(delete_file(bam_file)) # Delete H5 when test completes.
  withr::defer(delete_file(bam_bai_file)) # Delete H5 when test completes.

  # To use this test with other data files:
  # * Edit gff_file.
  # * If you have a SAM file, edit sam_file, set create_bam <- TRUE.
  # * If you have a BAM file, comment out sam_file, edit bam_file,
  #   edit_bam_bai_file, set create_bam <- FALSE.
  # * Edit variables to consistent with the configuration used to
  #   create the H5 file. 

  create_bam <- TRUE

  gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")
  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  bam_file <- here::here("test_bam_to_h5_data.bam")
  bam_bai_file <- here::here("test_bam_to_h5_data.bam.bai")
  dataset <- "Mok-tinysim"
  buffer <- 20
  min_read_length <- 10
  max_read_length <- 50
  primary_id <- "Name"
  secondary_id <- "NULL"
  is_riboviz_gff <- TRUE
  stop_in_cds <- FALSE

  h5_file <- here::here("test_bam_to_h5_data.h5")

  print(paste0("bam_to_h5.R: ", bam_to_h5))
  print(paste0("GFF: ", gff_file))

  if (create_bam)
  {
    print(paste0("SAM: ", sam_file))
    bam_cmd_template <- "samtools view -S -b {sam_file} > {bam_file}"
    bam_cmd <- glue(bam_cmd_template)
    exit_code <- system(bam_cmd)
    print(paste0("samtools view exit code: ", exit_code))
    expect_equal(exit_code, 0,
                 info = "Unexpected exit code from 'samtools view'")

    bai_cmd_template <- "samtools index {bam_file}"
    bai_cmd <- glue(bai_cmd_template)
    exit_code <- system(bai_cmd)
    print(paste0("samtools index exit code: ", exit_code))
    expect_equal(exit_code, 0,
        info = "Unexpected exit code from 'samtools index'")
  }

  print(paste0("BAM: ", bam_file))
  print(paste0("HDF5: ", h5_file))

  h5_cmd_template <- "Rscript --vanilla {bam_to_h5} --num-processes=1 --min-read-length={min_read_length} --max-read-length={max_read_length} --buffer={buffer} --primary-id={primary_id} --secondary-id={secondary_id} --dataset={dataset} --bam-file={bam_file} --hd-file={h5_file} --orf-gff-file={gff_file} --is-riboviz-gff={is_riboviz_gff} --stop-in-cds={stop_in_cds}"
  h5_cmd <- glue(h5_cmd_template)
  print(h5_cmd)

  exit_code <- system(h5_cmd)
  print(paste0("bam_to_h5.R exit code: ", exit_code))
  expect_equal(exit_code, 0,
    info = "Unexpected exit code from bam_to_h5.R")

  validate_h5(h5_file, gff_file, bam_file, dataset, buffer,
    min_read_length, max_read_length)
})
