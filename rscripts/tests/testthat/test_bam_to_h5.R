#' testthat tests for `bam_to_h5.R`.
#'
#' These tests run `bam_to_h5.R` using a GFF and SAM file (which is
#' converted into a BAM file as part of the tests), validates the HDF5
#' files created based upon its expected qualities  given those of the
#' input GFF and BAM files and the `bam_to_h5.R` command-line
#' parameters.
#'
#' The tests assumes the following files are in the path:
#'
#' ```
#' rscripts/bam_to_h5.R
#' ```
#'
#' To run the tests interactively from within R:
#'
#' ```
#' test_file("rscripts/tests/testthat/test_bam_to_h5.R")
#' ```
#'
#' To run the tests from the console:
#'
#' ```
#' Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_bam_to_h5.R")'
#' ```
#'
#' These tests assumes the following test data files exist:
#'
#' ```
#' data/Mok-tinysim-gffsam/A.sam
#' data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3
#' ```
#'
#' At present, the following behaviours are not tested:
#'
#' * BAM file specifies -ve strands.
#' * There are multiple exon genes.
#'
#' @export

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

#' Delete a file.
#'
#' Delete a file if it exists.
#' @param file_name File name (character).
#' @export
DeleteFile <- function(file_name) {
  print(paste0("Deleting ", file_name))
  if (file.exists(file_name)) {
      file.remove(file_name)
  }
}

#' Validate H5 data for a specific sequence.
#'
#' @param sequence Sequence name (character).
#' @param h5_file H5 file with sequence data to be validated
#' (character).
#' @param gff GFF data (tbl_df tbl data.frame).
#' @param bam_hdr_seq_info Data on sequences from BAM file header
#' (GenomeInfoDb::Seqinfo).
#' @param bam Data on alignments from BAM file
#' (GenomicAlignments::GAlignments).
#' @param dataset Human-readable name of the dataset (character).
#' @param buffer Length of flanking region around the CDS (integer).
#' @param min_read_length Minimum read length in H5 output (integer).
#' @param max_read_length Maximum read length in H5 output (integer).
#' @param is_riboviz_gff Does the GFF file contain 3 elements per gene
#' - UTR5, CDS, and UTR3? (logical).
#' @param stop_in_cds Are stop codons part of the CDS annotations in
#' GFF? (logical).
#'
#' @export
ValidateH5Sequence <- function(sequence, h5_file, gff,
  bam_hdr_seq_info, bam, dataset, buffer, min_read_length,
  max_read_length, is_riboviz_gff, stop_in_cds) {

  num_read_counts <- max_read_length - min_read_length + 1
  # Get sequence positions from GFF
  gff_cds_start <- GetCDS5start(sequence, gff, ftype = "CDS")
  gff_cds_end <- GetCDS3end(sequence, gff, ftype = "CDS")
  gff_cds_length <- gff_cds_end - gff_cds_start + 1
  if (is_riboviz_gff) {
    # Get sequence positions from GFF
    utr5_start <- GetCDS5start(sequence, gff, ftype = "UTR5")
    utr5_end <- GetCDS3end(sequence, gff, ftype = "UTR5")
    utr5_length <- utr5_end - utr5_start + 1
    utr3_start <- GetCDS5start(sequence, gff, ftype = "UTR3")
    utr3_end <- GetCDS3end(sequence, gff, ftype = "UTR3")
    utr3_length <- utr3_end - utr3_start + 1
    stop_codon_pos <- as.array(seq(gff_cds_end - 2, gff_cds_end))
    h5_buffer_left_info <-
      "Unexpected buffer_left, compared to GFF UTR5 length"
    h5_buffer_right_info <-
      "Unexpected buffer_right, compared to GFF UTR5 length"
    h5_stop_codon_info <-
      "Unexpected stop_codon_pos, compared to GFF CDS positions"
  } else {
    # Get sequence positions from GFF and buffer
    utr5_start <- 1
    utr5_length <- buffer
    utr5_end <- utr5_start + buffer - 1
    utr3_start <- gff_cds_end + 1
    utr3_length <- buffer
    utr3_end <- utr3_start + buffer - 1
    if (stop_in_cds) {
      offset <- 2
    } else {
      offset <- -1
    }
    stop_codon_loc <- utr3_end - buffer - offset
    stop_codon_pos <- as.array(seq(stop_codon_loc, stop_codon_loc + 2))
    h5_buffer_left_info <-
      "Unexpected buffer_left, compared to 'buffer' parameter"
    h5_buffer_right_info <-
      "Unexpected buffer_right, compared to 'buffer' parameter"
    h5_stop_codon_info <-
      "Unexpected stop_codon_pos, compared to that derived from 'buffer' parameter"
  }
  print(paste0("UTR5 start/length/end: ", utr5_start, " ",
    utr5_length, " ", utr5_end))
  print(paste0("CDS start/length/end: ", gff_cds_start, " ",
    gff_cds_length, " ", gff_cds_end))
  print(paste0("UTR3 start/length/end: ", utr3_start, " ",
    utr3_length, " ", utr3_end))

  # Get sequence length from BAM header
  bam_hdr_sequence <- bam_hdr_seq_info[sequence] # GenomeInfoDb::Seqinfo, S4
  bam_hdr_sequence_seq_length <- bam_hdr_sequence@seqlengths
  print(paste0("Sequence length: ", bam_hdr_sequence_seq_length))

  # Get sequence entries from BAM
  # GenomicAlignments::GAlignments, S4
  bam_sequence <- bam[(GenomicAlignments::seqnames(bam) == sequence)]
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
  expect_equal(h5_buffer_left, utr5_length,
    info = paste0(sequence, ": ", h5_buffer_left_info))

  # Validate buffer_right: number of nucleotides downstream of the
  # stop codon (TAA/TAG/TGA) (UTR3 length)
  h5_buffer_right <- GetGeneBufferRight(sequence, dataset, h5_file) # integer
  print(paste0("buffer_right: ", h5_buffer_right))
  buffer_right <- utr3_end - stop_codon_pos[3]
  expect_equal(h5_buffer_right, buffer_right,
    info = paste0(sequence, ": ", h5_buffer_right_info))

  # Validate start_codon_pos: Positions corresponding to start codon
  # of CDS in organism sequence
  gff_start_codon_pos <- as.array(seq(gff_cds_start, gff_cds_start + 2))
  h5_start_codon_pos <-
    GetGeneStartCodonPos(sequence, dataset, h5_file) # 1D array of 3 integer
  print(paste0("start_codon_pos: ", toString(h5_start_codon_pos)))
  expect_equal(length(h5_start_codon_pos), 3,
    info = paste0(sequence, ": Unexpected start_codon_pos length"))
  expect_equal(h5_start_codon_pos, gff_start_codon_pos,
    info = paste0(sequence,
      ": Unexpected start_codon_pos, compared to GFF CDS positions"))

  # Validate stop_codon_pos: Positions corresponding to stop codon of
  # CDS in organism sequence
  h5_stop_codon_pos <-
    GetGeneStopCodonPos(sequence, dataset, h5_file) # 1D array of 3 integer
  print(paste0("stop_codon_pos: ", toString(h5_stop_codon_pos)))
  expect_equal(length(h5_stop_codon_pos), 3,
    info = paste0(sequence, ": Unexpected stop_codon_pos length"))
  expect_equal(h5_stop_codon_pos, stop_codon_pos,
    info = paste0(sequence, ": ", h5_stop_codon_info))

  # Validate lengths: Lengths of mapped reads.
  lengths <- as.array(seq(min_read_length, max_read_length))
  # 1D array of <max_read_length - min_read_length + 1> integer
  h5_lengths <- GetGeneMappedReadLengths(sequence, dataset, h5_file)
  print(paste0("lengths: ", toString(h5_lengths)))
  expect_equal(length(h5_lengths), num_read_counts,
    info = paste0(sequence,
      ": lengths length != max_read_length - min_read_length + 1"))
  expect_equal(h5_lengths, lengths,
    info = paste0(sequence, ": Unexpected lengths"))

  # Validate reads_by_len: Counts of number of ribosome sequences of
  # each length
  # 1D array of <max_read_length - min_read_length + 1> double
  h5_reads_by_len <- GetGeneReadLength(sequence, dataset, h5_file)
  print(paste0("reads_by_len: ", toString(h5_reads_by_len)))
  expect_equal(length(h5_reads_by_len), num_read_counts,
    info = paste0(sequence,
      ": reads_by_len length != max_read_length - min_read_length + 1"))

  # Calculate expected reads_by_len based on information from BAM
  reads_by_len_bam <- as.array(integer(num_read_counts))
  for (width in sort(GenomicAlignments::qwidth(bam_sequence_kept))) {
    index <- width - min_read_length + 1
    reads_by_len_bam[index] <- reads_by_len_bam[index] + 1
  }
  print(paste0("reads_by_len_bam: ", toString(reads_by_len_bam)))
  expect_equal(h5_reads_by_len, reads_by_len_bam,
    info = paste0(sequence,
      ": Unexpected reads_by_len, compared to those computed from BAM"))

  # Validate reads_total: Total number of ribosome sequences
  h5_reads_total <-
    GetGeneReadsTotal(sequence, dataset, h5_file) # 1D array of 1 double
  print(paste0("reads_total: ", h5_reads_total))
  expect_equal(length(h5_reads_total), 1,
    info = paste0(sequence, ": Unexpected reads_total length"))
  h5_reads_len_total <- Reduce("+", h5_reads_total)
  expect_equal(h5_reads_total[1], h5_reads_len_total,
    info = paste0(sequence,
      ": reads_total != sum of totals in reads_by_len"))
  expect_equal(h5_reads_total[1], length(bam_sequence_kept),
    info = paste0(sequence,
      ": reads_total != number of BAM alignments with Flag = 0"))

  # Validate data: Positions and lengths of ribosome sequences within
  # the organism data
  h5_data <- GetGeneDatamatrix(sequence, dataset, h5_file) # matrix, integer
  print(paste0("data rows/columns: ", toString(dim(h5_data))))
  expect_equal(nrow(h5_data), num_read_counts,
    info = paste0(sequence,
      ": Number of data rows != max_read_length - min_read_length + 1"))
  expect_equal(ncol(h5_data), utr3_end,
    info = paste0(sequence,
      ": Number of data columns != GFF UTR3 final nt position"))
  expect_equal(ncol(h5_data), bam_hdr_sequence_seq_length,
    info = paste0(sequence,
      ": Number of data columns != BAM sequence length"))
  h5_reads_by_len_data <- rowSums(h5_data)
  expect_equal(h5_reads_by_len, as.array(h5_reads_by_len_data),
    info = paste0(sequence, ": reads_by_len is not consistent with data"))

  # Calculate expected data based on information from BAM
  data <- matrix(0, nrow = num_read_counts, ncol = utr3_end)
  if (sequence %in% GenomicAlignments::seqnames(bam)) {
    print("Sequence has alignments in BAM.")
    for (align in names(bam_sequence_kept)) {
      start <- GenomicAlignments::start(bam_sequence_kept[align])
      width <- GenomicAlignments::qwidth(bam_sequence_kept[align])
      width <- width - min_read_length + 1
      data[width, start] <- data[width, start] + 1
    }
    expect_equal(h5_data, data,
      info = paste0(sequence,
        ": Unexpected data, compared to that computed from BAM"))
  } else {
      print("Sequence has no alignments in BAM.")
      expect_equal(h5_data, data,
        info = paste0(sequence,
          ": Unexpected data, expected 0s as no alignments in BAM"))
  }
  reads_by_len_data <- rowSums(data)
  expect_equal(h5_reads_by_len_data, reads_by_len_data,
    info = paste0(sequence,
      ": Unexpected reads_by_len length, compared to those computed from BAM"))
  expect_equal(h5_reads_by_len, as.array(reads_by_len_data),
    info = paste0(sequence,
      ": Unexpected reads_by_len, compared to those computed from BAM"))
}

#' Validate H5 data within a RiboViz H5 data file.
#'
#' @param h5_file H5 file with data on sequence to be validated
#' (character).
#' @param gff_file GFF file (character).
#' @param bam_file BAM file (character).
#' @param dataset Human-readable name of the dataset (character).
#' @param buffer Length of flanking region around the CDS (integer).
#' @param min_read_length Minimum read length in H5 output (integer).
#' @param max_read_length Maximum read length in H5 output (integer).
#' @param is_riboviz_gff Does the GFF file contain 3 elements per gene
#' - UTR5, CDS, and UTR3? (logical).
#' @param stop_in_cds Are stop codons part of the CDS annotations in
#' GFF? (logical).
#' @export
ValidateH5 <- function(h5_file, gff_file, bam_file, dataset, buffer,
  min_read_length, max_read_length, is_riboviz_gff, stop_in_cds) {

  gff <- readGFFAsDf(gff_file) # tbl_df tbl data.frame, list
  gff_names <- unique(gff$seqnames) # factor, integer
  print(paste0("GFF sequence names (", length(gff_names), "):"))
  print(gff_names)

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

  for (sequence in h5_names) {
    print(paste0("Sequence: ", sequence))
    ValidateH5Sequence(sequence, h5_file, gff, bam_hdr_seq_info, bam,
      dataset, buffer, min_read_length, max_read_length,
      is_riboviz_gff, stop_in_cds)
  }
}

#' Run `bam_to_h5.R` on a SAM file.
#' 
#' `samtools view` and `samtools index` are first run to convert
#' the SAM file into a BAM file.
#'
#' @param bam_to_h5 `bam_to_h5.R` path (character).
#' @param sam_file SAM input file (character).
#' @param bam_file BAM file (character).
#' @param orf_gff_file GFF2/GFF3 Matched genome feature file,
#' specifying coding sequences locations (start and stop coordinates)
#' within the transcripts (GTF/GFF3 file) (character).
#' @param h5_file H5 output file (character).
#' @param min_read_length Minimum read length in H5 output (integer).
#' @param max_read_length Maximum read length in H5 output (integer).
#' @param buffer Length of flanking region around the CDS (integer).
#' @param primary_id Primary gene IDs to access the data (character).
#' @param secondary_id Secondary gene IDs to access the data (character).
#' @param dataset Human-readable name of the dataset (character).
#' @param is_riboviz_gff Does the GFF file contain 3 elements per gene
#' - UTR5, CDS, and UTR3? (logical).
#' @param stop_in_cds Are stop codons part of the CDS annotations in
#' GFF? (logical).
#' @param num_processes Number of processes to parallelize over
#' (integer).
#'
#' @export
RunSamToBamToH5 <- function(bam_to_h5, sam_file, bam_file,
  orf_gff_file, h5_file, min_read_length, max_read_length, buffer,
  primary_id, secondary_id, dataset, is_riboviz_gff, stop_in_cds,
  num_processes = 1) {

  print(paste0("bam_to_h5.R: ", bam_to_h5))
  print(paste0("GFF: ", orf_gff_file))
  print(paste0("SAM: ", sam_file))

  bam_cmd_template <- "samtools view -S -b {sam_file} > {bam_file}"
  bam_cmd <- glue(bam_cmd_template)
  print(bam_cmd)
  exit_code <- system(bam_cmd)
  print(paste0("'samtools view' exit code: ", exit_code))
  expect_equal(exit_code, 0,
               info = "Unexpected exit code from 'samtools view'")

  bai_cmd_template <- "samtools index {bam_file}"
  bai_cmd <- glue(bai_cmd_template)
  print(bai_cmd)
  exit_code <- system(bai_cmd)
  print(paste0("'samtools index' exit code: ", exit_code))
  expect_equal(exit_code, 0,
      info = "Unexpected exit code from 'samtools index'")

  print(paste0("BAM: ", bam_file))
  print(paste0("HDF5: ", h5_file))

  h5_cmd_template <- "Rscript --vanilla {bam_to_h5} --num-processes={num_processes} --min-read-length={min_read_length} --max-read-length={max_read_length} --buffer={buffer} --primary-id={primary_id} --secondary-id={secondary_id} --dataset={dataset} --bam-file={bam_file} --hd-file={h5_file} --orf-gff-file={orf_gff_file} --is-riboviz-gff={is_riboviz_gff} --stop-in-cds={stop_in_cds}" # nolint
  h5_cmd <- glue(h5_cmd_template)
  print(h5_cmd)

  exit_code <- system(h5_cmd)
  print(paste0("'bam_to_h5.R' exit code: ", exit_code))
  expect_equal(exit_code, 0,
    info = "Unexpected exit code from 'bam_to_h5.R'")
}

testthat::test_that(
  "Test bam_to_h5.R (is_riboviz_gff=TRUE)", {

  withr::defer(DeleteFile(h5_file)) # Delete H5 when test completes.
  withr::defer(DeleteFile(bam_file)) # Delete H5 when test completes.
  withr::defer(DeleteFile(bam_bai_file)) # Delete H5 when test completes.

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  bam_file <- here::here("test_bam_to_h5_data.bam")
  bam_bai_file <- here::here("test_bam_to_h5_data.bam.bai")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")

  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- "NULL"
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- TRUE
  stop_in_cds <- FALSE
  h5_file <- here::here("test_bam_to_h5_data.h5")
  num_processes <- 1

  RunSamToBamToH5(bam_to_h5, sam_file, bam_file, orf_gff_file, h5_file,
    min_read_length, max_read_length, buffer, primary_id, secondary_id,
    dataset, is_riboviz_gff, stop_in_cds, num_processes)

  ValidateH5(h5_file, orf_gff_file, bam_file, dataset, buffer,
    min_read_length, max_read_length, is_riboviz_gff, stop_in_cds)
})


testthat::test_that(
  "Test bam_to_h5.R (is_riboviz_gff=FALSE, stop_in_CDS=FALSE)", {

  withr::defer(DeleteFile(h5_file)) # Delete H5 when test completes.
  withr::defer(DeleteFile(bam_file)) # Delete H5 when test completes.
  withr::defer(DeleteFile(bam_bai_file)) # Delete H5 when test completes.

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  bam_file <- here::here("test_bam_to_h5_data.bam")
  bam_bai_file <- here::here("test_bam_to_h5_data.bam.bai")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")

  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- "NULL"
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- FALSE
  stop_in_cds <- FALSE
  h5_file <- here::here("test_bam_to_h5_data.h5")
  num_processes <- 1

  RunSamToBamToH5(bam_to_h5, sam_file, bam_file, orf_gff_file, h5_file,
    min_read_length, max_read_length, buffer, primary_id, secondary_id,
    dataset, is_riboviz_gff, stop_in_cds, num_processes)

  ValidateH5(h5_file, orf_gff_file, bam_file, dataset, buffer,
    min_read_length, max_read_length, is_riboviz_gff, stop_in_cds)
})

testthat::test_that(
  "Test bam_to_h5.R (is_riboviz_gff=FALSE, stop_in_CDS=TRUE)", {

  withr::defer(DeleteFile(h5_file)) # Delete H5 when test completes.
  withr::defer(DeleteFile(bam_file)) # Delete H5 when test completes.
  withr::defer(DeleteFile(bam_bai_file)) # Delete H5 when test completes.

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  bam_file <- here::here("test_bam_to_h5_data.bam")
  bam_bai_file <- here::here("test_bam_to_h5_data.bam.bai")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")

  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- "NULL"
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- FALSE
  stop_in_cds <- TRUE
  h5_file <- here::here("test_bam_to_h5_data.h5")
  num_processes <- 1

  RunSamToBamToH5(bam_to_h5, sam_file, bam_file, orf_gff_file, h5_file,
    min_read_length, max_read_length, buffer, primary_id, secondary_id,
    dataset, is_riboviz_gff, stop_in_cds, num_processes)

  ValidateH5(h5_file, orf_gff_file, bam_file, dataset, buffer,
    min_read_length, max_read_length, is_riboviz_gff, stop_in_cds)
})
