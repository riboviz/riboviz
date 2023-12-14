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
#' * BAM file with multiple exon genes.
#'
#' @export

suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(withr, quietly = T))

source(here::here("rscripts", "tests", "testthat", "bam_to_h5_helpers.R"))
bam_to_h5 <- here::here("rscripts/bam_to_h5.R")

context("test_bam_to_h5.R")

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
#' @param primary_id Primary gene IDs to access the data (character).
#' @param secondary_id Secondary gene IDs to access the data (character).
#' @param dataset Human-readable name of the dataset (character).
#' @param min_read_length Minimum read length in H5 output (integer).
#' @param max_read_length Maximum read length in H5 output (integer).
#' @param buffer Length of flanking region around the feature (integer).
#' @param is_riboviz_gff Does the GFF file contain 3 elements per gene
#' - UTR5, feature, and UTR3? (logical).
#' @param stop_in_feature Are stop codons part of the feature
#' annotations in GFF? (logical).
#' @param feature Feature e.g. `CDS`, `ORF`, or `uORF` (character).
#' @param num_processes Number of processes to parallelize over
#' (integer).
#' @param exit_code Expected exit code (integer).
#'
#' @export
RunSamToBamToH5 <- function(bam_to_h5, sam_file, bam_file,
  orf_gff_file, h5_file, primary_id, secondary_id, dataset,
  min_read_length, max_read_length, buffer, is_riboviz_gff,
  stop_in_feature, feature = "CDS", num_processes = 1,
  exit_code = 0) {

  bam_cmd_template <- "samtools view -S -b {sam_file} > {bam_file}"
  bam_cmd <- glue::glue(bam_cmd_template)
  print(bam_cmd)
  actual_exit_code <- system(bam_cmd)
  print(paste0("'samtools view' exit code: ", actual_exit_code))
  testthat::expect_equal(actual_exit_code, 0,
               info = "Unexpected exit code from 'samtools view'")

  bai_cmd_template <- "samtools index {bam_file}"
  bai_cmd <- glue::glue(bai_cmd_template)
  print(bai_cmd)
  actual_exit_code <- system(bai_cmd)
  print(paste0("'samtools index' exit code: ", actual_exit_code))
  testthat::expect_equal(actual_exit_code, 0,
      info = "Unexpected exit code from 'samtools index'")

  h5_cmd_template <- "Rscript --vanilla {bam_to_h5} --num-processes={num_processes} --min-read-length={min_read_length} --max-read-length={max_read_length} --buffer={buffer} --primary-id={primary_id} --dataset={dataset} --bam-file={bam_file} --hd-file={h5_file} --orf-gff-file={orf_gff_file} --is-riboviz-gff={is_riboviz_gff} --stop-in-feature={stop_in_feature} --feature={feature}" # nolint
  if (!is.na(secondary_id)) {
      h5_cmd_template <- paste(h5_cmd_template,
                               "--secondary-id={secondary_id}")
  }
  h5_cmd <- glue::glue(h5_cmd_template)
  print(h5_cmd)

  actual_exit_code <- system(h5_cmd)
  print(paste0("'bam_to_h5.R' exit code: ", actual_exit_code))
  testthat::expect_equal(actual_exit_code, exit_code,
    info = "Unexpected exit code from 'bam_to_h5.R'")
}

testthat::test_that("Default", {

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")
  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- NA
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- TRUE
  stop_in_feature <- FALSE

  withr::with_tempdir({
    bam_file <- file.path(getwd(), "test_bam_to_h5_data.bam")
    bam_bai_file <- file.path(getwd(), "test_bam_to_h5_data.bam.bai")
    h5_file <- file.path(getwd(), "test_bam_to_h5_data.h5")

    RunSamToBamToH5(bam_to_h5, sam_file, bam_file, orf_gff_file, h5_file,
      primary_id, secondary_id, dataset, min_read_length, max_read_length,
      buffer, is_riboviz_gff, stop_in_feature)
    ValidateH5(h5_file, orf_gff_file, bam_file, primary_id, secondary_id,
      dataset, min_read_length, max_read_length, buffer,
      is_riboviz_gff, stop_in_feature)
  })
})

testthat::test_that("--is-riboviz-gff=FALSE", {

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")
  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- NA
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- FALSE
  stop_in_feature <- FALSE

  withr::with_tempdir({
    bam_file <- file.path(getwd(), "test_bam_to_h5_data.bam")
    bam_bai_file <- file.path(getwd(), "test_bam_to_h5_data.bam.bai")
    h5_file <- file.path(getwd(), "test_bam_to_h5_data.h5")

    RunSamToBamToH5(bam_to_h5, sam_file, bam_file, orf_gff_file, h5_file,
      primary_id, secondary_id, dataset, min_read_length, max_read_length,
      buffer, is_riboviz_gff, stop_in_feature)
    ValidateH5(h5_file, orf_gff_file, bam_file, primary_id, secondary_id,
      dataset, min_read_length, max_read_length, buffer,
      is_riboviz_gff, stop_in_feature)
  })
})

testthat::test_that("--is-riboviz-gff=FALSE, --stop-in-feature=TRUE", {

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")
  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- NA
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- FALSE
  stop_in_feature <- TRUE

  withr::with_tempdir({
    bam_file <- file.path(getwd(), "test_bam_to_h5_data.bam")
    bam_bai_file <- file.path(getwd(), "test_bam_to_h5_data.bam.bai")
    h5_file <- file.path(getwd(), "test_bam_to_h5_data.h5")

    RunSamToBamToH5(bam_to_h5, sam_file, bam_file, orf_gff_file, h5_file,
      primary_id, secondary_id, dataset, min_read_length, max_read_length,
      buffer, is_riboviz_gff, stop_in_feature)
    ValidateH5(h5_file, orf_gff_file, bam_file, primary_id, secondary_id,
      dataset, min_read_length, max_read_length, buffer,
      is_riboviz_gff, stop_in_feature)
  })
})

testthat::test_that("--feature=ORF", {

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")
  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- NA
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- TRUE
  stop_in_feature <- FALSE
  feature <- "ORF"

  withr::with_tempdir({
    test_gff_file <- file.path(getwd(), "test_bam_to_h5_data.gff3")
    bam_file <- file.path(getwd(), "test_bam_to_h5_data.bam")
    bam_bai_file <- file.path(getwd(), "test_bam_to_h5_data.bam.bai")
    h5_file <- file.path(getwd(), "test_bam_to_h5_data.h5")

    gff <- rtracklayer::import(orf_gff_file)
    levels(gff$type)[levels(gff$type) == "CDS"] <- "ORF"
    # Set phase to avoid export warning:
    # The phase information is missing for some CDS. The written file
    # will contain some CDS with no phase information.
    gff$phase <- numeric(length(gff$Name))
    rtracklayer::export(gff, test_gff_file)

    RunSamToBamToH5(bam_to_h5, sam_file, bam_file, test_gff_file, h5_file,
      primary_id, secondary_id, dataset, min_read_length, max_read_length,
      buffer, is_riboviz_gff, stop_in_feature, feature)
    ValidateH5(h5_file, test_gff_file, bam_file, primary_id, secondary_id,
      dataset, min_read_length, max_read_length, buffer,
      is_riboviz_gff, stop_in_feature, feature)
  })
})

testthat::test_that("--secondary-id=ID", {

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")
  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- "ID"
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- TRUE
  stop_in_feature <- FALSE

  withr::with_tempdir({
    test_gff_file <- file.path(getwd(), "test_bam_to_h5_data.gff3")
    bam_file <- file.path(getwd(), "test_bam_to_h5_data.bam")
    bam_bai_file <- file.path(getwd(), "test_bam_to_h5_data.bam.bai")
    h5_file <- file.path(getwd(), "test_bam_to_h5_data.h5")

    gff <- rtracklayer::import(orf_gff_file)
    gff$Name <- replace(gff$Name, TRUE, paste0("X", gff$Name))
    gff$ID <- replace(gff$Name, TRUE, paste0("X", reverse(gff$Name)))
    # Set phase to avoid export warning:
    # The phase information is missing for some CDS. The written file
    # will contain some CDS with no phase information.
    gff$phase <- numeric(length(gff$Name))
    rtracklayer::export(gff, test_gff_file)

    RunSamToBamToH5(bam_to_h5, sam_file, bam_file, test_gff_file, h5_file,
      primary_id, secondary_id, dataset, min_read_length, max_read_length,
      buffer, is_riboviz_gff, stop_in_feature)
    ValidateH5(h5_file, test_gff_file, bam_file, primary_id, secondary_id,
      dataset, min_read_length, max_read_length, buffer,
      is_riboviz_gff, stop_in_feature)
  })
})

testthat::test_that("feature=Unknown", {

  sam_file <- here::here("data/Mok-tinysim-gffsam/A.sam")
  orf_gff_file <- here::here("data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3")
  min_read_length <- 10
  max_read_length <- 50
  buffer <- 20
  primary_id <- "Name"
  secondary_id <- NA
  dataset <- "Mok-tinysim"
  is_riboviz_gff <- TRUE
  stop_in_feature <- FALSE
  feature <- "Unknown"
  exit_code <- 1

  withr::with_tempdir({
    bam_file <- file.path(getwd(), "test_bam_to_h5_data.bam")
    bam_bai_file <- file.path(getwd(), "test_bam_to_h5_data.bam.bai")
    h5_file <- file.path(getwd(), "test_bam_to_h5_data.h5")

    RunSamToBamToH5(bam_to_h5, sam_file, bam_file, orf_gff_file, h5_file,
      primary_id, secondary_id, dataset, min_read_length, max_read_length,
      buffer, is_riboviz_gff, stop_in_feature, feature,
      exit_code = exit_code)
  })
})
