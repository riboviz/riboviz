#' testthat tests for `collate_tpms.R`.
#'
#' These tests run `collate_tpms.R` using a sample-specific TPMs files
#' then validate the collated TPMs file output.
#'
#' The tests assume the following files are in the path:
#'
#' ```
#' rscripts/collate_tpms.R
#' ```
#'
#' To run the tests interactively from within R:
#'
#' ```
#' test_file("rscripts/tests/testthat/test_collate_tpms.R")
#' ```
#'
#' To run the tests from the console:
#'
#' ```
#' Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_collate_tpms.R")'
#' ```
#'
#' @export

suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(readr, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(tibble, quietly = T))
suppressMessages(library(withr, quietly = T))

print(paste0("here: ", here()))
collate_tpms <- here::here(file.path("rscripts", "collate_tpms.R"))

context("test_collate_tpms.R")

#' Run `collate_tpms.R` on sample-specific TSV files.
#'
#' @param collate_tpms `collate_tpms.R` path (character).
#' @param output_dir Output directory in which to find sample-specific
#' TPMs files and into which to write the collated TPMs file (character).
#' @param tpms_file Name of collated TPMs file, relative to
#' `output_dir` (character).
#' @param sample_subdirs Are sample-specific TPMs files in
#' sample-specific sub-directories, in files
#' `<output_dir>/<sample>/tpms.tsv`? If not then it is assumed
#' they are in files `<output_dir>/<sample>_tpms.tsv` (logical).
#' @param orf_fasta ORF FASTA file that was aligned to and from which
#' ORF names are to be retrieved (character).
#' @param samples List of sample names (character).
#'
#' @export
RunCollateTpms <- function(collate_tpms, output_dir, tpms_file,
  sample_subdirs, orf_fasta, samples) {

  print(paste0("collate_tpms.R: ", collate_tpms))

  cmd_template <- "Rscript --vanilla {collate_tpms} --sample-subdirs={sample_subdirs} --output-dir={output_dir}"
  if (!is.na(tpms_file)) {
    cmd_template <- paste(cmd_template, "--tpms-file={tpms_file}")
  }
  if (!is.na(orf_fasta)) {
    cmd_template <- paste(cmd_template, "--orf-fasta={orf_fasta}")
  }
  # TODO equivalent of Python join?
  for (sample in samples) {
    cmd_template <- paste(cmd_template, sample)
  }
  cmd <- glue::glue(cmd_template)
  print(cmd)

  exit_code <- system(cmd)
  print(paste0("'collate_tpms.R' exit code: ", exit_code))
  testthat::expect_equal(exit_code, 0,
    info = "Unexpected exit code from 'collate_tpms.R'")
}

#' Get sample-specific TPMs.
#'
#' @param orfs List of ORFs (character).
#' @param samples List of sample names (character).
#' @return List of sample data (lists of integer of length `orfs`),
#' keyed by sample name (list).
#'
#' @export
GetSampleTpms <- function(orfs, samples) {
  offset <- 1
  sample_data <- list()
  for (sample in samples) {
    new_offset <- offset + length(orfs)
    data <- offset:(new_offset - 1)
    offset <- new_offset
    sample_data[[sample]] <- data
  }
  return(sample_data)
}

#' Get expected collated TPMs.
#'
#' @param orfs List of ORFs (character).
#' @param sample_data List of sample data (lists of integer of length
#' `orfs`),  keyed by sample name (list).
#' @return Tibble with an `ORF` column and one column for each
#' sample in `sample_data` with the sample name as column name
#' (tbl_df).
#'
#' @export
GetCollatedTpms <- function(orfs, sample_data) {
  expected <- tibble(ORF = orfs)
  for (sample in names(sample_data)) {
    expected <- add_column(expected, "{sample}" := sample_data[[sample]])
  }
  return(expected)
}

testthat::test_that("collate_tpms.R: Python workflow", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W")
  samples <- c("WTnone", "WT3AT")
  sample_data <- GetSampleTpms(orfs, samples)
  expected <- GetCollatedTpms(orfs, sample_data)
  withr::with_tempdir({
   for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_dir <- file.path(getwd(), sample)
      dir.create(sample_dir)
      sample_file <- file.path(sample_dir, "tpms.tsv")
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    RunCollateTpms(
      collate_tpms = collate_tpms,
      output_dir = output_dir,
      tpms_file = NA,
      sample_subdirs = TRUE,
      orf_fasta = NA,
      samples = samples)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: Nextflow workflow", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W")
  samples <- c("WTnone", "WT3AT")
  sample_data <- GetSampleTpms(orfs, samples)
  expected <- GetCollatedTpms(orfs, sample_data)
  withr::with_tempdir({
   for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    RunCollateTpms(
      collate_tpms = collate_tpms,
      output_dir = output_dir,
      tpms_file = NA,
      sample_subdirs = FALSE,
      orf_fasta = NA,
      samples = samples)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  # read_tsv returns spec_tbl_df subclass so cast.
  # Suggested by https://github.com/tidyverse/dplyr/issues/5126
  # Alternative is to subset with no arguments e.g. actual[]
  # https://www.tidyverse.org/blog/2018/12/readr-1-3-1/
  actual <- tibble::as_tibble(actual)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: Single sample", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W")
  samples <- c("WTnone")
  sample_data <- GetSampleTpms(orfs, samples)
  expected <- GetCollatedTpms(orfs, sample_data)
  withr::with_tempdir({
   for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    RunCollateTpms(
      collate_tpms = collate_tpms,
      output_dir = output_dir,
      tpms_file = NA,
      sample_subdirs = FALSE,
      orf_fasta = NA,
      samples = samples)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: tpms_file", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W")
  samples <- c("WTnone", "WT3AT")
  sample_data <- GetSampleTpms(orfs, samples)
  expected <- GetCollatedTpms(orfs, sample_data)
  withr::with_tempdir({
   for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    tpms_file <- "collated_tpms.tsv"
    RunCollateTpms(
      collate_tpms = collate_tpms,
      output_dir = output_dir,
      tpms_file = tpms_file,
      sample_subdirs = FALSE,
      orf_fasta = NA,
      samples = samples)
    actual_tpms_file <- file.path(output_dir, tpms_file)
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

