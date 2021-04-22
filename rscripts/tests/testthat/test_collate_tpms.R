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

suppressMessages(library(Biostrings, quietly = T))
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
#' @param samples List of sample names (character).
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
#'
#' @export
RunCollateTpms <- function(collate_tpms, samples, output_dir,
  tpms_file = NA, sample_subdirs = FALSE, orf_fasta = NA) {

  print(paste0("collate_tpms.R: ", collate_tpms))

  cmd_template <- "Rscript --vanilla {collate_tpms} --sample-subdirs={sample_subdirs} --output-dir={output_dir}"
  # TODO equivalent of Python join?
  for (sample in samples) {
    cmd_template <- paste(cmd_template, sample)
  }
  if (!is.na(tpms_file)) {
    cmd_template <- paste(cmd_template, "--tpms-file={tpms_file}")
  }
  if (!is.na(orf_fasta)) {
    cmd_template <- paste(cmd_template, "--orf-fasta={orf_fasta}")
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
  expected <- tibble::tibble(ORF = orfs)
  for (sample in names(sample_data)) {
    expected <- add_column(expected, "{sample}" := sample_data[[sample]])
  }
  return(expected)
}

testthat::test_that("collate_tpms.R: Default", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTnone", "WT3AT")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  # read_tsv returns spec_tbl_df subclass so cast.
  # Suggested by https://github.com/tidyverse/dplyr/issues/5126
  # Alternative is to subset with no arguments e.g. actual[]
  # https://www.tidyverse.org/blog/2018/12/readr-1-3-1/
  actual <- tibble::as_tibble(actual)
  expected <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: --sample-subdirs = TRUE", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTnone", "WT3AT")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_dir <- file.path(getwd(), sample)
      dir.create(sample_dir)
      sample_file <- file.path(sample_dir, "tpms.tsv")
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir, sample_subdirs = TRUE)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  expected <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: One sample", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTnone")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  expected <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: --tpms-file", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTnone", "WT3AT")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    tpms_file <- "collated_tpms.tsv"
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir, tpms_file = tpms_file)
    actual_tpms_file <- file.path(output_dir, tpms_file)
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  expected <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: Header-only sample files", {
  samples <- c("WTnone", "WT3AT")
  expected <- tibble(ORF = character())
  withr::with_tempdir({
    for (sample in samples) {
      data <- tibble(ORF = character(), tpm = double())
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
      # Expected column types are character as collated TPMs TSV file
      # will have only a header.
      expected <- add_column(expected, "{sample}" := character())
    }
    output_dir <- getwd()
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: Sample-specific ORFs in different order", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTnone", "WT3AT", "WTother")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    # For every 2nd sample reverse the order of ORFs/TPMs rows.
    # Results consistent with `expected` are sill expected.
    orf_decreasing <- FALSE
    for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      data <- data[order(data$ORF, decreasing = orf_decreasing), ]
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
      orf_decreasing <- ! orf_decreasing
    }
    output_dir <- getwd()
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  expected <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: --orf-fasta", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  fasta <- Biostrings::DNAStringSet(x = replicate(length(orfs), "GATTACCA"))
  names(fasta) <- orfs
  samples <- c("WTnone", "WT3AT", "WTother")
  sample_data <- GetSampleTpms(orfs, samples)
  expected <- GetCollatedTpms(orfs, sample_data)
  withr::with_tempdir({
    orf_fasta <- "test.fasta"
    Biostrings::writeXStringSet(fasta, orf_fasta, format = "fasta")
    # For every 2nd sample reverse the order of ORFs/TPMs rows.
    # Results consistent with `expected` are sill expected.
    orf_decreasing <- FALSE
    for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      data <- data[order(data$ORF, decreasing = orf_decreasing), ]
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
      orf_decreasing <- ! orf_decreasing
    }
    output_dir <- getwd()
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir, orf_fasta = orf_fasta)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  expected <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

testthat::test_that("collate_tpms.R: Missing sample file", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTnone")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    for (sample in names(sample_data)) {
      data <- tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    # No file exists for WT3AT. ORFs will be taken from WTnone.
    samples <- c("WTnone", "WT3AT")
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  expected <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})

"' bbb
testthat::test_that("collate_tpms.R: Missing sample files and --orf-fasta", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  fasta <- Biostrings::DNAStringSet(x = replicate(length(orfs), "GATTACCA"))
  names(fasta) <- orfs
  withr::with_tempdir({
    orf_fasta <- "test.fasta"
    Biostrings::writeXStringSet(fasta, orf_fasta, format = "fasta")
    output_dir <- getwd()
    # No file exists for WT3AT or WTNone. ORFs will be taken from orf_fasta.
    samples <- c("WTnone", "WT3AT")
    RunCollateTpms(samples = samples, collate_tpms = collate_tpms,
      output_dir = output_dir, orf_fasta = orf_fasta)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  expected <- tibble::tibble(ORF = orfs)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})
