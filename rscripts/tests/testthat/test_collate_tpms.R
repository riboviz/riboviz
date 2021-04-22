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
#' @return Tibble with collated TPMs (tbl_df).
#'
#' @export
RunCollateTpms <- function(collate_tpms, samples, output_dir,
  tpms_file = NA, sample_subdirs = FALSE, orf_fasta = NA) {

  print(paste0("collate_tpms.R: ", collate_tpms))

  cmd_template <- "Rscript --vanilla {collate_tpms} --sample-subdirs={sample_subdirs} --output-dir={output_dir}"
  cmd_template <- paste(cmd_template, paste(samples, collapse = " "))
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
  if (!is.na(tpms_file)) {
    tpms_path <- file.path(output_dir, tpms_file)
  } else {
    tpms_path <- file.path(output_dir, "TPMs_collated.tsv")
  }
  testthat::expect_true(file.exists(tpms_path),
      info = paste(tpms_path, "does not exist"))
  # read_tsv returns spec_tbl_df subclass so cast.
  # Suggested by https://github.com/tidyverse/dplyr/issues/5126
  # Alternative is to subset with no arguments e.g. actual[]
  # https://www.tidyverse.org/blog/2018/12/readr-1-3-1/
  collated_tpms <- tibble::as_tibble(readr::read_tsv(tpms_path,
                                                     comment = "#"))
  return(collated_tpms)
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
  expected_tpms <- tibble::tibble(ORF = orfs)
  for (sample in names(sample_data)) {
    expected_tpms <- add_column(expected_tpms,
                                "{sample}" := sample_data[[sample]])
  }
  return(expected_tpms)
}

#' Save sample-specific TPMs. Sample-specific TPM's data is saved
#' in files `<output_dir>/<sample>_tpms.tsv`.
#'
#' @param orfs List of ORFs (character).
#' @param sample_data List of sample data (lists of integer of length
#' `orfs`),  keyed by sample name (list).
#' @param output_dir Output directory (character).
#' @param shuffle Reorder rows of every second sample file? (logical)
#'
#' @export
SaveTpms <- function(orfs, sample_data, output_dir, shuffle = FALSE) {
  orf_decreasing <- FALSE
  for (sample in names(sample_data)) {
    data <- tibble::tibble(ORF = orfs, tpm = sample_data[[sample]])
    if (shuffle) {
      data <- data[order(data$ORF, decreasing = orf_decreasing), ]
    }
    sample_file <- file.path(output_dir, paste0(sample, "_", "tpms.tsv"))
    readr::write_tsv(data, sample_file, col_names = TRUE)
      orf_decreasing <- ! orf_decreasing
  }
}

testthat::test_that("Default", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    output_dir <- getwd()
    SaveTpms(orfs, sample_data, output_dir)
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms, output_dir = output_dir)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("--sample-subdirs = TRUE", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    for (sample in names(sample_data)) {
      data <- tibble::tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_dir <- file.path(getwd(), sample)
      dir.create(sample_dir)
      sample_file <- file.path(sample_dir, "tpms.tsv")
      readr::write_tsv(data, sample_file, col_names = TRUE)
    }
    output_dir <- getwd()
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms, output_dir = output_dir, sample_subdirs = TRUE)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("One sample", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    output_dir <- getwd()
    SaveTpms(orfs, sample_data, output_dir)
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms, output_dir = output_dir)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("--tpms-file", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    output_dir <- getwd()
    SaveTpms(orfs, sample_data, output_dir)
    tpms_file <- "collated_tpms.tsv"
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms, output_dir = output_dir, tpms_file = tpms_file)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Header-only sample files", {
  samples <- c("WTone", "WTtwo")
  expected_tpms <- tibble::tibble(ORF = character())
  withr::with_tempdir({
    for (sample in samples) {
      data <- tibble::tibble(ORF = character(), tpm = double())
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
      # Expected_tpms column types are `character` as collated TPMs TSV
      # file will have only a header.
      expected_tpms <- add_column(expected_tpms, "{sample}" := character())
    }
    output_dir <- getwd()
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms, output_dir = output_dir)
  })
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Sample-specific ORFs in different orders", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    output_dir <- getwd()
    SaveTpms(orfs, sample_data, output_dir, TRUE)
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms, output_dir = output_dir)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("--orf-fasta", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  fasta <- Biostrings::DNAStringSet(x = replicate(length(orfs), "GATTACCA"))
  names(fasta) <- orfs
  samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    orf_fasta <- "test.fasta"
    Biostrings::writeXStringSet(fasta, orf_fasta, format = "fasta")
    output_dir <- getwd()
    SaveTpms(orfs, sample_data, output_dir, TRUE)
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms, output_dir = output_dir, orf_fasta = orf_fasta)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Missing sample file", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    output_dir <- getwd()
    SaveTpms(orfs, sample_data, output_dir)
    # No file exists for WTtwo,...,WTfour. ORFs will be taken from WTone.
    samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms, output_dir = output_dir)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Missing sample files and --orf-fasta", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  fasta <- Biostrings::DNAStringSet(x = replicate(length(orfs), "GATTACCA"))
  names(fasta) <- orfs
  withr::with_tempdir({
    orf_fasta <- "test.fasta"
    Biostrings::writeXStringSet(fasta, orf_fasta, format = "fasta")
    output_dir <- getwd()
    # No file exists for WTone,...,WTfour. ORFs will be taken from orf_fasta.
    samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
    actual_tpms <- RunCollateTpms(samples = samples, collate_tpms =
      collate_tpms,  output_dir = output_dir, orf_fasta = orf_fasta)
  })
  expected_tpms <- tibble::tibble(ORF = orfs)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})
