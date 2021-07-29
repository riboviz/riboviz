#' testthat tests for `collate_tpms_functions.R`.
#'
#' These tests run `CollateTpms` using a sample-specific TPMs files
#' then validate the collated TPMs file output.
#'
#' The tests assume the following files are in the path:
#'
#' ```
#' rscripts/collate_tpms_functions.R
#' ```
#'
#' To run the tests interactively from within R:
#'
#' ```
#' test_file("rscripts/tests/testthat/test_collate_tpms_functions.R")
#' ```
#'
#' To run the tests from the console:
#'
#' ```
#' Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_collate_tpms_functions.R")'
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

source(here::here("rscripts", "collate_tpms_functions.R"))

context("test_collate_tpms_functions.R")

#' Run `CollateTpms` on sample-specific TSV files.
#'
#' @param samples List of sample files (where `names` attribute of
#' `samples` are the sample names) (list).
#' @param tpms_file Name of collated TPMs file (character).
#' @param orf_fasta ORF FASTA file that was aligned to and from which
#' ORF names are to be retrieved (character).
#' @param sort_orfs sort ORF list lexicographically (logical)
#' @param digits Number of decimal places to be used in table in file
#' output (integer).
#' @return Tibble with collated TPMs (tbl_df).
#'
#' @export
RunCollateTpms <- function(samples, tpms_file = "tpms.tsv",
  orf_fasta = NA, sort_orfs = FALSE, digits = 1) {

  if (!is.na(tpms_file)) {
    tpms_path <- file.path(tpms_file)
  } else {
    tpms_path <- file.path("TPMs_collated.tsv")
  }
  CollateTpms(tpms_file, orf_fasta, samples, sort_orfs, digits)
  testthat::expect_true(file.exists(tpms_path),
      info = paste(tpms_path, "does not exist"))
  # read_tsv returns spec_tbl_df subclass so cast.
  # Suggested by https://github.com/tidyverse/dplyr/issues/5126
  # Alternative is to subset with no arguments e.g. actual[]
  # https://www.tidyverse.org/blog/2018/12/readr-1-3-1/
  collated_tpms <- tibble::as_tibble(readr::read_tsv(tpms_path,
                                                     comment = "#"))
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
#' @return List of sample files (where `names` attribute of
#' `samples` are the sample names) (list).
#'
#' @export
SaveTpms <- function(orfs, sample_data, output_dir) {
  orf_decreasing <- FALSE
  sample_files <- character()
  for (sample in names(sample_data)) {
    data <- tibble::tibble(ORF = orfs, tpm = sample_data[[sample]])
    sample_file <- file.path(output_dir, paste0(sample, "_", "tpms.tsv"))
    readr::write_tsv(data, sample_file, col_names = TRUE)
      orf_decreasing <- ! orf_decreasing
    sample_files <- append(sample_files, sample_file)
  }
  names(sample_files) <- names(sample_data)
  return(sample_files)
}

testthat::test_that("Defaults", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    actual_tpms <- RunCollateTpms(sample_files)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("digits=0", {
  orfs <- c("YAL004C", "YAL003W", "YAL002W", "YAL001W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    actual_tpms <- RunCollateTpms(sample_files, digits = 0)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  expected_tpms <- dplyr::mutate_if(expected_tpms, is.numeric, round,
                                    digits = 0)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Sample files in directories", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- character()
    for (sample in names(sample_data)) {
      data <- tibble::tibble(ORF = orfs, tpm = sample_data[[sample]])
      sample_dir <- file.path(getwd(), sample)
      dir.create(sample_dir)
      sample_file <- file.path(sample_dir, "tpms.tsv")
      readr::write_tsv(data, sample_file, col_names = TRUE)
      sample_files <- append(sample_files, sample_file)
    }
    names(sample_files) <- samples
    actual_tpms <- RunCollateTpms(sample_files)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("One sample only", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    actual_tpms <- RunCollateTpms(sample_files)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Header-only sample files", {
  samples <- c("WTone", "WTtwo")
  expected_tpms <- tibble::tibble(ORF = character())
  withr::with_tempdir({
    sample_files <- character()
    for (sample in samples) {
      data <- tibble::tibble(ORF = character(), tpm = double())
      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(data, sample_file, col_names = TRUE)
      # Expected_tpms column types are `character` as collated TPMs TSV
      # file will have only a header.
      expected_tpms <- add_column(expected_tpms, "{sample}" := character())
      sample_files <- append(sample_files, sample_file)
    }
    names(sample_files) <- samples
    actual_tpms <- RunCollateTpms(sample_files)
  })
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("sort_orfs=TRUE", {
  orfs <- c("YAL004C", "YAL003W", "YAL002W", "YAL001W")
  samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    actual_tpms <- RunCollateTpms(sample_files, sort_orfs = TRUE)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  expected_tpms <- dplyr::arrange(expected_tpms, ORF)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("sort_orfs=FALSE", {
  orfs <- c("YAL004C", "YAL003W", "YAL002W", "YAL001W")
  samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    actual_tpms <- RunCollateTpms(sample_files)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Sample files with one non-existent sample file", {
  # FIXME  
  skip("FIXME")
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    for (sample in c("WTtwo", "WTthree", "WTfour")) {
      sample_files[[sample]] <- "nosuchfile.tsv"
    }
    # ORFs will be taken from WTone.
    actual_tpms <- RunCollateTpms(sample_files)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("tpms_file", {
  orfs <- c("YAL004C", "YAL003W", "YAL002W", "YAL001W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    tpms_file <- "collated_tpms.tsv"
    actual_tpms <- RunCollateTpms(sample_files, tpms_file)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("orf_fasta (sort_orfs=TRUE)", {
  orfs <- c("YAL004C", "YAL003W", "YAL002W", "YAL001W")
  fasta <- Biostrings::DNAStringSet(x = replicate(length(orfs), "GATTACCA"))
  names(fasta) <- orfs
  samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    orf_fasta <- "test.fasta"
    Biostrings::writeXStringSet(fasta, orf_fasta, format = "fasta")
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    actual_tpms <- RunCollateTpms(sample_files,
      orf_fasta = orf_fasta, sort_orfs = TRUE)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  expected_tpms <- dplyr::arrange(expected_tpms, ORF)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("orf_fasta (sort_orfs=FALSE)", {
  orfs <- c("YAL004C", "YAL003W", "YAL002W", "YAL001W")
  fasta <- Biostrings::DNAStringSet(x = replicate(length(orfs), "GATTACCA"))
  names(fasta) <- orfs
  samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    orf_fasta <- "test.fasta"
    Biostrings::writeXStringSet(fasta, orf_fasta, format = "fasta")
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    actual_tpms <- RunCollateTpms(sample_files,
      orf_fasta = orf_fasta)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Non-existent orf_fasta file raises error", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    expect_error(
      actual_tpms <- RunCollateTpms(sample_files,
        orf_fasta = "nosuchfile.fa")
    )
  })
})

testthat::test_that("Non-existent samples file for ORFs raises error", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    sample_files[["WTone"]] <- "nosuchfile.tsv"
    expect_error(
      actual_tpms <- RunCollateTpms(sample_files)
    )
  })
})

testthat::test_that("Non-existent samples files and orf_fasta", {
  # FIXME  
  skip("FIXME")
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  fasta <- Biostrings::DNAStringSet(x = replicate(length(orfs), "GATTACCA"))
  names(fasta) <- orfs
  withr::with_tempdir({
    orf_fasta <- "test.fasta"
    Biostrings::writeXStringSet(fasta, orf_fasta, format = "fasta")
    sample_files <- character()
    for (sample in c("WTone", "WTtwo", "WTthree", "WTfour")) {
      sample_files[[sample]] <- "nosuchfile.tsv"
    }
    actual_tpms <- RunCollateTpms(sample_files, orf_fasta = orf_fasta)
  })
  expected_tpms <- tibble::tibble(ORF = orfs)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Missing sample file for sample raises error", {
  # FIXME  
  skip("FIXME")
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    sample_files[["WTtwo"]] <- ""
    expect_error(
      actual_tpms <- RunCollateTpms(sample_files)
    )
  })
})
