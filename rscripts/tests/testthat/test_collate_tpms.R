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
#' @param samples List of sample files (where `names` attribute of
#' `samples` are the sample names) (list).
#' @param tpms_file Name of collated TPMs file (character).
#' @param orf_fasta ORF FASTA file that was aligned to and from which
#' ORF names are to be retrieved (character).
#' @param expected_exit_code Expected exit code (integer).
#' @return Tibble with collated TPMs (`NULL` if `expected_exit_code`
#' is not 0) (tbl_df).
#'
#' @export
RunCollateTpms <- function(collate_tpms, samples, tpms_file = NA,
  orf_fasta = NA, expected_exit_code = 0) {

  print(paste0("collate_tpms.R: ", collate_tpms))
  cmd_template <- "Rscript --vanilla {collate_tpms}"
  if (!is.na(tpms_file)) {
    cmd_template <- paste(cmd_template, "--tpms-file={tpms_file}")
  }
  if (!is.na(orf_fasta)) {
    cmd_template <- paste(cmd_template, "--orf-fasta={orf_fasta}")
  }
  cmd_template <- paste(cmd_template,
    paste(as.vector(rbind(names(samples), samples)),
    collapse = " "))
  cmd <- glue::glue(cmd_template)
  print(cmd)

  exit_code <- system(cmd)
  print(paste0("'collate_tpms.R' exit code: ", exit_code))
  testthat::expect_equal(expected_exit_code, exit_code,
    info = "Unexpected exit code from 'collate_tpms.R'")
  if (!is.na(tpms_file)) {
    tpms_path <- file.path(tpms_file)
  } else {
    tpms_path <- file.path("TPMs_collated.tsv")
  }
  if ((expected_exit_code) != 0) {
    return(NULL)
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
#' @return List of sample files (where `names` attribute of
#' `samples` are the sample names) (list).
#'
#' @export
SaveTpms <- function(orfs, sample_data, output_dir, shuffle = FALSE) {
  orf_decreasing <- FALSE
  sample_files <- character()
  for (sample in names(sample_data)) {
    data <- tibble::tibble(ORF = orfs, tpm = sample_data[[sample]])
    if (shuffle) {
      data <- data[order(data$ORF, decreasing = orf_decreasing), ]
    }
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
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
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
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files)
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
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files)
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
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files)
  })
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Sample files with ORFs in different orders", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo", "WTthree", "WTfour")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd(), TRUE)
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Sample files with one non-existent sample file", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    for (sample in c("WTtwo", "WTthree", "WTfour")) {
      sample_files[[sample]] <- "nosuchfile.tsv"
    }
    # ORFs will be taken from WTone.
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("--tpms-file", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    tpms_file <- "collated_tpms.tsv"
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files, tpms_file)
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
    sample_files <- SaveTpms(orfs, sample_data, getwd(), TRUE)
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files,
      orf_fasta = orf_fasta)
  })
  expected_tpms <- GetCollatedTpms(orfs, sample_data)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Non-existent --orf-fasta file raises error", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files,
      orf_fasta = "nosuchfile.fa", expected_exit_code = 1)
  })
})

testthat::test_that("Non-existent samples file for ORFs raises error", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    sample_files[["WTone"]] <- "nosuchfile.tsv"
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files,
      expected_exit_code = 1)
  })
})

testthat::test_that("Non-existent samples files and --orf-fasta", {
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
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files,
      orf_fasta = orf_fasta)
  })
  expected_tpms <- tibble::tibble(ORF = orfs)
  testthat::expect_equal(expected_tpms, actual_tpms, info = "TPMs differ")
})

testthat::test_that("Missing sample file for sample raises error", {
  orfs <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  samples <- c("WTone", "WTtwo")
  sample_data <- GetSampleTpms(orfs, samples)
  withr::with_tempdir({
    sample_files <- SaveTpms(orfs, sample_data, getwd())
    sample_files[["WTtwo"]] <- ""
    actual_tpms <- RunCollateTpms(collate_tpms, sample_files,
      expected_exit_code = 1)
  })
})
