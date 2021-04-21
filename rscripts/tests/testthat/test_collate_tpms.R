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

testthat::test_that("collate_tpms.R: Python workflow", {
  sample_subdirs <- TRUE
  tpms_file <- NA
  orf_fasta <- NA
  orfs <- c("YAL001C", "YAL002W", "YAL003W")
  samples <- c("WTnone", "WT3AT")
  expected <- tibble(ORF = orfs)
  offset <- 1
  withr::with_tempdir({
   for (sample in samples) {
      new_offset <- offset + length(orfs)
      data <- offset:(new_offset - 1)
      offset <- new_offset
      sample_tpms <- tibble(ORF = orfs, tpm = data)

      sample_dir <- file.path(getwd(), sample)
      dir.create(sample_dir)
      sample_file <- file.path(sample_dir, "tpms.tsv")
      readr::write_tsv(sample_tpms, sample_file, col_names = TRUE)

      expected <- add_column(expected, "{sample}" := data)
    }
    output_dir <- getwd()
    RunCollateTpms(collate_tpms, output_dir, tpms_file,
      sample_subdirs, orf_fasta, samples)

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

testthat::test_that("collate_tpms.R: Nextflow workflow", {
  sample_subdirs <- FALSE
  tpms_file <- NA
  orf_fasta <- NA
  orfs <- c("YAL001C", "YAL002W", "YAL003W")
  samples <- c("WTnone", "WT3AT")
  expected <- tibble(ORF = orfs)
  offset <- 1
  withr::with_tempdir({
   for (sample in samples) {
      new_offset <- offset + length(orfs)
      data <- offset:(new_offset - 1)
      offset <- new_offset
      sample_tpms <- tibble(ORF = orfs, tpm = data)

      sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
      readr::write_tsv(sample_tpms, sample_file, col_names = TRUE)

      expected <- add_column(expected, "{sample}" := data)
    }
    output_dir <- getwd()
    RunCollateTpms(collate_tpms, output_dir, tpms_file,
      sample_subdirs, orf_fasta, samples)

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

testthat::test_that("collate_tpms.R: Single column", {
  sample_subdirs <- FALSE
  tpms_file <- NA
  orf_fasta <- NA

  sample <- "WTnone"
  orfs <- c("YAL001C", "YAL002W", "YAL003W")
  data <- seq_len(length(orfs))
  sample_tpms <- tibble(ORF = orfs, tpm = data)
  expected <- tibble(ORF = orfs, "{sample}" := data)
  withr::with_tempdir({
    sample_file <- file.path(getwd(), paste0(sample, "_", "tpms.tsv"))
    readr::write_tsv(sample_tpms, sample_file, col_names = TRUE)

    output_dir <- getwd()
    RunCollateTpms(collate_tpms, output_dir, tpms_file,
      sample_subdirs, orf_fasta, c(sample))

    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    actual <- readr::read_tsv(actual_tpms_file, comment = "#")
  })
  actual <- tibble::as_tibble(actual)
  testthat::expect_equal(expected, actual, info = "Data differs'")
})
