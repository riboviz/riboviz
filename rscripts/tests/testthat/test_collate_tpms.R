#' testthat tests for `collate_tpms.R`.
#'
#' These tests run `bam_to_h5.R` using a sample-specific TPMs files
#' then validate the collated TPMs file output.
#'
#' The tests assumes the following files are in the path:
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
#' These tests assumes the following test data files exist:
#'
#' ```
#' data/
#'   tpms/
#'     TPMs_collated.tsv
#'     nextflow/
#'       WT3AT_tpms.tsv
#'       WTnone_tpms.tsv
#'     python/
#'       WT3AT/
#'         tpms.tsv
#'       WTnone/
#'         tpms.tsv
#' ```
#'
#' @export

suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(withr, quietly = T))

print(paste0("here: ", here()))
collate_tpms <- here::here(file.path("rscripts", "collate_tpms.R"))

context("test_collate_tpms.R")

#' Run `collate_tpms.R` on sample-specific TSV files.
#'
#' @param collate_tpms `collate_tpms.R` path (character).
#' @param output_dir Output directory in which to find sample-specific
#' TPMs files and into which to write collated TPMs file (character).
#' @param tpms_file Name of collated TPMs file, relative to
#' `output_dir` (character).
#' @param sample_subdirs Are samples-specific TPMs files in
#' sample-specific sub-directories, in files
#' `<output_dir>/<sample>/tpms.tsv`? If not then it is assumed
#' they are in files `<output_dir>/<sample>_tpms.tsv` (logical).
#' @param orf_fasta FASTA file from which ORF names are to be
#' retrieved (character).
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

#' Compare TSV files for equality.
#'
#' @param expected_tsv TSV with expected data (character).
#' @param actual_tsv TSV with actual data (character).
#;
#' @export
CompareTsv <- function(expected_tsv, actual_tsv) {
  expected <- read.table(file = expected_tsv)
  actual <- read.table(file = actual_tsv)
  testthat::expect_equal(expected, actual, info = "Data differs'")
}

testthat::test_that("collate_tpms.R: Python workflow", {
  output_dir <- here::here(file.path("data", "tpms", "python"))
  sample_subdirs <- TRUE
  tpms_file <- NA # "TPMs_collated.tsv"
  orf_fasta <- NA
  samples <- c("WTnone", "WT3AT")
  expected_tpms_file <- here::here(
    file.path("data", "tpms", "TPMs_collated.tsv"))
  withr::with_tempdir({
    RunCollateTpms(collate_tpms, output_dir, tpms_file,
      sample_subdirs, orf_fasta, samples)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    CompareTsv(expected_tpms_file, actual_tpms_file)
  })
})

testthat::test_that("collate_tpms.R: Nextflow workflow", {
  output_dir <- here::here(file.path("data", "tpms", "nextflow"))
  sample_subdirs <- FALSE
  tpms_file <- NA # "TPMs_collated.tsv"
  orf_fasta <- NA
  samples <- c("WTnone", "WT3AT")
  expected_tpms_file <- here::here(
    file.path("data", "tpms", "TPMs_collated.tsv"))
  withr::with_tempdir({
    RunCollateTpms(collate_tpms, output_dir, tpms_file,
      sample_subdirs, orf_fasta, samples)
    actual_tpms_file <- file.path(output_dir, "TPMs_collated.tsv")
    CompareTsv(expected_tpms_file, actual_tpms_file)
  })
})
