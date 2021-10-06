#' testthat tests for `YAL5-single-codons.R`.
#'
#' These tests run `YAL5-single-codons.R`.
#'
#' The tests assumes the following files are in the path:
#'
#' ```
#' rscripts/YAL5-single-codons.R
#' ```
#'
#' To run the tests interactively from within R:
#'
#' ```
#' test_file("rscripts/tests/testthat/test_YAL5_single_codons.R")
#' ```
#'
#' To run the tests from the console:
#'
#' ```
#' Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_YAL5_single_codons.R")'
#' ```
#'
#' These tests assumes the following test data files exist:
#'
#' ```
#' data/Mok-simYAL5/A.h5
#' data/Mok-simYAL5/Scer_YAL_5genes_w_250utrs.gff3
#' data/yeast_standard_asite_disp_length.txt
#' data/yeast_codon_table.tsv
#' data/codons.tsv
#' ```
#'
#' At present, the following behaviours are not tested:
#'
#' * Default values for file-based optional parameters. These are
#'   provided explicitly as default paths are assumed relative to
#'   the riboviz home directory whereas tests are run in a temporary
#'   directory.
#' * `--expand_width`
#' * `--frame`
#' * `--minreadlen`
#' * `--filter_for_frame`
#' * `--snapdisp`
#'
#' @export

suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(withr, quietly = T))

yal5_single_codons <- here::here("rscripts/YAL5-single-codons.R")

context("test_YAL5_single_codons.R")

#' Run `YAL5-single-codons.R`.
#'
#' @param yal5_single_codons `YAL5-single-codons.R` path (character).
#' @param h5_file H5 input file (character).
#' @param dataset Name of dataset of sample being studied (character).
#' @param gff_file GFF feature file of organism being studied (character).
#' @param asite_file Species-specific A-site displacement length file
#' (character).
#' @param annotation_file Codon table file (character).
#' @param feature Feature of interest (feature itself or file of
#' features (character).
#' @param output_dir Output directory (character).
#' @param exit_code Expected exit code (integer).
#'
#' @export
RunYal5SingleCodons <- function(yal5_single_codons, h5_file, dataset,
  gff_file, asite_file, annotation_file, feature, output_dir,
  exit_code = 0) {

  cmd_template <- "Rscript --vanilla {yal5_single_codons} --input={h5_file} --dataset={dataset} --gff={gff_file} --asite_length={asite_file} --annotation={annotation_file} --feature='{feature}' --output={output_dir}"
  cmd <- glue::glue(cmd_template)
  print(cmd)
  actual_exit_code <- system(cmd)
  print(paste0("Exit code: ", actual_exit_code))
  testthat::expect_equal(actual_exit_code, exit_code,
    info = "Unexpected exit code")
}

testthat::test_that("--feature codon", {

  h5_file <- here::here("data/Mok-simYAL5/A.h5")
  dataset <- "Mok-simYAL5"
  gff_file <- here::here("data/Mok-simYAL5/Scer_YAL_5genes_w_250utrs.gff3")
  asite_file <- here::here("data/yeast_standard_asite_disp_length.txt")
  annotation_file <- here::here("data/yeast_codon_table.tsv")
  feature <- "CGA"
  output_dir <- "."

  withr::with_tempdir({
    RunYal5SingleCodons(yal5_single_codons, h5_file, dataset,
      gff_file, asite_file, annotation_file, feature, output_dir)

    output_file <- "Meta_feature_plot_CGA_Mok-simYAL5.pdf"
    output_path <- file.path(".", output_file)
    testthat::expect_true(file.exists(output_path),
         info = paste(output_path, "does not exist"))
  })
})

testthat::test_that("--feature file", {

  h5_file <- here::here("data/Mok-simYAL5/A.h5")
  dataset <- "Mok-simYAL5"
  gff_file <- here::here("data/Mok-simYAL5/Scer_YAL_5genes_w_250utrs.gff3")
  asite_file <- here::here("data/yeast_standard_asite_disp_length.txt")
  annotation_file <- here::here("data/yeast_codon_table.tsv")
  feature <- here::here("data/codons.tsv")
  output_dir <- "."

  withr::with_tempdir({
    RunYal5SingleCodons(yal5_single_codons, h5_file, dataset,
      gff_file, asite_file, annotation_file, feature, output_dir)

    output_file <- "Feature_Relative_use_Mok-simYAL5.tsv"
    output_path <- file.path(".", output_file)
    testthat::expect_true(file.exists(output_path),
         info = paste(output_path, "does not exist"))
  })
})
