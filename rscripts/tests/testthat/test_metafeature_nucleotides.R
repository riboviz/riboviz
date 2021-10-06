#' testthat tests for `metafeature_nucleotides.R`.
#'
#' These tests run `metafeature_nucleotides.R`.
#'
#' The tests assumes the following files are in the path:
#'
#' ```
#' rscripts/metafeature_nucleotides.R
#' ```
#'
#' To run the tests interactively from within R:
#'
#' ```
#' test_file("rscripts/tests/testthat/test_metafeature_nucleotides.R")
#' ```
#'
#' To run the tests from the console:
#'
#' ```
#' Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_metafeature_nucleotides.R")'
#' ```
#'
#' These tests assumes the following test data files exist:
#'
#' ```
#' data/Mok-tinysim/A.h5
#' data/Mok-tinysim/tiny_2genes_20utrs.gff3
#' data/Mok-tinysim/tiny_2genes_20utrs.fa
#' data/yeast_standard_asite_disp_length.txt
#' data/feature_pos.tsv
#' ```
#'
#' At present, the following behaviours are not tested:
#'
#' * Default values for file-based optional parameters. These are
#'   provided explicitly as default paths are assumed relative to
#'   the riboviz home directory whereas tests are run in a temporary
#'   directory.
#' * `--expand_width`
#' * `--minreadlen`
#'
#' @export

suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(withr, quietly = T))

metafeature_nucleotides <- here::here("rscripts/metafeature_nucleotides.R")

context("test_metafeature_nucleotides.R")

#' Run `metafeature_nucleotides.R`.
#'
#' @param metafeature_nucleotides `metafeature_nucleotides.R` path (character).
#' @param h5_file H5 input file (character).
#' @param dataset Name of dataset of sample being studied (character).
#' @param gff_file GFF feature file of species being studied (character).
#' @param fasta_file FASTA file of species being studied (character).
#' @param asite_file Species-specific A-site displacement length file
#' (character).
#' @param feature_pos_file Genes and positions to normalize over (character).
#' @param output_dir Output directory (character).
#' @param expand_width Desired range either side of feature (integer).
#' @param exit_code Expected exit code (integer).
#'
#' @export
RunMetafeatureNucleotides <- function(metafeature_nucleotides,
  h5_file, dataset, gff_file, fasta_file, asite_file, feature_pos_file,
  output_dir, expand_width, exit_code = 0) {

  cmd_template <- "Rscript --vanilla {metafeature_nucleotides} --input={h5_file} --dataset={dataset} --gff={gff_file} --fasta={fasta_file} --asite_length={asite_file} --feature_pos='{feature_pos_file}' --output={output_dir} --expand_width={expand_width}"
  cmd <- glue::glue(cmd_template)
  print(cmd)
  actual_exit_code <- system(cmd)
  print(paste0("Exit code: ", actual_exit_code))
  testthat::expect_equal(actual_exit_code, exit_code,
    info = "Unexpected exit code")
}

testthat::test_that("Default", {

  h5_file <- here::here("data/Mok-tinysim/A.h5")
  dataset <- "Mok-tinysim"
  gff_file <- here::here("data/Mok-tinysim/tiny_2genes_20utrs.gff3")
  fasta_file <- here::here("data/Mok-tinysim/tiny_2genes_20utrs.fa")
  asite_file <- here::here("data/yeast_standard_asite_disp_length.txt")
  feature_pos_file <- here::here("data/feature_pos.tsv")
  output_dir <- "."
  expand_width <- 2

  withr::with_tempdir({
    RunMetafeatureNucleotides(metafeature_nucleotides, h5_file, dataset,
      gff_file, fasta_file, asite_file, feature_pos_file, output_dir,
      expand_width)

    output_file <- "Meta_feature_plot_positons_of_interest_Mok-tinysim.pdf"
    output_path <- file.path(".", output_file)
    testthat::expect_true(file.exists(output_path),
         info = paste(output_path, "does not exist"))
  })
})
