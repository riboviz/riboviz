#' testthat tests for metafeature_nucleotides.R
#'
#' These tests run metafeature_nucleotides.R using #### file then
#' validates ##### based upon its expected qualities
#' given those of the input ##### files and the metafeature_nucleotides.R
#' command-line parameters.
#'
#' The tests assumes the following files are in the path:
#'
#' * `rscripts/metafeature_nucleotides.R`
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
#' * `data/Mok-tinysim-gffsam/A.sam`
#' * `data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3`
#'
#' To use thse tests with other data files see comments on
#' test_that("Run metafeature_nucleotides.R and validate ##### "...).
#'
#' At present, the following behaviours are not tested:
#'  #######
#'
#' @export

suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(withr, quietly = T))
suppressMessages(library(GenomicAlignments, quietly = T))
suppressMessages(library(Rsamtools, quietly = T))

#source(here::here("rscripts", "read_count_functions.R"))
print(paste0("here: ", here()))
#metafeature_nucleotides <- here::here("rscripts/metafeature_nucleotides.R")

context("test_metafeature_nucleotides.R")

