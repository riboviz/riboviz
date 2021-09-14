#' testthat tests for `YAL5-codon-pairs.R`.
#' 
#' These tests run `YAL5-codon-pairs.R` using #### file 
#' then validates #### based upon its expected qualities 
#' given those of the input #### files and the `YAL5-codon-pairs.R`
#' command-line parameters.
#' 
#' The tests assume the following files are in the path:
#' 
#' ```
#' rscripts/YAL5-codon-pairs.R
#' ```
#' 
#' To run the tests interactively from within R:
#' 
#' ```
#' test_file("rscripts/tests/testthat/test_YAL5-codon-pairs.R")
#' ```
#' 
#' To run the tests from the console:
#' 
#' ```
#' Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_YAL5-codon-pairs.R")'
#' ```
#' 
#' These tests assumes the following test data files exist:
#' 
#' `data/Mok-tinysim-gffsam/A.sam`
#' `data/Mok-tinysim-gffsam/tiny_2genes_20utrs.gff3`
#' 
#' To use these tests with other data files see comments on 
#' test_that("Run YAL5-codon-pairs.R and validate #### " ...).
#' 
#' At present, the following behaviours are not tested:
#' 
#'   ####
#' 
#' @export

suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(plotly, quietly = T))
suppressMessages(library(purrr, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(optparse, quietly = T))
suppressMessages(library(tidyverse, quietly = T))
suppressMessages(library(stringr, quietly = T))

# (source(here::here("rscripts", "read_count_functions.R")
print(paste0("here: ", here()))
# codon_pairs <- here::here(file.path("rscripts", "YAL5-codon-pairs.R"))

context("test_YAL5-codon-pairs.R")

#' Run `test_YAL5-codon-pairs.R` on ####.
#' 
#' 