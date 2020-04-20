# test_generate_stats_figs.R
# tests for code in generate_stats_figs.R using {testthat}
# https://github.com/riboviz/riboviz/

# libraries

# Install package from CRAN only if not installed, and load the library
#if (!require(testthat)) install.packages('testthat') # TODO FLIC be aware of this addition
library(testthat)
library(here)

# So, if you have a file with code, for example in R/code.R:
#   
#   add_one <- function(x) x + 1
#
# your test file should be tests/testthat/test_code.R:
#   
# source("../../R/code.R") # This is only needed if your project is not a package
# 
# test_that("Add one to 99", {
#   expect_equal(add_one(99), 100)
# })

# source generate_stats_figs.R:
source(here::here("rscripts", "generate_stats_figs.R"))

# checking options / opt set-up is fine, we expect 22 items 
#  (21 riboviz options passed to generate_stats_figs.R + help = 22)
expect_equal(length(opt), 22)

