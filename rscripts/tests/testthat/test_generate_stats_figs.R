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

context("Set-up & imports")

test_that("passed-in options all present", {
  
  # checking options / opt set-up is fine, we expect 22 items 
  #  (21 riboviz options passed to generate_stats_figs.R + help = 22)
  expect_equal(length(opt), 22)

})

test_that("all required files exist, are readable", {

  expect_equal(file.exists(hd_file), as.logical("TRUE"))
  
  # str(h5ls(here(opt$hd_file), recursive = 1))
  # 'data.frame':	68 obs. of  5 variables:
  #   $ group : chr  "/" "/" "/" "/" ...
  # $ name  : chr  "YAL001C" "YAL002W" "YAL003W" "YAL005C" ...
  # $ otype : chr  "H5I_GROUP" "H5I_GROUP" "H5I_GROUP" "H5I_GROUP" ...
  # $ dclass: chr  "" "" "" "" ...
  # $ dim   : chr  "" "" "" "" ...
  
})