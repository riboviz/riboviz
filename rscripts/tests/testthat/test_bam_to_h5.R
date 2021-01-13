library(glue)
library(here)
library(withr)

print(here())
bam_to_h5 <- here::here("rscripts/bam_to_h5.R")
print(bam_to_h5)
bam_file = here::here("test-bam-to-h5/WTnone/WTnone.bam")
print(bam_file)
gff_file = here::here("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
print(gff_file)
h5_file = here::here("WTnone.h5")
print(h5_file)
expected_h5_file = here::here("test-bam-to-h5/WTnone/WTnone.h5")
print(expected_h5_file)

context("test_bam_to_h5.R")

delete_file <- function(file_name) {
  print("Invoking delete_file fixture...")
  if (file.exists(file_name)) {
    file.remove(file_name)
  }
}

test_that("Run bam_to_h5.R and validate H5 file", {
  withr::defer(delete_file(h5_file))
  bam_to_h5_cmd_template <- "Rscript --vanilla {bam_to_h5} --num-processes=1 --min-read-length=10 --max-read-length=50 --buffer=250 --primary-id=Name --secondary-id=NULL --dataset=vignette --bam-file={bam_file} --hd-file={h5_file} --orf-gff-file={gff_file} --is-riboviz-gff=true --stop-in-cds=false"
  print(bam_to_h5_cmd_template)
  cmd <- glue(bam_to_h5_cmd_template)
  print(cmd)
  exit_code <- system(cmd)
  print(glue("bam_to_h5.R exit code: {exit_code}"))
  expect_equal(exit_code, 0, info = "Unexpected exit code from bam_to_h5.R")
  h5diff_cmd_template <- "h5diff -q {h5_file} {expected_h5_file}"
  cmd <- glue(h5diff_cmd_template)
  print(cmd)		
  exit_code <- system(cmd)
  print(glue("h5diff exit code: {exit_code}"))
  expect_equal(exit_code, 0, info = "Unexpected exit code from h5diff")
})
