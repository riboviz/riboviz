#' Collate TPMs
#'
#' Collate sample-specific TPMs files into a single file with data
#' for all the samples.
#'
#' Each sample-specific TPMs file is assumed to be a tab-separated
#' values (TSV) file with `ORF` and `tpm` columns (other columns
#' are ignored.
#'
#' The collated TPMs file is a TSV file with an `ORF` column and
#' sample-specific columns, named by sample name, with the `tpm`
#' values from the sample-specific files.
#'
#' This script accepts the following command-line parameters:
#'
#' | Parameter | Description |
#' | --------- | ----------- |
#' | `--tpms-file` | Name of collated TPMs file |
#' | `--orf-fasta` | ORF FASTA file that was aligned to and from which ORF names are to be retrieved |
#' | `--sort-orfs` | Sort ORF list lexicographically (default `FALSE`) |
#' | `--digits` | Number of decimal places to be used in output (default 1) |
#' | `<sample> <sample_tpms>` | Space-delimited list of one or more pairs of sample names and sample-specific TPMs files |
#'
#' @export

suppressMessages(library(getopt, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(optparse, quietly = T))

# Load local dependencies.
if (interactive()) {
  # Use hard-coded script name and assume script is in "rscripts"
  # directory. This assumes that interactive R is being run within
  # the parent of rscripts/ but imposes no other constraints on
  # where rscripts/ or its parents are located.
  self <- "collate_tpms.R"
  path_to_self <- here("rscripts", self)
  source(here::here("rscripts", "provenance.R"))
  source(here::here("rscripts", "collate_tpms_functions.R"))
} else {
  # Deduce file name and path using reflection as before.
  self <- getopt::get_Rscript_filename()
  path_to_self <- self
  source(file.path(dirname(self), "provenance.R"))
  source(file.path(dirname(self), "collate_tpms_functions.R"))
}

option_list <- list(
  make_option("--tpms-file", type = "character",
    default = "TPMs_all_CDS_all_samples.tsv",
    help = "Name of collated TPMs file"),
  make_option("--orf-fasta", type = "character", default = NA,
    help = "ORF FASTA file that was aligned to and from which ORF names are to be retrieved"),
  make_option("--sort-orfs", type = "logical", default = FALSE,
    help = "Sort ORF list lexicographically"),
  make_option("--digits", type = "integer", default = 1,
    help = "Number of decimal places to be used in output"))


print_provenance(get_Rscript_filename())

opt <- parse_args(OptionParser(option_list = option_list),
                  positional_arguments = TRUE,
                  convert_hyphens_to_underscores = TRUE)
print("collate_tpms.R running with parameters:")
print(opt)

if (length(opt$args) == 0) {
  stop("No <sample_name> <sample_file> list provided")
}
if ((length(opt$args) %% 2) != 0) {
  stop("Invalid <sample_name> <sample_file> list provided")
}
sample_names_files <- split(opt$args, ceiling(seq_along(opt$args) %% 2))
sample_files <- sample_names_files[[1]]
names(sample_files) <- sample_names_files[[2]]

CollateTpms(opt$options$tpms_file,
            opt$options$orf_fasta,
            sample_files,
            sort_orfs = opt$options$sort_orfs,
            digits = opt$options$digits)

print("collate_tpms.R done")
