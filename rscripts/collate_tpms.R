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
#' | `--output-dir ` | Output directory in which to find sample-specific TPMs files and into which to write the collated TPMs file |
#' | `--tpms-file` | Name of collated TPMs file, relative to `output-dir` |
#' | `--sample-subdirs` | Are sample-specific TPMs files in sample-specific sub-directories, in files `<output_dir>/<sample>/tpms.tsv`? If not then it is assumed they are in files `<output_dir>/<sample>_tpms.tsv` |
#' | `--orf-fasta` | ORF FASTA file that was aligned to and from which ORF names are to be retrieved |
#' | `<sample>` | Space-delimited list of one or more sample names |
#'
#' @export

suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(getopt, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(optparse, quietly = T))
suppressMessages(library(purrr, quietly = T))
suppressMessages(library(readr, quietly = T))
suppressMessages(library(tibble, quietly = T))
suppressMessages(library(tidyr, quietly = T))

# Load local dependencies.
if (interactive()) {
  # Use hard-coded script name and assume script is in "rscripts"
  # directory. This assumes that interactive R is being run within
  # the parent of rscripts/ but imposes no other constraints on
  # where rscripts/ or its parents are located.
  self <- "collate_tpms.R"
  path_to_self <- here("rscripts", self)
  source(here::here("rscripts", "provenance.R"))
} else {
  # Deduce file name and path using reflection as before.
  self <- getopt::get_Rscript_filename()
  path_to_self <- self
  source(file.path(dirname(self), "provenance.R"))
}

#' Load TPMs, sort by `ORF` column, check that gene names in `ORF`
#' column is consistent with given `orfs`, and, if so, return `tpm`
#' column.
#'
#' Warnings are printed if `tpms_file` does not exist or if the `ORF`
#' column is inconsistent with `orfs`.
#'
#' @param tpms_file TPMs file (character).
#' @param orfs Sorted list of ORFs (character).
#' @return List of TPMs, consistent with order of `orfs` (double).
#'
#' @export
LoadTpms <- function(tpms_file, orfs) {
  print(paste0("Loading TPMs from: ", tpms_file))
  if (!file.exists(tpms_file)) {
    warning(paste(tpms_file, "does not exist, returning empty list"))
    return(NULL)
  }
  features <- readr::read_tsv(tpms_file, comment = "#")
  features <- features[order(features$ORF), ]
  if (isFALSE(all.equal(features$ORF, orfs))) {
    warning(paste("Inconsistent ORF names in ", tpms_file))
  }
  return(features$tpm)
}

#' Get name of sample-specific TPMs file.
#'
#' @param samples_dir Directory in which to find sample-specific TPMs
#' files (character).
#' @param sample Sample ID, used to access sample-specific TPMs files
#' (character).
#' @param tpms_suffix Suffix for sample-specific TPMs files (character).
#' @param sample_subdirs Are sample-specific TPMs files in
#' sample-specific sub-directories, in files
#' `<samples_dir>/<sample>/<tpms_suffix>`? If not then it is assumed
#' they are in files `<samples_dir>/<sample>_<tpms_suffix>` (logical).
#' @return `<samples_dir>/<sample>/<tpms_suffix>` (if `sample_subdirs`
#' is true) else `<samples_dir>/<sample>_<tpms_suffix>` (character).
#'
#' @export
GetTpmsFileName <- function(
  samples_dir, sample, tpms_suffix, sample_subdirs) {
  if (sample_subdirs) {
    file_name <- file.path(samples_dir, sample, tpms_suffix)
  } else {
    file_name <- file.path(samples_dir, paste0(sample, "_", tpms_suffix))
  }
  return(file_name)
}

#' Get TPMs for a specific sample.
#'
#' @param sample Sample ID, used to access sample-specific TPMs files
#' (character).
#' @param samples_dir Directory in which to find sample-specific TPMs
#' files (character).
#' @param tpms_suffix Suffix for sample-specific TPMs files (character).
#' @param sample_subdirs Are sample-specific TPMs files in
#' sample-specific sub-directories, in files
#' `<samples_dir>/<sample>/<tpms_suffix>`? If not then it is assumed
#' they are in files `<samples_dir>/<sample>_<tpms_suffix>` (logical).
#' @param orfs Sorted list of ORFs (character).
#' @return List of TPMs, consistent with order of `orfs` (double).
#'
#' @export
GetTpms <- function(sample, samples_dir, tpms_suffix, sample_subdirs, orfs) {
  tpms <- LoadTpms(
    GetTpmsFileName(samples_dir, sample, tpms_suffix, sample_subdirs), orfs)
  return(tpms)
}

#' Collate TPMs from sample-specific files and return collated TPMs.
#'
#' If `orf_fasta` is provided then gene names are loaded from this
#' file, else they are loaded from the first sample-specific TPMs
#' file's `ORF` column. ORFs are sorted after loading.
#'
#' @param samples_dir Directory in which to find sample-specific TPMs
#' files (character).
#' @param sample_subdirs Are sample-specific TPMs files in
#' sample-specific sub-directories, in files
#' `<samples_dir>/<sample>/tpms.tsv`? If not then it is assumed
#' they are in files `<samples_dir>/<sample>_tpms.tsv` (logical).
#' @param orf_fasta ORF FASTA file that was aligned to and from which
#' ORF names are to be retrieved (character).
#' @param samples List of sample names (character).
#' @param tpms_suffix Suffix for sample-specific TPMs files (character).
#' @return Collated TPMs with an `ORF` column and sample-specific
#' columns, named by sample name, with the TPMs values for each sample
#' (data.frame).
#'
#' @export
MakeTpmTable <- function(samples_dir, sample_subdirs, orf_fasta,
                         samples, tpms_suffix="tpms.tsv") {
  if (is.na(orf_fasta)) {
    tpms_file <- GetTpmsFileName(samples_dir,
                                 samples[1],
                                 tpms_suffix,
                                 sample_subdirs)
    print(paste("Loading ORFs from:", tpms_file))
    orfs <- tpms_file %>% readr::read_tsv(comment = "#") %>% .$ORF
  } else {
    suppressMessages(library(Biostrings, quietly = T))
    print(paste("Loading ORFs from:", orf_fasta))
    orfs <- Biostrings::readDNAStringSet(orf_fasta) %>% names
  }
  orfs <- sort(orfs)
  tpm_list <- lapply(samples,
                     GetTpms,
                     samples_dir = samples_dir,
                     tpms_suffix = tpms_suffix,
                     sample_subdirs = sample_subdirs,
                     orfs = orfs)
  non_null_elts <- sapply(tpm_list, function(elt) !is.null(elt))
  names(tpm_list) <- samples
  collated_tpms <- dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])
  return(collated_tpms)
}

#' Collate TPMs from sample-specific files and saved collated TPMs.
#'
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
CollateTpms <- function(
  output_dir, tpms_file, sample_subdirs, orf_fasta, samples) {

  round1 <- function(x) round(x, digits = 1)
  tpms_file_path <- file.path(output_dir, tpms_file)
  write_provenance_header(get_Rscript_filename(), tpms_file_path)
  MakeTpmTable(output_dir, sample_subdirs, orf_fasta, samples) %>%
    dplyr::mutate_if(is.numeric, round1) %>%
    readr::write_tsv(tpms_file_path, col_names = TRUE, append = TRUE)
}

option_list <- list(
  make_option("--output-dir", type = "character", default = "./",
    help = "Output directory"),
  make_option("--tpms-file", type = "character",
    default = "TPMs_collated.tsv",
    help = "Name of collated TPMs file, relative to `output-dir`"),
  make_option("--sample-subdirs", type = "logical", default = FALSE,
    help = "Are samples in sample-specific subdirectories of `output-dir`?"),
  make_option("--orf-fasta", type = "character", default = NA,
    help = "ORF FASTA file that was aligned to and from which ORF names are to be retrieved"))

print_provenance(get_Rscript_filename())

opt <- parse_args(OptionParser(option_list = option_list),
                  positional_arguments = TRUE,
                  convert_hyphens_to_underscores = TRUE)
print("collate_tpms.R running with parameters:")
print(opt)

CollateTpms(opt$options$output_dir,
  opt$options$tpms_file,
  opt$options$sample_subdirs,
  opt$options$orf_fasta,
  opt$args)

print("collate_tpms.R done")
