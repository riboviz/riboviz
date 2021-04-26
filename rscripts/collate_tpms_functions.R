#' Collate TPMs functions.
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
#' See `CollateTpms` below.
#'
#' @export

suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(here, quietly = T))
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
  self <- "collate_tpms_functions.R"
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
    warning(paste("Inconsistent ORF names in", tpms_file))
  }
  return(features$tpm)
}

#' Collate TPMs from sample-specific files and return collated TPMs.
#'
#' If `orf_fasta` is provided then gene names are loaded from this
#' file, else they are loaded from the first sample-specific TPMs
#' file's `ORF` column. ORFs are sorted after loading.
#'
#' @param orf_fasta ORF FASTA file that was aligned to and from which
#' ORF names are to be retrieved (character).
#' @param samples List of sample files (where `names` attribute of
#' `samples` are the sample names) (list).
#' @return Collated TPMs with an `ORF` column and sample-specific
#' columns, named by sample name, with the TPMs values for each sample
#' (data.frame).
#'
#' @export
MakeTpmTable <- function(orf_fasta, samples) {
  if (is.na(orf_fasta)) {
    tpms_file <- samples[[1]]
    print(paste("Loading ORFs from:", tpms_file))
    orfs <- tpms_file %>% readr::read_tsv(comment = "#") %>% .$ORF
  } else {
    suppressMessages(library(Biostrings, quietly = T))
    print(paste("Loading ORFs from:", orf_fasta))
    orfs <- Biostrings::readDNAStringSet(orf_fasta) %>% names
  }
  orfs <- sort(orfs)
  tpm_list <- lapply(samples,
                     LoadTpms,
                     orfs = orfs)
  non_null_elts <- sapply(tpm_list, function(elt) !is.null(elt))
  names(tpm_list) <- names(samples)
  collated_tpms <- dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])
  return(collated_tpms)
}

#' Collate TPMs from sample-specific files and saved collated TPMs.
#'
#' @param tpms_file Name of collated TPMs file (character).
#' @param orf_fasta ORF FASTA file that was aligned to and from which
#' ORF names are to be retrieved (character).
#' @param samples List of sample files (where `names` attribute of
#' `samples` are the sample names) (list).
#'
#' @export
CollateTpms <- function(tpms_file, orf_fasta, samples) {
  round1 <- function(x) round(x, digits = 1)
  write_provenance_header(get_Rscript_filename(), tpms_file)
  MakeTpmTable(orf_fasta, samples) %>%
    dplyr::mutate_if(is.numeric, round1) %>%
    readr::write_tsv(tpms_file, col_names = TRUE, append = TRUE)
}
