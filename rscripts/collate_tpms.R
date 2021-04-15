#' Collate TPMs

suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(getopt, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(optparse, quietly = T))
suppressMessages(library(purrr, quietly = T))
suppressMessages(library(readr, quietly = T))
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

LoadTpms <- function(ffile, orfs) {
  # Load data from ffile, check that gene names in 'ORF' column
  # equal orfs and return 'tpm' column.
  print(paste0("Loading TPMs from: ", ffile))
  if (!file.exists(ffile)) {
    warning(paste(ffile, "does not exist, returning empty list"))
    return(NULL)
  }
  features <- readr::read_tsv(ffile, comment = "#")
  if (!all.equal(features$ORF, orfs)) {
    warning(paste("ORF names are not right in ", ffile))
  }
  return(features$tpm)
}

GetTpmsFileName <- function(ddir, fstem, fend, sample_subdirs) {
  if (sample_subdirs) {
    file_name <- file.path(ddir, fstem, fend)
  } else {
    file_name <- file.path(ddir, paste0(fstem, "_", fend))
  }
  return(file_name)
}

GetTpms <- function(fstem, ddir, fend, sample_subdirs, orfs) {
  LoadTpms(GetTpmsFileName(ddir, fstem, fend, sample_subdirs), orfs)
}

MakeTpmTable <- function(output_dir, sample_subdirs, samples,
                         orf_fasta, fend="tpms.tsv") {
  # Collate TPMs into a table
  if (is.na(orf_fasta)) {
    orf_file <- GetTpmsFileName(output_dir,
                                samples[1],
                                fend,
                                sample_subdirs)
    print(paste("Loading ORFs from:", orf_file))
    orfs <- orf_file %>% readr::read_tsv(comment = "#") %>% .$ORF
  } else {
    print(paste("Loading ORFs from:", orf_fasta))
    # TODO untested
    library(Biostrings)
    orfs <- Biostrings::readDNAStringSet(orf_fasta) %>% names
  }
  tpm_list <- lapply(samples,
                     GetTpms,
                     ddir = output_dir,
                     fend = fend,
                     sample_subdirs = sample_subdirs,
                     orfs = orfs)
  non_null_elts <- sapply(tpm_list, function(elt) !is.null(elt))
  names(tpm_list) <- samples
  dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])
}

CollateTpms <- function(
  output_dir, tpms_file, sample_subdirs, orf_fasta, samples) {

  round1 <- function(x) round(x, digits = 1)
  tpms_file_path <- file.path(output_dir, tpms_file)
  write_provenance_header(get_Rscript_filename(), tpms_file_path)
  MakeTpmTable(output_dir, sample_subdirs, samples, orf_fasta) %>%
    mutate_if(is.numeric, round1) %>%
    write_tsv(tpms_file_path, col_names = TRUE, append = TRUE)
}

option_list <- list(
  make_option("--output-dir", type = "character", default = "./",
              help = "Output directory"),
  make_option("--tpms-file", type = "character",
              default = "TPMs_collated.tsv",
              help = "Output file, relative to output directory"),
  make_option("--sample-subdirs", type = "logical", default = FALSE,
              help = "Are samples in sample-specific subdirectories of output directory?"),
  make_option("--orf-fasta", type = "character", default = NA,
              help = "ORF file that was aligned to")
)

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
