suppressMessages(library(getopt, quietly=T))
# Determine location of	provenance.R relative to current file
source(file.path(dirname(getopt::get_Rscript_filename()), "provenance.R"))
suppressMessages(library(tidyr, quietly=T))
suppressMessages(library(dplyr, quietly=T))
suppressMessages(library(readr, quietly=T))
suppressMessages(library(purrr, quietly=T))
suppressMessages(library(optparse, quietly=T))

option_list <- list(
  make_option("--output-dir",
              type = "character",
              default = "./",
              help = "Output directory"),
  make_option("--tpms-file",
              type = "character",
              default = "TPMs_collated.tsv",
              help = "Output file, relative to output directory"),
  make_option("--sample-subdirs",
              type = "logical",
              default = FALSE,
              help = "Are samples in sample-specific subdirectories of output directory?"),
  make_option("--orf-fasta",
              type = "character",
              default = NULL,
              help = "ORF file that was aligned to")
)

print(get_version(get_Rscript_filename()))
parser <- OptionParser(option_list = option_list)

opts <- parse_args(parser,
                   positional_arguments = TRUE,
                   convert_hyphens_to_underscores = TRUE)

print("collate_tpms.Rrunning with parameters:")
opts

output_dir <- opts$options$output_dir
tpms_file <- opts$options$tpms_file
sample_subdirs <- opts$options$sample_subdirs
orf_fasta <- opts$options$orf_fasta
samples <- opts$args

load_tpms <- function(ffile, orfs) {
    # Load data from ffile, check that gene names in 'ORF' column
    # equal orfs and return 'tpm' column.
    print(paste0("Loading TPMs from: ", ffile))
    if (!file.exists(ffile)) {
        warning(paste(ffile, "does not exist, returning empty list"))
        return(NULL)
    }
    features <- read_tsv(ffile)
    if (!all.equal(features$ORF, orfs)) {
        warning(paste("ORF names are not right in ", ffile))
    }
    return(features$tpm)
}

get_tpms_file_name <- function(ddir, fstem, fend, sample_subdirs) {
    if (sample_subdirs) {
        file_name <- file.path(ddir, fstem, fend)
    } else {
        file_name <- file.path(ddir, paste0(fstem, "_", fend))
    }
    return(file_name)
}

get_tpms <- function(fstem, ddir, fend, sample_subdirs, orfs) {
    load_tpms(get_tpms_file_name(ddir, fstem, fend, sample_subdirs),
              orfs)
}

make_tpm_table <- function(output_dir,
                           sample_subdirs,
                           samples,
                           orf_fasta,
                           fend="tpms.tsv") {
    # Collate TPMs into a table
    if (is.null(orf_fasta)) {
         orf_file <- get_tpms_file_name(output_dir,
                                        samples[1],
                                        fend,
                                        sample_subdirs)
        print(paste("Loading ORFs from:", orf_file))
        orfs <- orf_file %>% read_tsv() %>% .$ORF
    } else {
        print(paste("Loading ORFs from:", orf_fasta))
        # TODO untested
        library(Biostrings)
        orfs <- readDNAStringSet(orf_fasta) %>% names
    }
    tpm_list <- lapply(samples,
                       get_tpms,
                       ddir = output_dir,
                       fend = fend,
                       sample_subdirs = sample_subdirs,
                       orfs = orfs)
    non_null_elts <- sapply(tpm_list, function(elt) !is.null(elt))
    names(tpm_list) <- samples
    bind_cols(ORF = orfs, tpm_list[non_null_elts])
}

round1 <- function(x) round(x, digits = 1)

make_tpm_table(output_dir, sample_subdirs, samples, orf_fasta) %>%
    mutate_if(is.numeric, round1) %>%
    write_tsv(file.path(output_dir, tpms_file))
