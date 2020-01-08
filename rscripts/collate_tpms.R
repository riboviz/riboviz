library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(optparse)

option_list <- list( 
  make_option("--output-dir", type="character", default="./",
              help="Output directory"),
  make_option("--tpms-file", type="character", default="TPMs_collated.tsv",
              help="Output file, relative to output directory"),
  make_option("--sample-subdirs", type="logical", default=FALSE,
              help="Are samples in sample-specific subdirectories of output directory?"),
  make_option("--orf-fasta", type="character", default=NULL,
              help="ORF file that was aligned to")
)

parser <- OptionParser(option_list=option_list)

opts <- parse_args(parser, positional_arguments=TRUE, convert_hyphens_to_underscores=TRUE)

output_dir <- opts$options$output_dir
tpms_file <- opts$options$tpms_file
sample_subdirs <- opts$options$sample_subdirs
orf_fasta <- opts$options$orf_fasta
samples <- opts$args

get_tpms <- function(ffile, ORFs) {
    # get tpm column from ffile
    # checking that gene names are as expected
    if(!file.exists(ffile)) {
        warning(paste(ffile, "does not exist, returning empty list"))
        return(NULL)
    }
    features_tab <- read_tsv(ffile)
    if(!all.equal(features_tab$ORF,ORFs)) {
        warning(paste("ORF names are not right in ", ffile))
    }
    return(features_tab$tpm)
}

get_tpms_file_name <- function(ddir, fstem, fend, sample_subdirs) {
    print(paste0("sample_subdirs: ", sample_subdirs))
    if (sample_subdirs)
    {
        file_name = paste0(ddir, "/", fstem, "/", fstem, fend)
    }
    else
    {
        file_name = paste0(ddir, "/", fstem, fend)
    }
    print(paste0("file_name: ", file_name))
    return(file_name)
}

get_tpms_bits <- function(fstem, ddir, fend, sample_subdirs, ORFs) {
    # get_tpms but putting the filename together
    get_tpms(get_tpms_file_name(ddir, fstem, fend, sample_subdirs),
             ORFs)
}

# collate
make_tpm_table <- function(output_dir, sample_subdirs,
                           samples, orf_fasta,
                           fend="_tpms.tsv") {
    if(is.null(orf_fasta))
    {
         ORFs <- get_tpms_file_name(output_dir,
                                    samples[1],
                                    fend,
                                    sample_subdirs) %>%
            read_tsv() %>%
            .$ORF
    }
    else
    {
        # TODO untested
        print(paste("Taking ORFs from", orf_fasta))
        library(Biostrings)
        ORFs <- readDNAStringSet(orf_fasta) %>% names
    }
    tpm_list <- lapply(samples,
                       get_tpms_bits,
                       ddir=output_dir,
                       fend=fend,
                       sample_subdirs=sample_subdirs,
                       ORFs=ORFs)
    non_null_elts <- sapply(tpm_list,function(elt) !is.null(elt))
    names(tpm_list) <- samples
    bind_cols(ORF=ORFs, tpm_list[non_null_elts])
}

round1 <- function(x) round(x,digits=1)

make_tpm_table(output_dir, sample_subdirs, samples, orf_fasta) %>%
    mutate_if(is.numeric, round1) %>%
    write_tsv(paste0(output_dir, "/", tpms_file))
