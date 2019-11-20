library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(optparse)

option_list <- list( 
  make_option("--dir_out", type="character", default="./",
              help="Output directory"),
  make_option("--file_out", type="character", default="TPMs_collated.tsv",
              help="Output file, relative to output directory"),
  make_option("--orf_fasta", type="character", default=NULL,
              help="ORF file that was aligned to")
)

parser <- OptionParser(option_list=option_list)
opts <- parse_args(parser, positional_arguments=TRUE)

dir_out <- opts$options$dir_out
file_out <- opts$options$file_out
samples <- opts$args
orf_fasta  <- opts$options$orf_fasta

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

get_tpms_bits <- function(fstem, ddir, fend, ORFs) {
    # get_tpms but putting the filename together
    get_tpms(paste0(ddir, "/", fstem, fend), ORFs)
}

# collate
make_tpm_table <- function(dir_out, samples, orf_fasta, fend="_tpms.tsv") {
    if(is.null(orf_fasta))
    {
         ORFs <- paste0(dir_out, "/", samples[1], fend) %>%
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
                       ddir=dir_out,
                       fend=fend,
                       ORFs=ORFs)
    non_null_elts <- sapply(tpm_list,function(elt) !is.null(elt))
    names(tpm_list) <- samples
    bind_cols(ORF=ORFs, tpm_list[non_null_elts])
}

round1 <- function(x) round(x,digits=1)

make_tpm_table(dir_out, samples, orf_fasta) %>%
    mutate_if(is.numeric, round1) %>%
    write_tsv(paste0(dir_out, "/", file_out))
