library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(yaml)
library(optparse)
# library(Biostrings) # if we wanted to take ORF list from fasta

# Load yaml
option_list <- list( 
  make_option("--yaml", type="character", default=NULL,
              help="config file in yaml format")
)
opt <- parse_args(OptionParser(option_list=option_list))
yamlparams <- read_yaml(opt$yaml)

## for debugging in vignette
# setwd("~/Repos/RiboViz/")
# yamlparams <- read_yaml("vignette/vignette_config.yaml")

get_tpms <- function(ffile,ORFs) {
    # get tpm column from ffile
    # checking that gene names are as expected
    if(!file.exists(ffile)) {
        warning( paste(ffile, "does not exist, returning empty list") )
        return(NULL)
    }
    features_tab <- read_tsv(ffile)
    if(!all.equal(features_tab$ORF,ORFs)) {
        warning(paste("ORF names are not right in ", ffile))
    }
    return(features_tab$tpm)
}

get_tpms_bits <- function(fstem,ddir,fend,ORFs) {
    # get_tpms but putting the filename together
    get_tpms(paste0(ddir,"/",fstem,fend),ORFs)
}

# ffilename <- function(sample,dir_out) {}
# 
# make list of filenames

# collate
make_tpm_table <- function(yps,fend="_tpms.tsv") {
    # ORFs <- readDNAStringSet(yps$orf_fasta) %>% names # if ORF list from fasta
    ORFs <- paste0(yps$dir_out,"/",names(yps$fq_files)[1],fend) %>%
        read_tsv() %>%
        .$ORF
    tpm_list <- lapply(names(yps$fq_files), get_tpms_bits,
                         ddir=yps$dir_out,fend=fend,
                         ORFs=ORFs)
    non_null_elts <- sapply(tpm_list,function(elt) !is.null(elt))
    names(tpm_list) <-  names(yps$fq_files)
    bind_cols(ORF=ORFs,
              tpm_list[non_null_elts])
}

round1 <- function(x) round(x,digits=1)

make_tpm_table(yamlparams) %>%
    mutate_if(is.numeric,round1) %>%
    write_tsv(paste0(yamlparams$dir_out,"/","TPMs_collated.tsv"))