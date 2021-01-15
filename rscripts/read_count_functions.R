# read_count_functions.R
# R functions for dealing with data in h5 "ribogrid" format
# https://github.com/riboviz/riboviz/

# read in dependent packages
suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))
suppressMessages(library(rhdf5))
suppressMessages(library(parallel))
suppressMessages(library(optparse))
suppressMessages(library(RcppRoll))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(purrr))
suppressMessages(library(here))
#suppressMessages(library(devtools)) # TODO not sure if this needs to be enabled all the time, or just when developing?
#suppressMessages(library(roxygen2)) # TODO not sure if this needs to be enabled all the time, or just when developing?

#####

### Functions to read data from GFF ###


#' readGFFAsDf(): Read GFF file as a tibble (nicer dataframe)
#' 
#' Read in positions of all reatures in GFF format and convert to tibble data frame
#' 
#' @param orf_gff_file A filepath to a riboviz generated GFF2/GFF3 annotation file
#' 
#' @return Tidy data frame (tibble) of GFF data from GFF file 
#' 
#' @examples 
#' readGFFAsDf(orf_gff_file="vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#' 
#' @export
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)
#TEST: readGFFAsDf(x): creates tidy dataframe (tibble) = TRUE
#TEST: names(readGFFAs(x)) in: c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "phase", "Name"), but this may differ between GFFs?
#TEST: UTR names in 'type' match expected? (according to GFF definitions this 
 # should be something like threeprimeUTR or similar? But this isn't what we have in this file, which is UTR3)

#' GetCDS5start(): Get start coordinate for one named feature from GFF data frame.
#' 
#' If a named feature has multiple entries (e.g. multiplke exons), this returns the minimum.
#' 
#' @param name character; feature name.
#' @param gffdf data.frame or tibble; riboviz-format GFF in tidy data format, as created by readGFFAsDf()
#' @param ftype character; feature type to extract start location from GFF object for. Default: "CDS"
#' @param fstrand character; which strand of GFF file data to use. Default: "+"
#' 
#' @return Single 'start' value of transcript coordinate for named gene from CDS.
#' 
#' @examples 
#' gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#' GetCDS5start("YAL068C", gff_df)
#' 
#' @export
GetCDS5start <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(start) %>%  # pull() pulls out single variable
    min 
}
#TEST: GetCDS5start(): returns 1 integer value = TRUE

#' GetCDS3end(): Get end coordinate for one named feature from GFF data frame.
#' 
#' If a named feature has multiple entries (e.g. multiplke exons), this returns the maximum.
#' 
#' @param name character; feature name (usually, a gene)
#' @param gffdf data.frame or tibble; riboviz-format GFF in tidy data format, as created by readGFFAsDf()
#' @param ftype character; feature type to extract end location from GFF object for. Default: "CDS"
#' @param fstrand character; which strand of GFF file data to use. Default: "+"
#' 
#' @return single 'end' value of transcript coordinate for named gene from CDS
#' 
#' @examples 
#' gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#' GetCDS3end("YAL068C", gff_df)
#' 
#' @export
GetCDS3end <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(end) %>% 
    max 
}
#TEST: GetCDS3end(): returns 1 integer value = TRUE
#TEST: could compare 'width' is same as GetCDS3end - GetCDS5start ?

#####

### Functions to read data from h5 file ###


#' GetGeneDatamatrix(): Get matrix of read counts  by length for one gene and dataset from .h5 file
#' 
#' This accesses the attribute `reads/data` in .h5 file.
#' 
#' @param gene gene name to pull out read counts for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' 
#' @return numeric matrix of read count data for given gene in given dataset
#' 
#' @examples  
#' GetGeneDatamatrix(gene="YAL068C", dataset="vignette", hd_file="vignette/output/WTnone/WTnone.h5")
#' 
#' @export
GetGeneDatamatrix <- function(gene, dataset, hd_file){
  rhdf5::h5read(file = hd_file, name = paste0("/", gene, "/", dataset, "/reads/data")) %>%
    return()
}
#TEST: GetGeneDatamatrix(): returns matrix TRUE

#' GetGeneReadLength(): Get vector of read lengths for one gene and dataset from .h5 file
#' 
#' This accesses the attribute `reads/reads_by_len` in .h5 file.
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file  name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' 
#' @return vector of integer read lengths for specified gene in dataset from hd_file
#' 
#' @examples 
#' GetGeneReadLength(gene="YAL068C", dataset="vignette", hd_file="vignette/output/WTnone/WTnone.h5")
#' 
#' @export
GetGeneReadLength <- function(gene, dataset, hd_file){
  rhdf5::h5readAttributes(file = hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["reads_by_len"]]
}
#TEST: GetGeneReadLength(): returns numeric vector

#' GetGeneDatamatrix5start(): Get matrix of read counts by length between specific positions at 5' start for one gene and dataset from .h5 file and GFF data
#' 
#' Get matrix of read counts between specific positions 
#' (from n_buffer before start codon to nnt_gene after start codon)
#' for gene and dataset from hd5 file hd_file, using UTR5 annotations in gff
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file  name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param posn_5start numeric value, transcript-centric coordinate value for 
#' 5' location gene feature (e.g. CDS) starts from; ~equivalent to output of GetCDS5start()
#' @param n_buffer numeric value, number 'n' nucleotides of UTR buffer to include in metagene plots; riboviz default (set in generate_stats_figs.R): 25
#' @param nnt_gene numeric value, n nucleotides of gene to include in metagene plots; riboviz default (set in generate_stats_figs.R): 50
#' @param posn_3end numeric value, the 3'-end of the protein-coding sequence for rare situations where nnt_gene may lead to indexing error; default: -Inf
#' 
#' @return matrix of read counts for specific gene using .h5 and .gff information
#' 
#' @examples 
#' gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#' GetGeneDatamatrix5start(
#'   gene="YAL068C", 
#'   dataset="vignette",
#'   hd_file="vignette/output/WTnone/WTnone.h5",
#'   posn_5start = GetCDS5start("YAL068C", gff_df),
#'   n_buffer = 25,
#'   nnt_gene = 50)
#' 
#' @export
GetGeneDatamatrix5start <- function(gene, dataset, hd_file, 
                                    posn_5start, n_buffer, nnt_gene,posn_3end=Inf) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hd_file)
  
  # if n_buffer bigger than length n_utr5, pad with zeros:
  if (posn_5start > n_buffer) {
    # if posn_5start bigger than n_buffer
    n_left5 <- posn_5start - n_buffer # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  } else {
    # if length n_utr5 less than n_buffer
    n_left5 <- 1 # column to start from (5'end)
    zeropad5_mat <- matrix(
      0, 
      nrow = nrow(data_mat_all), 
      ncol = (n_buffer - posn_5start + 1 )
    )
  }
  n_right3 <- posn_5start + nnt_gene - 1 # column to end with (3'end)
  if (n_right3 > posn_3end) {
  
    zeropad3_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = n_right3 - posn_3end)
    n_right3 <- posn_3end
  } else{
    zeropad3_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  }

  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  x<-do.call("cbind",list(zeropad5_mat, data_mat_5start,zeropad3_mat))
  return(x)
}
#TEST: GetGeneDatamatrix5start(): returns a numeric matrix: TRUE
#TEST: GetGeneDatamatrix5start() : number of columns in matrix should be same as nnt_gene + n_buffer

#' GetGeneDatamatrix3end(): Get matrix of read counts by length between specific positions at 3' end for one gene and dataset from .h5 file and GFF data
#' 
#' get data matrix of read counts from nnt_gene before stop codon to n_buffer after
#' for gene and dataset from hd5 file hd_file, using UTR3 annotations in gff
#' if n_buffer bigger than length n_utr3, pad with zeros.
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file  name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param posn_3end numeric value, transcript-centric coordinate value for 
#' 3' location gene feature (e.g. CDS) ends at; ~equivalent to output of GetCDS3end()
#' @param n_buffer numeric value, number 'n' nucleotides of UTR buffer to include in metagene plots; riboviz default (set in generate_stats_figs.R): 25
#' @param nnt_gene numeric value, 'n' nucleotides of gene to include in metagene plots; riboviz default (set in generate_stats_figs.R): 50 
#' @param posn_5start numeric value, the 5'-end of the protein-coding sequence for rare situations where nnt_gene may lead to indexing error; default: -Inf
#' @return matrix of read counts for specific gene using .h5 and .gff information
#' 
#' @examples 
#' gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#' GetGeneDatamatrix5start(
#'   gene="YAL068C", 
#'   dataset="vignette",
#'   hd_file="vignette/output/WTnone/WTnone.h5",
#'   posn_3end = GetCDS3end("YAL068C", gff_df),
#'   n_buffer = 25,
#'   nnt_gene = 50
#'  )
#' 
#' @export
GetGeneDatamatrix3end <- function(gene, dataset, hd_file, 
                                  posn_3end,

                                  n_buffer, nnt_gene,posn_5start=-Inf) {

  # TODO: CHECK startpos/off-by-one
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hd_file)
  n_all <- ncol(data_mat_all)
  n_left5 <- posn_3end - nnt_gene + 1 # column to start from (5'end)
  if (n_left5 < posn_5start) {
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = posn_5start - n_left5)
    n_left5 <- posn_5start
  } else {
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  }
  n_utr3  <- n_all - posn_3end
  if (n_utr3 >= n_buffer) {
    # length n_utr3 bigger than n_buffer
    n_right3 <- posn_3end + n_buffer # column to end with (3'end)
    zeropad3_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  } else {
    # length n_utr3 less than n_buffer
    n_right3 <- n_all # column to end with (3'end)
    zeropad3_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = n_buffer - n_utr3)
  }
  data_mat_3end <- data_mat_all[, n_left5:n_right3]
  x<-do.call("cbind",list(zeropad5_mat, data_mat_3end,zeropad3_mat))
  return(x)
}
#TEST: GetGeneDatamatrix3end(): returns a numeric matrix: TRUE
#TEST: GetGeneDatamatrix3end() : number of columns in matrix should be same as nnt_gene + n_buffer


#' TidyDatamatrix(): Convert gene data matrix into tidy dataframe
#' 
#' Converts data matrix into readable tidy format with columns: `ReadLen, Pos, Counts`
#' to hold read lengths, position (in transcript-centric coordinates), and number of reads 
#' 
#' @param data_mat single data matrix of gene data in format created by GetGeneDatamatrix5start() and GetGeneDatamatrix3end().
#' @param startpos numeric value, start position along transcript-centric alignment for specific gene 
#' @param startlen numeric value, value from which to begin counting lengths from (ie equivalent to '10' in read_lengths= 10:50)
#' 
#' @return tidy format data frame (tibble), with columns: ReadLen, Pos and Counts
#' 
#' @examples 
#' gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#' datamatrix_YAL068C_5start <- GetGeneDatamatrix5start(
#'   gene="YAL068C", 
#'   dataset="vignette",
#'   hd_file="vignette/output/WTnone/WTnone.h5",
#'   posn_3end = GetCDS3end("YAL068C", gff_df),
#'   n_buffer = 25,
#'   nnt_gene = 50
#' )
#' TidyDatamatrix(data_mat=datamatrix_YAL068C_5start, startpos=1, startlen=1)
#'  
#' @export 
TidyDatamatrix <- function(data_mat, startpos = 1, startlen = 1) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(data_mat) - 1)
  readlengths <- startlen:(startlen + nrow(data_mat) - 1)
  data_mat %>%
    set_colnames(positions) %>%
    as_tibble() %>%
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts))
}
#TEST: TidyDatamatrix(): returns tidy format data frame (tibble)
#TEST: TidyDatamatrix() number of rows of output tibble = nrow(data_mat) * ncol(data_mat)
#TEST: TidyDatamatrix(): column names are %in% c("ReadLen", "Pos", "Counts")

#' AllGenes5StartPositionLengthCountsTibble(): Calculate sum of position- and read-length specific total counts over all genes at 5' start
#' 
#' TODO Variables needed for internal function GetGeneDatamatrix5start() in AllGenes5StartPositionLengthCountsTibble() aren't passed through using `...` function, but perhaps could be. 
#' 
#' @param gene_names vector of gene names to pull out data for (created early in generate_stats_figs.R by line: "gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name")
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data format, as created by readGFFAsDf()
#' 
#' @return tibble (tidy data frame) of metagene data (summed across all genes), 
#' with three columns: "ReadLen", "Pos", "Counts"
#' 
#' @examples
#' # Variables needed for internal function GetGeneDatamatrix5start() in AllGenes5StartPositionLengthCountsTibble() aren't passed through using `...` function, but perhaps could be. See issue riboviz/#248
#'  min_read_length <- 10
#'  max_read_length <- 50
#'  n_buffer <- 25
#'  nnt_buffer <- 25
#'  nnt_gene <- 50
#'  
#'  gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#'  gene_names <- rhdf5::h5ls(hd_file="vignette/output/WTnone/WTnone.h5", recursive = 1)$name
#'  AllGenes5StartPositionLengthCountsTibble(
#'    gene_names=gene_names, 
#'    dataset="vignette", 
#'    hd_file="vignette/output/WTnone/WTnone.h5", 
#'    gff_df
#'  )
#' 
#' @export
AllGenes5StartPositionLengthCountsTibble <- function(gene_names, dataset, hd_file, gff_df){
gene_poslen_counts_5start_df <-
  lapply(gene_names,
         function(gene) 
           GetGeneDatamatrix5start(gene,
                                   dataset,
                                   hd_file,
                                   posn_5start = GetCDS5start(gene, gff_df),
                                   posn_3end = GetCDS3end(gene,gff_df),
                                   n_buffer = nnt_buffer,
                                   nnt_gene = nnt_gene)
  ) %>%
  Reduce("+", .) %>% # sums the list of data matrices
  TidyDatamatrix(startpos = -nnt_buffer + 1, startlen = min_read_length) 

  return(gene_poslen_counts_5start_df)
} 
#TEST: AllGenes5StartPositionLengthCountsTibble(): returns tidy format data frame (tibble)
#TEST: AllGenes5StartPositionLengthCountsTibble(): tibble has 3 columns
#TEST: AllGenes5StartPositionLengthCountsTibble() number of rows of output tibble = nrow(data_mat) * ncol(data_mat)
#TEST: AllGenes5StartPositionLengthCountsTibble(): the column names are %in% c("ReadLen", "Pos", "Counts")
# gives: 
# > str(gene_poslen_counts_5start_df)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	3075 obs. of  3 variables:
#   $ ReadLen: int  10 11 12 13 14 15 16 17 18 19 ...
#   $ Pos    : int  -24 -24 -24 -24 -24 -24 -24 -24 -24 -24 ...
#   $ Counts : int  0 0 0 0 0 0 5 4 6 0 ...

#' AllGenes3EndPositionLengthCountsTibble(): Calculate sum of position- and read-length specific total counts over all genes at 3' end
#' 
#' @param gene_names vector of gene names to pull out data for (created early in generate_stats_figs.R by line: "gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name")
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data format, as created by readGFFAsDf()
#' 
#' @return tibble (tidy data frame) of metagene data (summed across all genes), 
#' with three columns: "ReadLen", "Pos", "Counts"
#' 
#' @examples
#'  min_read_length <- 10
#'  max_read_length <- 50
#'  n_buffer <- 25
#'  nnt_buffer <- 25
#'  nnt_gene <- 50
#'  
#'  gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#'  gene_names <- rhdf5::h5ls(hd_file="vignette/output/WTnone/WTnone.h5", recursive = 1)$name
#'  AllGenes3EndPositionLengthCountsTibble(
#'    gene_names=gene_names, 
#'    dataset="vignette", 
#'    hd_file="vignette/output/WTnone/WTnone.h5", 
#'    gff_df
#'  )
#' 
#' @export
AllGenes3EndPositionLengthCountsTibble <- function(gene_names, dataset, hd_file, gff_df){
  gene_poslen_counts_3end_df <-
    lapply(gene_names,
           function(gene)
             GetGeneDatamatrix3end(
               gene,
               dataset,
               hd_file,
               posn_3end = GetCDS3end(gene, gff_df),
               posn_5start = GetCDS5start(gene,gff_df),
               n_buffer = nnt_buffer,
               nnt_gene = nnt_gene
             )) %>%
    Reduce("+", .) %>% # sums the list of data matrices
    TidyDatamatrix(startpos = -nnt_gene + 1, startlen = min_read_length)
  
  return(gene_poslen_counts_3end_df)
} 
#TEST: AllGenes3EndPositionLengthCountsTibble(): returns tidy format data frame (tibble)
#TEST: AllGenes3EndPositionLengthCountsTibble(): tibble has 3 columns
#TEST: AllGenes3EndPositionLengthCountsTibble(): number of rows of output tibble = nrow(data_mat) * ncol(data_mat)
#TEST: AllGenes3EndPositionLengthCountsTibble(): the column names are %in% c("ReadLen", "Pos", "Counts")
# # gives:
# # > str(gene_poslen_counts_3end_df)
# # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	3075 obs. of  3 variables:
# #   $ ReadLen: int  10 11 12 13 14 15 16 17 18 19 ...
# #   $ Pos    : int  -49 -49 -49 -49 -49 -49 -49 -49 -49 -49 ...
# #   $ Counts : int  0 0 0 0 0 49 62 27 219 50 ...

#' plot_ribogrid(): Plot a ribogrid heatmap using position, read length and count data
#' 
#' @param tidymat tibble (tidy data frame) of metagene data (summed across all genes), 
#' with three columns: "ReadLen", "Pos", "Counts" ; as created by TidyDatamatrix().
#' Used for output of AllGenes5StartPositionLengthCountsTibble() and AllGenes3EndPositionLengthCountsTibble()
#' 
#' @return ggplot object; e.g start_codon_ribogrid_plot
#' 
#' @examples 
#'  min_read_length <- 10
#'  max_read_length <- 50
#'  n_buffer <- 25
#'  nnt_buffer <- 25
#'  nnt_gene <- 50
#'  
#'  gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#'  gene_names <- rhdf5::h5ls(hd_file="vignette/output/WTnone/WTnone.h5", recursive = 1)$name
#'  gene_poslen_counts_5start_df <- AllGenes5StartPositionLengthCountsTibble(
#'    gene_names=gene_names, 
#'    dataset="vignette", 
#'    hd_file="vignette/output/WTnone/WTnone.h5", 
#'    gff_df
#'  )
#'  
#' plot_ribogrid(gene_poslen_counts_5start_df)
#' 
#' @export
plot_ribogrid <- function(tidymat) {
  ggplot(data = tidymat, aes(x = Pos, y = ReadLen, fill = Counts)) +
    geom_tile() +
    scale_fill_gradient("count", low = "white", high = "darkblue") +
    theme(panel.grid = element_blank()) +
    labs(x = "position of read 5' end", y = "read length")
}
#TEST: plot_ribogrid(): does this return a ggplot object?
#TEST: plot_ribogrid(): ? are there any other good 'tests' for plot objects? TODO: Investigate this.

#' barplot_ribogrid(): plot ribogrid barplot using position, read length and count data
#' 
#' @param tidymat tibble (tidy data frame) of metagene data (summed across all genes), 
#' with three columns: "ReadLen", "Pos", "Counts" ; as created by TidyDatamatrix
#' Used for output of AllGenes5StartPositionLengthCountsTibble() or AllGenes3EndPositionLengthCountsTibble()
#' @param small_read_range numeric range of read lengths to plot barplots for (default: 26:32)
#' 
#' @return ggplot object; e.g start_codon_ribogrid_bar_plot
#' 
#' @examples 
#'  min_read_length <- 10
#'  max_read_length <- 50
#'  n_buffer <- 25
#'  nnt_buffer <- 25
#'  nnt_gene <- 50
#'  
#'  gff_df <- readGFFAsDf("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#'  gene_names <- rhdf5::h5ls(hd_file="vignette/output/WTnone/WTnone.h5", recursive = 1)$name
#'  gene_poslen_counts_5start_df <- AllGenes5StartPositionLengthCountsTibble(
#'    gene_names=gene_names, 
#'    dataset="vignette", 
#'    hd_file="vignette/output/WTnone/WTnone.h5", 
#'    gff_df
#'  )
#'  
#' barplot_ribogrid(gene_poslen_counts_5start_df)
#' 
#' @export
barplot_ribogrid <- function(tidymat, small_read_range = 26:32) {
  ggplot(
    data = filter(tidymat, ReadLen %in% small_read_range),
    aes(x = Pos, y = Counts)
  ) +
    geom_col() +
    facet_grid(ReadLen ~ ., scales = "free_y") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "position of read 5' end", y = "count")
}
#TEST: barplot_ribogrid(): does this return a ggplot object?
#TEST: barplot_ribogrid(): ? are there any other good 'tests' for plot objects? TODO: Investigate this.

#####
## functions for position specific distribution of reads

#' GetCodonPositionReads(): Get codon-position-specific reads
#' 
#' For Ribosome profiling datasets, assigns reads of lengths 28, 29, 30 to codons based on heuristic approach. NOTE: this function needs to be addressed/improved and replaced by using a-site calculations.
#' TODO: consider whether buffer should be an argument of this function, if left and right are ALWAYS based on "buffer +/- NN" formula as in CalculatePositionSpecificDistributionOfReads()
#' TODO: replace by better function
#'  
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param left integer, used in subsetting codon positions; e.g. "left = (buffer - 15)" as in CalculatePositionSpecificDistributionOfReads() usage
#' @param right integer, used in subsetting codon positions; e.g. "right = (buffer + 11)" as in CalculatePositionSpecificDistributionOfReads() usage
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' 
#' @return matrix of codon-specific reads
#' 
#' @examples 
#' buffer <- 250 # set in yaml arguments
#' min_read_length <- 10  # set in yaml arguments, default: 10
#'
#' GetCodonPositionReads(gene = "YAL068C", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", left = (buffer - 15), right = (buffer + 11), min_read_length = min_read_length)
#' # returns matrix: num [1:122] 0 0 0 0 0 0 0 0 0 0 ...
#'
#' @export
GetCodonPositionReads <- function(gene, dataset, hd_file = hd_file, left, right, min_read_length) {
  length_id <- 28 - min_read_length + 1
  reads_pos <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
  reads_pos_subset <- reads_pos[, left:(dim(reads_pos)[2] - right)] # Subset positions such that only CDS codon-mapped reads are considered
  end_reads_pos_subset <- ncol(reads_pos_subset) # Number of columns of the subset

  l28 <- RcppRoll::roll_suml(reads_pos_subset[length_id, 2:end_reads_pos_subset], n = 3, fill = NULL)[seq(1, length(reads_pos_subset[14, 2:end_reads_pos_subset]), 3)] # Map reads of length 28 to codons
  l29 <- RcppRoll::roll_suml(reads_pos_subset[(length_id + 1), 2:end_reads_pos_subset], n = 3, fill = NULL)[seq(1, length(reads_pos_subset[15, 2:end_reads_pos_subset]), 3)] # Map reads of length 29 to codons
  l30 <- RcppRoll::roll_suml(reads_pos_subset[(length_id + 2), 1:end_reads_pos_subset], n = 3, fill = NULL)[seq(1, length(reads_pos_subset[16, 1:end_reads_pos_subset]), 3)] # Map reads of length 30 to codons

  cod_sp_counts <- l28 + l29 + l30 # Sum of reads of lengths 28-30 at each codon
  cod_sp_counts <- cod_sp_counts[1:(length(cod_sp_counts) - 1)]
  return(cod_sp_counts)
}
#TEST: GetCodonPositionReads(): does it return a matrix? TRUE

#' GetMRNACoverage(): Calculate nucleotide-specific coverage 
#' 
#' Nt-specific coverage, suitable for mRNA datasets as does not include codon information.
#' TODO: consider whether buffer should be an argument of this function, if left and right are ALWAYS based on "buffer +/- NN" formula as in CalculateNucleotideBasedPositionSpecificReadsMRNA()
#' TODO: consider if we need this function.
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param left integer, used in subsetting codon positions; e.g. "left = (buffer - 49)" as in CalculateNucleotideBasedPositionSpecificReadsMRNA() usage
#' @param right integer, used in subsetting codon positions; e.g. "right = (buffer - 3)" as in CalculateNucleotideBasedPositionSpecificReadsMRNA() usage
#' @param read_range numeric range of read length values, set in generate_stats_figs.R from "min_read_length:max_read_length" values from yaml; Default: 10:50
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param buffer numeric, length of flanking region around the CDS in nucleotides, set in yaml config file. Default: 250
#' 
#' @return matrix of nucleotide-specific reads for mRNA datasets
#' 
#' @examples 
#' buffer <- 250 # set in yaml arguments
#' min_read_length <- 10  # set in yaml arguments, default: 10
#' read_range <- 10:50 # read_range set towards start of generate_stats_figs.R from yaml
#' 
#' GetMRNACoverage(gene = "YAL068C", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", left = (buffer - 49), right = (buffer - 3), read_range, min_read_length, buffer)
#' # returns matrix: num [1:368] 0 0 0 0 0 0 0 0 0 0 ...
#' 
#' @export 
GetMRNACoverage <- function(gene, dataset, hd_file, left, right, read_range, min_read_length, buffer) {
  # Get the matrix of read counts
  reads_pos <- GetGeneDatamatrix(gene, dataset, hd_file) 
  # Subset positions such that only CDS mapped reads are considered
  reads_pos_subset <- reads_pos[, left:(dim(reads_pos)[2] - right)] 

  nt_IR_list <- lapply(read_range, function(w) {
    IRanges::IRanges(start = rep(1:ncol(reads_pos_subset), reads_pos_subset[(w - min_read_length + 1), ]), width = w)
  }) # Create list of IRanges for position-specific reads of all length
  nt_IR <- unlist(as(nt_IR_list, "IRangesList")) # Combine IRanges from different read lengths
  nt_cov <- IRanges::coverage(nt_IR) # Estimate nt-specific coverage of mRNA reads

  # Subset coverage to only CDS
  nt_counts <- rep.int(S4Vectors::runValue(nt_cov), S4Vectors::runLength(nt_cov))
  if (length(nt_counts) >= (buffer - left)) {
    nt_counts <- nt_counts[(buffer - left):length(nt_counts)]
  } else {
    nt_counts <- 0
  }

  cds_length <- ncol(reads_pos_subset) - (buffer - left - 1) # Length of CDS
  nt_sp_counts <- rep(0, cds_length)

  if (length(nt_counts) < cds_length) {
    if (length(nt_counts) > 0) {
      nt_sp_counts[1:length(nt_counts)] <- nt_counts
    }
  } else {
    nt_sp_counts <- nt_counts[1:cds_length]
  }
  return(nt_sp_counts)
}
#TEST: GetMRNACoverage(): does it return a matrix? TRUE

#####

### functions to calculate read frame etc. ###

#' CalcAsiteFixedOneLength(): Calculate read A-site using a fixed displacement for a single read length
#' 
#' This is a helper function for `CalcAsiteFixed`.
#' 
#' @param reads_pos_length matrix of read lengths and positions (e.g. as given by GetGeneDatamatrix(gene, dataset, hd_file) )
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param read_range numeric range of read length values, set in generate_stats_figs.R from `min_read_length:max_read_length` values from yaml; Default: 10:50
#' @param asite_displacement numeric value giving read frame displacement from 5' end to A-site
#' 
#' @return matrix
#' 
#' @examples 
#' reads_pos_length <- GetGeneDatamatrix(gene = "YAL068C", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#'  # int [1:41, 1:863] 0 0 0 0 0 0 0 0 0 0 ...
#' 
#' CalcAsiteFixedOneLength(reads_pos_length, min_read_length = 10, read_length = 10:50, asite_displacement = 15)
#'  # num [1:41, 1:863] 0 0 0 0 0 0 0 0 0 0 ...
#' 
#' @export 
CalcAsiteFixedOneLength <- function(reads_pos_length, min_read_length,
                                    read_length, asite_displacement) {
  length_row_choose <- read_length - min_read_length + 1
  reads_pos_length[length_row_choose, ] %>%
    dplyr::lag(n = asite_displacement, default = 0)
}
#TEST: CalcAsiteFixedOneLength(): TODO

#' CalcAsiteFixed(): Calculate read A-site using a fixed displacement for fixed read lengths
#' 
#' The assignment rules are specified in a user-supplied data frame, `asite_displacement_length`.
#' 
#' @param reads_pos_length matrix of read lengths and positions (e.g. as given by GetGeneDatamatrix(gene, dataset, hd_file) )
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param asite_displacement_length data frame with columns `read_length` and `asite_displacement`
#'  default: read_length = c(28, 29, 30), and asite_displacement = c(15, 15, 15). 
#' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
#' 
#' @return numeric vector if colsum_out = TRUE; matrix with number of rows equivalent to number of rows in asite_displacement_length if colsum_out=FALSE
#' 
#' @examples 
#' reads_pos_length <- GetGeneDatamatrix(gene = "YAL068C", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#'  # int [1:41, 1:863] 0 0 0 0 0 0 0 0 0 0 ...
#'  
#' CalcAsiteFixed(reads_pos_length, min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE) 
#'  # num [1:863] 0 0 0 0 0 0 0 0 0 0 ...
#'  
#' @export 
CalcAsiteFixed <- function(reads_pos_length, min_read_length,
                           asite_displacement_length = data.frame(
                             read_length = c(28L, 29L, 30L),
                             asite_displacement = c(15L, 15L, 15L)
                           ),
                           colsum_out = TRUE) {
  npos <- ncol(reads_pos_length)
  Asite_counts_bylength <-
    purrr::map2(
      asite_displacement_length$read_length, asite_displacement_length$asite_displacement,
      function(read_length, asite_displacement) {
        CalcAsiteFixedOneLength(
          reads_pos_length,
          min_read_length,
          read_length,
          asite_displacement
        )
      }
    )
  if (colsum_out) {
    Asite_counts <- purrr::reduce(Asite_counts_bylength, `+`)
    return(Asite_counts)
  } else {
    # this has only as many columns as asite_displacement_length,
    # probably LESS than data_mat
    Asite_counts_bylengthmat <- unlist(Asite_counts_bylength) %>%
      matrix(ncol = npos, byrow = TRUE)
    return(Asite_counts_bylengthmat)
    # nrow(Asite_counts_bylengthmat) represents the rows in read_length column of asite_displacement_length. 
    # TODO: perhaps add row naming to make this ^ clear?
  }
}
#TEST: CalcAsiteFixed(): if col_sum=TRUE, return numeric vector? TRUE
#TEST: CalcAsiteFixed(): if col_sum=FALSE, return matrix? TRUE

#' SumByFrame(): Calculate sums of vector by 3nt frames 0, 1, 2
#' 
#' Note: this function will not operate correctly if used on outputs for CalcAsiteFixed() where col_sum=FALSE, as this would create a matrix, and input to x should be a vector.
#' 
#' @param x vector
#' @param left integer for starting position, frame 0; (e.g. 251 for gene YAL003W in S.cerevisiae in vignette dataset, as this is first NT of CDS)
#' @param right integer for ending position; (e.g. 871 for gene YAL003W in S.cerevisiae in vignette dataset, as this is last NT of CDS)
#'  
#' @return numeric vector of 3 values, sum of values for each of frame 0, 1, 2.
#'  
#' @examples 
#' SumByFrame(rep(1,9),1,9)
#' SumByFrame(rep(c(1,0,0),3),1,9)
#' SumByFrame(1:9,1,9)
#' 
#' reads_pos_length <- GetGeneDatamatrix(gene = "YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#'  # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#' reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE) 
#'  # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#' SumByFrame(reads_asitepos, left =251, right =871)
#'  # [1] 638  96 157
#'  
#' @export    
SumByFrame <- function(x, left, right) {
  positions_frame0 <- seq(left, right, 3) # positions used to pick out frame 0 reads
  sums_byframe <- c(
    x[ positions_frame0 ] %>% sum(),
    x[ positions_frame0 + 1 ] %>% sum(),
    x[ positions_frame0 + 2 ] %>% sum()
  )
  return(sums_byframe)
}
#TEST: SumByFrame(): returns numeric vector; TRUE
#TEST: SumByFrame(): length of numeric vector returned is 3 (ie frame 0, 1, 2); TRUE
#TODO: perhaps create a test on argument input for x to check it is vector, not matrix, (Note: this function will not operate correctly if used on outputs for CalcAsiteFixed() where col_sum=FALSE, as this would create a matrix, and input to x should be a vector.)

#' SnapToCodon(): Snap nucleotide-aligned reads to codon position
#' 
#' This takes a numeric vector (usually of read counts), and returns a numeric vector summed in groups of 3.
#' It is a thin wrapper of the rolling sum function, `RcppRoll::roll_suml`.
#' 
#' #TODO: clarify inclusion/exclusion of start/stop codon in left/right documentation, especially right. 
#' 
#' @param x numeric vector
#' @param left integer for starting position, frame 0; (e.g. 251 for gene YAL003W in S.cerevisiae in vignette dataset, as this is first NT of CDS)
#' @param right integer for ending position; (e.g. 871 for gene YAL003W in S.cerevisiae in vignette dataset, as this is last NT of CDS)
#' @param snapdisp integer any additional displacement in the snapping
#' 
#' @return numeric vector of summed read values per codon for specified region (left to right)
#' 
#' @examples 
#' SnapToCodon(rep(1,9),1,9)
#' SnapToCodon(rep(c(1,0,0),3),1,9)
#' SnapToCodon(1:9,1,9)
#'
#' reads_pos_length <- GetGeneDatamatrix(gene = "YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#'  # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#' reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE) 
#'  # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#' SnapToCodon(reads_asitepos, left=251, right=871, snapdisp=0L)
#'  # num [1:207] 0 9 0 0 0 3 0 0 0 0 ...
#' 
#' @export
SnapToCodon <- function(x, left, right, snapdisp=0L) {
  RcppRoll::roll_suml(x[(left:right) + snapdisp], n=3L, by=3L, fill = NULL)
}
#TEST: SnapToCodon(): returns numeric vector; TRUE
#TEST: SnapToCodon(): number of values (length of vector) returned is ~ (right - left)/3; TRUE
#TEST: SnapToCodon(rep(1,6),1,6) == c(3,3)

#' GetGeneCodonPosReads1dsnap(): Get codon positions and reads, snapped to codon, for one gene
#' 
#' This uses fixed A-site displacement from `CalcAsiteFixed`, then wraps `SnapToCodon`.
#' 
#' TODO: explain this function better
#' 
#' Note: GetGeneCodonPosReads1dsnap is not used in the riboviz pipeline (as of Jan 2021), 
#' but might be useful to test CalcAsiteFixed/SnapToCodon or for alternative visualizations.
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param left integer for starting position, frame 0; (e.g. 251 for gene YAL003W in S.cerevisiae in vignette dataset, as this is first NT of CDS)
#' @param right integer for ending position; (e.g. 871 for gene YAL003W in S.cerevisiae in vignette dataset, as this is last NT of CDS)
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param asite_displacement_length data frame of 2 columns: read_length default: c(28, 29, 30), and asite_displacement, default: c(15, 15, 15) 
#' @param snapdisp integer any additional displacement in the snapping
#' 
#' @return vector of summed read values per codon for a gene
#' 
#' @examples  
#' GetGeneCodonPosReads1dsnap(gene = "YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", left=251, right=871, min_read_length=10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), snapdisp=0L)
#'  # num [1:207] 0 9 0 0 0 3 0 0 0 0 ...
#' 
#' @export
GetGeneCodonPosReads1dsnap <- function(gene, dataset, hd_file, left, right, 
                         min_read_length, 
                         asite_displacement_length = data.frame(
                             read_length = c(28, 29, 30),
                             asite_displacement = c(15, 15, 15)
                           ), 
                         snapdisp=0L) {
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
  reads_asitepos <- CalcAsiteFixed(
    reads_pos_length, 
    min_read_length,
    asite_displacement_length
  )
  SnapToCodon(reads_asitepos, left, right, snapdisp)
}
#TEST: GetGeneCodonPosReads1dsnap(): returns numeric vector; TRUE
#TEST: GetGeneCodonPosReads1dsnap(): number of values (length of vector) returned is ~ (right - left)/3; TRUE

#' GatherByFrameCodon(): Gather a vector by 3nt frames 0, 1, 2, for each codon position
#' 
#' Note: there may be a better tidyverse-style implementation
#' 
#' @param x vector (e.g. reads_asitepos, as created by CalcAsiteFixed())
#' @param left integer for starting position, frame 0; (e.g. 251 for gene YAL003W in S.cerevisiae in vignette dataset, as this is first NT of CDS)
#' @param right integer for ending position; (e.g. 871 for gene YAL003W in S.cerevisiae in vignette dataset, as this is last NT of CDS)
#' 
#' @return tidy data frame (tibble) of codon position, counts in frame 0, 1, 2; columns: CodonPos, Ct_fr0, Ct_fr1, Ct_fr2
#' 
#' @examples 
#' GatherByFrameCodon(1:9,1,9)
#' 
#' reads_pos_length <- GetGeneDatamatrix(gene = "YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#'  # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#' reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE) 
#'  # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#' GatherByFrameCodon(reads_asitepos, left = 251, right = 871)
#'  # A tibble: 207 x 4
#' 
#' @export 
GatherByFrameCodon <- function(x, left, right) {
  positions_frame0 <- seq(left, right, 3) # positions to pick out frame 0 reads
  tibble(
    CodonPos = seq_len(length(positions_frame0)),
    Ct_fr0 = x[ positions_frame0 ],
    Ct_fr1 = x[ positions_frame0 + 1 ],
    Ct_fr2 = x[ positions_frame0 + 2 ]
  )
}
#TEST: GatherByFrameCodon(): returns tibble; TRUE
#TEST: GatherByFrameCodon(): columns of tibble are %in% c("CodonPos", "Ct_fr0", "Ct_fr1", "Ct_fr2"); TRUE
#TEST: GatherByFrameCodon(1:9,1,9): output is expected data frame.

#' CombinePValuesFisher(): Combine p-values using Fisher's 1-sided method 
#' 
#' This function is not used in riboviz codebase currently, we instead use the higher-powered `CombinePValuesStouffer`.
#'  
#' @param p numeric vector of input p-values.
#' 
#' @return combined p-value, numeric vector of length 1.
#'  
#' @examples
#' CombinePValuesFisher(0.1) 
#' CombinePValuesFisher(c(0.5,0.5)) 
#' CombinePValuesFisher(c(0.1,0.1)) 
#' CombinePValuesFisher(rep(0.1,10)) 
#' 
#' @export
#' @seealso CombinePValuesStouffer
CombinePValuesFisher <- function(p) {
  # Fisher's method (1-sided) to combine p-values
  pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = FALSE)
}
#TEST: CombinePValuesFisher(0.1) == 0.1
#TEST: CombinePValuesFisher(rep(0.5,10)) == 0.5

#' CombinePValuesStouffer(): Combine p-values using Stouffer's "inverse normal" 1-sided method.
#' 
#' @param p numeric vector of input p-values.
#' 
#' @return combined p-value, numeric vector of length 1.
#' 
#' @examples 
#' CombinePValuesStouffer(0.1) 
#' CombinePValuesStouffer(c(0.5,0.5)) 
#' CombinePValuesStouffer(c(0.1,0.1)) 
#' CombinePValuesStouffer(rep(0.1,10)) 
#' 
#' @export
#' @seealso CombinePValuesFisher
CombinePValuesStouffer <- function(p) {
  pnorm(sum(qnorm(p)) / sqrt(length(p)))
}
#TEST: CombinePValuesStouffer(0.1) == 0.1
#TEST: CombinePValuesStouffer(rep(0.5,10)) == 0.5

#' WilcoxTestFrame(): Calculate Wilcoxon rank-sum paired test of reads in frame 0 compared to 1 and 2.
#' 
#' This calculates p-values for a Wilcoxon rank-sum paired test of reads.
#' We compare reads in frame 0 with same-codon-paired frame 1, and respectively frame 2.
#' Then we combine p-values for frame 0 vs both frames 1 and 2, using Stoufer's method.
#' 
#' The test is based on that presented in the RiboCode paper:
#' 
#'   Xiao Z., Huang R., Xing X., Chen Y., Deng H., Yang X. De novo annotation and characterization of the translatome with ribosome profiling data. Nucleic Acids Res. 2018; 46:e61.
#' 
#' Unlike RiboCode, this function calculates the test values for fixed left/start and right/end positions.
#' RiboCode additionally searches for alternative start codons, which this function does not do.
#' 
#' @param x vector (e.g. reads_asitepos, as created by CalcAsiteFixed())
#' @param left integer for starting position, frame 0; (e.g. 251 for gene YAL003W in S.cerevisiae in vignette dataset, as this is first NT of CDS)
#' @param right integer for ending position; (e.g. 871 for gene YAL003W in S.cerevisiae in vignette dataset, as this is last NT of CDS)
#' 
#' @return named numeric vector with 3 values: "pval_fr0vs1" "pval_fr0vs2" "pval_fr0vsboth"
#' 
#' @examples 
#' reads_pos_length <- GetGeneDatamatrix(gene = "YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#'  # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#' reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE) 
#'  # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...#' 
#' WilcoxTestFrame(reads_asitepos, left=251, right =871)
#'  #  Named num [1:3] 3.65e-18 6.57e-13 6.03e-29
#'   # pval_fr0vs1    pval_fr0vs2 pval_fr0vsboth 
#'   # 3.649777e-18   6.574149e-13   6.025504e-29 
#' 
#' @export 
WilcoxTestFrame <- function(x, left, right) {
  gathered_by_frame <- GatherByFrameCodon(x, left, right)

  wtresults_fr0vs1 <-
    wilcox.test(
      x = gathered_by_frame$Ct_fr0,
      y = gathered_by_frame$Ct_fr1,
      alternative = "greater", paired = TRUE, exact = FALSE
    )
  wtresults_fr0vs2 <-
    wilcox.test(
      x = gathered_by_frame$Ct_fr0,
      y = gathered_by_frame$Ct_fr2,
      alternative = "greater", paired = TRUE, exact = FALSE
    )

  return(c(
    pval_fr0vs1 = wtresults_fr0vs1$p.value,
    pval_fr0vs2 = wtresults_fr0vs2$p.value,
    pval_fr0vsboth =
      CombinePValuesStouffer(c(
        wtresults_fr0vs1$p.value,
        wtresults_fr0vs2$p.value
      ))
  ))
}
#TEST: WilcoxTestFrame(): check for named numeric vector; TRUE
#TEST: WilcoxTestFrame(): length of returned vector is 3; TRUE

#' GetGeneReadFrame(): Calculate read frame statistics for one CDS or gene.
#' 
#' This:
#' * Assigns reads to A-sites by `CalcAsiteFixed`, using rules supplied in `asite_displacement_length`.
#' * Counts of reads in each frame, by calling `SumByFrame`.
#' * Calculates Wilcoxon rank-sum paired test p-values of reads in frame 0 vs 1, 2, or both, by calling `WilcoxTestFrame`
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param left integer for starting position, frame 0; (e.g. 251 for gene YAL003W in S.cerevisiae in vignette dataset, as this is first NT of CDS)
#' @param right integer for ending position; (e.g. 871 for gene YAL003W in S.cerevisiae in vignette dataset, as this is last NT of CDS)
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param asite_displacement_length data frame of 2 columns: read_length default: c(28, 29, 30), and asite_displacement, default: c(15, 15, 15) 
#' 
#' @return tibble of counts in each frame (0, 1, 2) per gene with combined (via Stouffer's inverse normal method) p-values from Wilcox test
#' 
#' @examples 
#' GetGeneReadFrame("YAL003W", dataset= "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", left=251, right=871, min_read_length=10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)))
#' 
#' @export 
#' @seealso [WilcoxTestFrame()], [CalcAsiteFixed()], [SumByFrame()]
GetGeneReadFrame <- function(gene, dataset, hd_file, left, right, min_read_length,
                             asite_displacement_length = data.frame(
                               read_length = c(28, 29, 30),
                               asite_displacement = c(15, 15, 15)
                             )) {
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  reads_asitepos <- CalcAsiteFixed(
    reads_pos_length, min_read_length,
    asite_displacement_length
  )
  sum_by_frame <- SumByFrame(reads_asitepos, left, right)
  wt_frame <- WilcoxTestFrame(reads_asitepos, left, right)
  tibble(
    gene = gene,
    Ct_fr0 = sum_by_frame[1],
    Ct_fr1 = sum_by_frame[2],
    Ct_fr2 = sum_by_frame[3],
    pval_fr0vs1 = wt_frame[1],
    pval_fr0vs2 = wt_frame[2],
    pval_fr0vsboth = wt_frame[3]
  )
}
#TEST: GetGeneReadFrame(): returns tibble; TRUE
#TEST: GetGeneReadFrame(): column names in c("gene", "Ct_fr0", "Ct_fr1", "Ct_fr2", "pval_fr0vs1", "pval_fr0vs2", "pval_fr0vsboth")

#' CalcReadFrameProportion(): Calculate read frame proportions from read frame counts
#' 
#' @param read_frame_df data frame of read frame counts per gene or feature, as output by GetGeneReadFrame().
#' This must have columns of the counts in each frame, `Ct_fr0, Ct_fr1, Ct_fr2.
#' It may additionally have any other columns, to desribe gene, CDS, or other features.
#' 
#' @return data frame, with additional columns for total counts and proportions.
#' * `Ct_all = Ct_fr0 + Ct_fr1 + Ct_fr2`
#' * `p_fr0 = Ct_fr0 / Ct_all`
#' * `p_fr1 = Ct_fr1 / Ct_all`
#' * `p_fr2 = Ct_fr2 / Ct_all`
#' 
#' @examples 
#' read_frame_df <- GetGeneReadFrame("YAL003W", dataset= "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", left=251, right=871, min_read_length=10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)))
#' CalcReadFrameProportion(read_frame_df)
#'  # A tibble: 1 x 11
#' 
#' @export
CalcReadFrameProportion <- function(read_frame_df) {
  stopifnot(all(c("Ct_fr0", "Ct_fr1", "Ct_fr2") %in% names(read_frame_df)))
  read_frame_df %>%
    mutate(
      Ct_all = Ct_fr0 + Ct_fr1 + Ct_fr2,
      p_fr0 = Ct_fr0 / Ct_all,
      p_fr1 = Ct_fr1 / Ct_all,
      p_fr2 = Ct_fr2 / Ct_all
    ) %>%
    return()
}
#TEST: CalcReadFrameProportion(): returns a data frame; TRUE
#TEST: CalcReadFrameProportion(): returns a tibble (TODO: test this, probably only returns tibble if read_frame_df is tibble?); TRUE
#TEST: CalcReadFrameProportion(): returned data frame contains these columns: c(Ct_all, p_fr0, p_fr1, p_fr2); TRUE

#' BoxplotReadFrameProportion(): Plot boxplot of proportion of reads in each frame, for each feature.
#' 
#' @param read_frame_df data frame of read frame proportions, with columns `p_fr0,p_fr1,p_fr2` 
#' and the value of `feat_names`.
#' For example, the output of [GetGeneReadFrame()].
#' @param feat_names character, name of feature being plotted, default: "gene"
#' 
#' @return ggplot object
#' 
#' @examples 
#' read_frame_df <- GetGeneReadFrame("YAL003W", dataset= "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", left=251, right=871, min_read_length=10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)))
#' BoxplotReadFrameProportion(read_frame_df, feat_names = "gene")
#' 
#' @export
BoxplotReadFrameProportion <- function(read_frame_df, feat_names = "gene") {
  rf_prop_long <- read_frame_df %>%
    CalcReadFrameProportion() %>%
    select(c(feat_names, "p_fr0", "p_fr1", "p_fr2")) %>%
    gather(-feat_names, key = "Frame", value = "Proportion") %>%
    mutate(Frame = factor(Frame,
      levels = c("p_fr0", "p_fr1", "p_fr2"),
      labels = 0:2
    ))
  ggplot(data = rf_prop_long, aes(x = Frame, colour = Frame, y = Proportion)) +
    geom_boxplot() +
    scale_y_continuous("Proportion, by feature", limits = c(0, 1), expand = c(0, 0)) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}
#TEST: BoxplotReadFrameProportion(): returns ggplot object: TRUE

#####
## biases in nucleotide composition along mapped read lengths

#' GetNTReadPosition(): Get a range of nucleotide read positions as IRanges object.
#' 
#' This function is only used as input to:
#' * `PositionSpecificConsensusMatrix`
#' * `CalculateBiasesInNucleotideComposition`
#' 
#' TODO: rewriting `CalculateBiasesInNucleotideComposition` in tidyverse will probably make this obsolete.
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param length_id integer vector of read length ids, allows subsetting for reads of particular length(s).
#' This is "length_id" not literal "length", because read lengths start at `min_read_length` not at 1. 
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' 
#' @return IRanges class object (from IRanges R package)
#' 
#' @examples 
#' GetNTReadPosition(gene ="YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", length_id = 30, min_read_length = 10)
#' 
#' @export
GetNTReadPosition <- function(gene, dataset, hd_file, length_id, min_read_length) {
  reads_pos_len <- GetGeneDatamatrix(gene, dataset, hd_file)[length_id, ] # Get reads of a particular length
  reads_pos_len <- reads_pos_len[1:(length(reads_pos_len) - (length_id + min_read_length - 1))] # Ignore reads whose 5' ends map close to the end of the 3' buffer
  pos <- rep(1:length(reads_pos_len), reads_pos_len) # nt positions weighted by number of reads mapping to it
  pos_IR <- IRanges::IRanges(start = pos, width = (length_id + min_read_length - 1)) # Create an IRanges object for position-specific reads of a particular length
  return(pos_IR)
}
#TEST: GetNTReadPosition(): returns IRanges object: TRUE (`class(IRanges::IRanges())`) 

#' PositionSpecificConsensusMatrix(): Create consensus matrix returning counts or frequencies based on position-specific reads of genes
#' 
#' This function is only used by `CalculateBiasesInNucleotideComposition`.
#' 
#' TODO: rewriting `CalculateBiasesInNucleotideComposition` in tidyverse will probably make this obsolete.
#' 
#' @param gene gene name 
#' @param pos_IR IRanges object (from IRanges R package) as output by GetNTReadPosition()
#' @param type character, type of operation to perform; "count" for counts (default), or "freq" for frequencies
#' @param cframe integer in c(0, 1, 2), coding frame 
#' @param length_id numeric value of length id to allow subsetting for reads of particular length
#' 
#' @return matrix of counts or frequencies by position for gene
#' 
#' @examples 
#' pos_IR <- GetNTReadPosition(gene ="YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", length_id = 30, min_read_length = 10)
#' PositionSpecificConsensusMatrix(gene ="YAL003W", pos_IR, type = "count", cframe = 0, length_id = 30)
#' 
#' @export
PositionSpecificConsensusMatrix <- function(gene, pos_IR, type = "count", cframe = 0L, length_id) {
  pos_IR_frame <- pos_IR[start(pos_IR) %% 3L == cframe] # Get position-specific reads of a particular length and ORF frame
  if (length(pos_IR_frame)) {
    # read in coding sequences
    coding_seqs <- readDNAStringSet(orf_fasta_file)
    pos_nt <- Biostrings::consensusMatrix(Biostrings::extractAt(coding_seqs[[gene]], pos_IR_frame))[1:4, ] # Get position-specific nucleotide counts
    if (type == "freq") {
      pos_nt <- pos_nt / colSums(pos_nt) # Select frequencies instead of counts
    }
  } else {
    pos_nt <- matrix(0, ncol = read_range[length_id], nrow = 4, dimnames = list(c("A", "C", "G", "T"), NULL))
  }
  return(pos_nt)
}
#TEST: PositionSpecificConsensusMatrix(): returns matrix; TRUE
#TEST: PositionSpecificConsensusMatrix(): number of rows = 4 (A, C, G, T); TRUE

#' CombineFrequencies(): Calculate normalized position-specific counts/freq of 
#' nucleotides A,C,G,T, across all genes, for reads of a particular length and frame.
#' 
#' TODO: Test and rewrite.
#' 
#' @param allfr matrix, multi-gene position-specific consensus matrix
#' 
#' @return TODO
#' 
#' @examples 
#' TODO: complete & test this example: may need to specify length_id?
#' gene_names <- rhdf5::h5ls(hd_file = "vignette/output/WTnone/WTnone.h5", recursive = 1)$name
#' fr0 <- mclapply(gene_names, function(gene) {PositionSpecificConsensusMatrix(gene = gene, pos_IR = out[[gene]], cframe = 0, length_id = length_id)}, mc.cores = num_processes)
#' allfr0 <- do.call(rbind, fr0)
#' CombineFrequencies(allfr0)
#' 
#' @export
CombineFrequencies <- function(allfr) {
  alph <- c("A", "C", "G", "T")
  nt_sp_freq <- c()
  for (i in alph) {
    nt_sp_freq <- rbind(nt_sp_freq, colSums(allfr[rownames(allfr) == i, ])) # Get position-specific counts/freq of a nt across all genes for reads of a particular length and frame
  }
  nt_sp_freq <- t(nt_sp_freq / colSums(nt_sp_freq)) # Convert total counts to freq OR freq to normalized frequencies
  colnames(nt_sp_freq) <- alph
  return(nt_sp_freq)
}
#TEST: CombineFrequencies(): TODO

#####
# read length and density functions

#' GetGeneLength(): Get read coding sequence (CDS) length for a gene, as stored in H5 file
#' 
#' TODO: this value is redundant with CDS length provided by the gff3 file.
#' 
#' @param gene gene name to get CDS length for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' 
#' @return integer, length of gene in nucleotides
#' 
#' @examples GetGeneLength(gene ="YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#'
#' @export
GetGeneLength <- function(gene, dataset, hd_file) {
  start_codon_pos <- rhdf5::h5readAttributes(file = hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["start_codon_pos"]][1]
  stop_codon_pos <- rhdf5::h5readAttributes(file = hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["stop_codon_pos"]][1]
  return(stop_codon_pos - start_codon_pos)
}
#TEST: GetGeneLength(): returns integer > 0; TRUE
#TEST: GetGeneLength(): could possibly check gff 'width' value for gene named against value returned by GetGeneLength()?

#' GetGeneReadsTotal(): Get total reads per gene from .h5 file
#'
#' Note, this just accesses the total value as stored in the .h5 file, it does not filter for a specific region.
#'
#' @param gene gene name to get reads counts for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' 
#' @return Total number of reads on `gene`, integer vector of length one.
#' 
#' @examples 
#' GetGeneReadsTotal(gene ="YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#' 
#' @export
GetGeneReadsTotal <- function(gene, dataset, hd_file) {
  rhdf5::h5readAttributes(file = hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["reads_total"]]
}
#TEST: GetGeneReadsTotal(): returns numeric; TRUE
#TEST: GetGeneReadsTotal(): returns integer >= 0; TRUE

#' GetGeneReadDensity(): Calculate density of reads per gene
#' 
#' TODO: It is ambiguous if this calculates reads per CDS, or reads per transcript. 
#' The argument `other_buffer` is ambiguous. This needs fixing.
#' 
#' @param gene gene name to get density of reads for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param other_buffer integer, default: 50. 
#' TODO: Fix this. Previously this was named 'buffer' but is not the same as yaml 'buffer' parameter with default of 250. Relevant to riboviz/#83 issue 
#' 
#' @return The density of reads on `gene`, total counts / length, a numeric vector of length one.
#' 
#' @examples 
#' GetGeneReadDensity(gene ="YAL003W", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5", other_buffer = 50)
#' 
#' @export 
GetGeneReadDensity <- function(gene, dataset, hd_file, other_buffer = 50L) {
  GetGeneReadsTotal(gene, dataset, hd_file) / (GetGeneLength(gene, dataset, hd_file) + other_buffer)
}
#TEST: GetGeneReadDensity(): returns numeric; TRUE