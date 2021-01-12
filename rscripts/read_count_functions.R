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


#' readGFFAsDf(): Convert GFF to tidy dataframe
#' 
#' Read in positions of all exons/genes in GFF format and convert to tibble data frame
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

#' GetCDS5start(): Get start coordinate for gene in CDS from GFF
#' 
#' Extract start locations for each gene from GFF tidy dataframe for CDS only
#' 
#' @param name character; gene name
#' @param gffdf data.frame or tibble; riboviz-format GFF in tidy data format, as created by readGFFAsDf()
#' @param ftype character; gene feature to extract start location from GFF object for. Default: "CDS"
#' @param fstrand character; which strand of GFF file data to use. Default: "+"
#' 
#' @return Single 'start' value of transcriptome coordinate for named gene from CDS.
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

#' GetCDS3end(): Get end coordinate for each gene in CDS from GFF
#' 
#' Extract end locations for each gene from GFF tidy dataframe for CDS only
#' 
#' @param name character; gene name
#' @param gffdf data.frame or tibble; riboviz-format GFF in tidy data format, as created by readGFFAsDf()
#' @param ftype character; gene feature to extract end location from GFF object for. Default: "CDS"
#' @param fstrand character; which strand of GFF file data to use. Default: "+"
#' 
#' @return single 'end' value of transcriptome coordinate for named gene from CDS
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

# TODO: could GetCDS5start() and GetCDS3end() be rewritten/renamed/combined to allow 
# pulling out other features such as UTR5/UTR3 (though would need to ensure 
# feature naming in the GFF was always consistent!), or using 1 function to do 
# both ends depending on an extra argument?
#TEST: could compare 'width' is same as GetCDS3end - GetCDS5start ?

#####

### Functions to read data from h5 file ###


#' GetGeneDatamatrix(): Get matrix of read counts for gene from .h5 file
#' 
#' function to get data matrix of read counts for given gene and specified dataset from hd_file
#' 
#' @param gene gene name to pull out read counts for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' 
#' @return numeric matrix of read count data for given gene in given dataset (e.g. vignette)
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

#' GetGeneReadLength(): Get read length from .h5 file
#' 
#' function to get read length stored as attribute "reads_by_len" of 'reads' in H5 file
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

#' GetGeneDatamatrix5start(): Get matrix of read counts between specific positions at 5' start of gene using .h5 and .gff data
#' 
#' function to get matrix of read counts between specific positions 
#' (from n_buffer before start codon to nnt_gene after start codon)
#' for gene and dataset from hd5 file hd_file, using UTR5 annotations in gff
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file  name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param posn_5start numeric value, transcriptome-centric coordinate value for 
#' 5' location gene feature (e.g. CDS) starts from; ~equivalent to output of GetCDS5start()
#' @param n_buffer numeric value, number 'n' nucleotides of UTR buffer to include in metagene plots; riboviz default (set in generate_stats_figs.R): 25
#' @param nnt_gene numeric value, n nucleotides of gene to include in metagene plots; riboviz default (set in generate_stats_figs.R): 50 
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
                                    posn_5start, n_buffer, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hd_file)
  # @ewallace: replace this by gff_df?
  # n_utr5 <- BiocGenerics::width(gff[gff$type == "UTR5" & gff$Name == gene])
  
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
  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  return(cbind(zeropad5_mat, data_mat_5start))
}
#TEST: GetGeneDatamatrix5start(): returns a numeric matrix: TRUE
#TEST: GetGeneDatamatrix5start() : number of columns in matrix should be same as nnt_gene + n_buffer

#' GetGeneDatamatrix3end(): Get matrix of read counts between specific positions at 3' end of gene using .h5 and .gff data
#' 
#' get data matrix of read counts from nnt_gene before stop codon to n_buffer after
#' for gene and dataset from hd5 file hd_file, using UTR3 annotations in gff
#' if n_buffer bigger than length n_utr3, pad with zeros.
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file  name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param posn_3end numeric value, transcriptome-centric coordinate value for 
#' 3' location gene feature (e.g. CDS) ends at; ~equivalent to output of GetCDS3end()
#' @param n_buffer numeric value, number 'n' nucleotides of UTR buffer to include in metagene plots; riboviz default (set in generate_stats_figs.R): 25
#' @param nnt_gene numeric value, 'n' nucleotides of gene to include in metagene plots; riboviz default (set in generate_stats_figs.R): 50 
#' 
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
                                  n_buffer, nnt_gene) {
  # CHECK startpos/off-by-one
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hd_file)
  n_all <- ncol(data_mat_all)
  # @ewallace: replace this by gff_df?
  # n_utr3 <- BiocGenerics::width(gff[gff$type == "UTR3" & gff$Name == gene])
  n_left5 <- posn_3end - nnt_gene + 1 # column to start from (5'end)
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
  return(cbind(data_mat_3end, zeropad3_mat))
}
#TEST: GetGeneDatamatrix3end(): returns a numeric matrix: TRUE
#TEST: GetGeneDatamatrix3end() : number of columns in matrix should be same as nnt_gene + n_buffer


#' TidyDatamatrix(): Function to take gene data matrix and return tidy dataframe
#' 
#' converts data matrix into readable tidy format with columns: ReadLen, Pos and Counts
#' to hold read lengths, position (in transcriptome-centric coordinates), and number of reads 
#' 
#' @param data_mat single data matrix of gene data in format created by GetGeneDatamatrix5start() and GetGeneDatamatrix3end()
#' @param startpos numeric value, start position along transcriptome-centric alignment for specific gene 
#' @param startlen numeric value, value from which to begin counting lengths from (TODO: NOT CERTAIN OF THIS, CONFIRM!)
#' 
#' @return 
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

#' AllGenes5StartPositionLengthCountsTibble(): function to get gene/position-specific total counts for all read lengths for ALL genes for 5' end
#' 
#' get gene and position specific total counts for all read lengths, from 5' start
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

#' AllGenes3EndPositionLengthCountsTibble(): function to get gene/position-specific total counts for all read lengths for ALL genes for 3' end
#' 
#' get gene and position specific total counts for all read lengths, from 3' end
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

#' plot_ribogrid(): function to plot a ribogrid using position, read length and count data
#' 
#' take tidy format data across all genes and plot metagene ribogrid 
#' 
#' @param tidymat tibble (tidy data frame) of metagene data (summed across all genes), 
#' with three columns: "ReadLen", "Pos", "Counts" ; as created by 
#' AllGenes5StartPositionLengthCountsTibble() or AllGenes3EndPositionLengthCountsTibble()
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

#' barplot_ribogrid(): function to plot ribogrid barplot using position, read length and count data
#' 
#' take tidy format data across all genes and plot metagene ribogrid split by read length
#' 
#' @param tidymat tibble (tidy data frame) of metagene data (summed across all genes), 
#' with three columns: "ReadLen", "Pos", "Counts" ; as created by 
#' AllGenes5StartPositionLengthCountsTibble() or AllGenes3EndPositionLengthCountsTibble()
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

#' GetCodonPositionReads(): function to get codon-specific reads for RPF datasets
#' 
#' TODO: more info about what this function does
#' TODO: consider whether buffer should be an argument of this function, if left and right are ALWAYS based on "buffer +/- NN" formula as in CalculatePositionSpecificDistributionOfReads()
#' #' 
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

#' GetMRNACoverage(): function to return nucleotide-specific coverage for mRNA datasets
#' 
#' Nt-specific coverage for mRNA datasets
#' TODO: consider whether buffer should be an argument of this function, if left and right are ALWAYS based on "buffer +/- NN" formula as in CalculateNucleotideBasedPositionSpecificReadsMRNA()
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

#' CalcAsiteFixedOneLength(): function to find a-site for single length read using fixed displacement
#' 
#' Calculate read A-site using a fixed displacement for a single read length
#' 
#' @param reads_pos_length matrix of read lengths and positions
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param read_range numeric range of read length values, set in generate_stats_figs.R from "min_read_length:max_read_length" values from yaml; Default: 10:50
#' @param asite_displacement numeric value giving read frame displacement from 5' end to A-site
#' 
#' @return TODO
#' 
#' @examples TODO
#' 
#' @export 
CalcAsiteFixedOneLength <- function(reads_pos_length, min_read_length,
                                    read_length, asite_displacement) {
  length_row_choose <- read_length - min_read_length + 1
  reads_pos_length[length_row_choose, ] %>%
    dplyr::lag(n = asite_displacement, default = 0)
}
#TEST: CalcAsiteFixedOneLength(): TODO

#' CalcAsiteFixed(): function to find a-site for single length read using fixed displacement
#' 
#' Calculate read A-site using a fixed displacement for fixed read lengths
#' 
#' @param reads_pos_length matrix of read lengths and positions
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param asite_displacement_length data frame of 2 columns: read_length default: c(28, 29, 30), and asite_displacement, default: c(15, 15, 15) 
#' @param read_range numeric range of read length values, set in generate_stats_figs.R from "min_read_length:max_read_length" values from yaml; Default: 10:50
#' @param asite_displacement numeric value giving read frame displacement from 5' end to A-site
#' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
#' 
#' @return TODO
#' 
#' @examples TODO
#' 
#' @export 
CalcAsiteFixed <- function(reads_pos_length, min_read_length,
                           asite_displacement_length = data.frame(
                             read_length = c(28, 29, 30),
                             asite_displacement = c(15, 15, 15)
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
  }
}
#TEST: CalcAsiteFixed(): TODO

#' SumByFrame(): function to sum vector by 3nt frames 0, 1, 2
#' 
#'  @param x
#'  @param left
#'  @param right
#'  
#'  @return TODO
#'  
#'  @examples TODO
#'  
#'  @export    
SumByFrame <- function(x, left, right) {
  positions_frame0 <- seq(left, right, 3) # positions used to pick out frame 0 reads
  sums_byframe <- c(
    x[ positions_frame0 ] %>% sum(),
    x[ positions_frame0 + 1 ] %>% sum(),
    x[ positions_frame0 + 2 ] %>% sum()
  )
  return(sums_byframe)
}
#TEST: SumByFrame(): 

#' SnapToCodon(): function to snap nucleotide-aligned reads to codon position
#' 
#' @param x vector
#' @param left integer for starting position, frame 0
#' @param right integer for ending position
#' @param snapdisp integer any additional displacement in the snapping
#' 
#' @return TODO
#' 
#' @examples TODO
#' 
#' @export
SnapToCodon <- function(x, left, right, snapdisp=0L) {
  RcppRoll::roll_suml(x[(left:right) + snapdisp], n=3L, by=3L, fill = NULL)
}
#TEST: SnapToCodon(): TODO

#' NormByMean(): function to normalise data by dividing by mean
#' 
#' @param x vector
#' @param ... other arguments
#' 
#' @return TODO
#' 
#' @examples TODO
#'
#' @export
NormByMean <- function(x, ...) {
  x / mean(x, ...)
}
#TEST: NormByMean(): min value not below zero, TRUE
#TEST: NormByMean(): max value not above 1, TRUE

#' GetGeneCodonPosReads1dsnap(): function to get gene codon positions and reads, snapped to codon
#' 
#' TODO: explain this function better
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param left integer for starting position, frame 0
#' @param right integer for ending position
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' 
#' @return TODO
#' 
#' @examples TODO
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
#TEST: GetGeneCodonPosReads1dsnap(): TODO

#' GatherByFrameCodon(): function to gather a vector by 3nt frames 0,1,2 for each codon position
#' 
#' @param x vector
#' @param left integer for starting position, frame 0
#' @param right integer for ending position
#' 
#' @return tidy data frame (tibble) of codon position, counts in frame 0, 1, 2; columns: CodonPos, Ct_fr0, Ct_fr1, Ct_fr2
#' 
#' @examples TODO
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

#' CombinePValuesFisher(): function using Fisher's 1-sided method to combine p-values
#' 
#' TODO: this function IS NOT USED in riboviz codebase currently: Remove?
#' 
#' @param p numeric, p-value
#' 
#' @return TODO
#'  
#' @examples TODO
#' 
#' @export
CombinePValuesFisher <- function(p) {
  # Fisher's method (1-sided) to combine p-values
  pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = FALSE)
}
#TEST: CombinePValuesFisher(): TODO

#' CombinePValuesStouffer(): function using Stouffer's "inverse normal" 1-sided method to combine p-values
#' 
#' @param p numeric, p-value
#' 
#' @return TODO
#' 
#' @examples TODO
#' 
#' @export
CombinePValuesStouffer <- function(p) {
  pnorm(sum(qnorm(p)) / sqrt(length(p)))
}
#TEST: CombinePValuesStouffer(): TODO

#' WilcoxTestFrame(): function to run Wilcox rank-sum paired test to check frame 0 has more reads
#' 
#' @param x vector
#' @param left integer for starting position, frame 0
#' @param right integer for ending position
#' 
#' @return TODO
#' 
#' @examples TODO
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
#TEST: WilcoxTestFrame(): TODO

#' GetGeneReadFrame(): function to find frame for gene
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' 
#' @return tibble of counts in each frame (0, 1, 2) per gene with combined (via Stouffer's inverse normal method) p-values from Wilcox test
#' 
#' @examples TODO
#' 
#' @export 
GetGeneReadFrame <- function(gene, dataset, hd_file, left, right, min_read_length,
                             asite_displacement_length = data.frame(
                               read_length = c(28, 29, 30),
                               asite_displacement = c(15, 15, 15)
                             )) {
  # example from vignette:
  #   GetGeneReadFrame(hd_file, "YAL003W", dataset, 251, 871, min_read_length)
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

#' CalcReadFrameProportion(): function to calculate read frame proportions from read frame counts
#' 
#' @param read_frame_df data frame of read frame counts
#' 
#' @return data frame
#' 
#' @examples TODO
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

#' BoxplotReadFrameProportion(): function to plot proportion of read frames as boxplot
#' 
#' @param read_frame_df data frame of read frame counts
#' @param feat_names character, name of feature being plotted
#' 
#' @return TODO
#' 
#' @examples TODO
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

#' GetNTReadPosition(): function to get nucleotide read positions
#' 
#' @param gene gene name to get read lengths for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param length_id numeric value of length id to allow subsetting for reads of particular length
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' 
#' @return IRanges object (from IRanges R package)
#' 
#' @examples TODO
#' 
#' @export
GetNTReadPosition <- function(gene, dataset, hd_file, length_id, min_read_length) {
  reads_pos_len <- GetGeneDatamatrix(gene, dataset, hd_file)[length_id, ] # Get reads of a particular length
  reads_pos_len <- reads_pos_len[1:(length(reads_pos_len) - (length_id + min_read_length - 1))] # Ignore reads whose 5' ends map close to the end of the 3' buffer
  pos <- rep(1:length(reads_pos_len), reads_pos_len) # nt positions weighted by number of reads mapping to it
  pos_IR <- IRanges::IRanges(start = pos, width = (length_id + min_read_length - 1)) # Create an IRanges object for position-specific reads of a particular length
  return(pos_IR)
}
#TEST: GetNTReadPosition(): returns IRanges object: TRUE (not sure if this is easily testable?)

#' PositionSpecificConsensusMatrix(): function to create consensus matrix returning counts or frequencies based on position-specific reads of genes
#' 
#' @param gene gene name 
#' @param pos_IR IRanges object (from IRanges R package) as output by GetNTReadPosition()
#' @param type character, type of operation to perform; "count" for counts (default), or "freq" for frequencies
#' @param cframe numeric in c(0, 1, 2), coding frame 
#' @param length_id numeric value of length id to allow subsetting for reads of particular length
#' 
#' @return matrix of counts or frequencies by position for gene
#' 
#' @examples TODO
#' 
#' @export
PositionSpecificConsensusMatrix <- function(gene, pos_IR, type = "count", cframe = 0, length_id) {
  pos_IR_frame <- pos_IR[start(pos_IR) %% 3 == cframe] # Get position-specific reads of a particular length and ORF frame
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

#' CombineFrequencies(): function to get position-specific frequencies across all genes
#' 
#' @param allfr matrix, multi-gene position-specific consensus matrix
#' 
#' @return TODO
#' 
#' @examples TODO
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

#' GetGeneLength(): function to read length-specific total counts stored as attributes of 'reads_total' in H5 file
#' 
#' @param gene gene name to get reads counts for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' 
#' @return numeric, length of gene in nucleotides
#' 
#' @examples TODO
#'
#' @export
GetGeneLength <- function(gene, dataset, hd_file) {
  start_codon_pos <- rhdf5::h5readAttributes(file = hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["start_codon_pos"]][1]
  stop_codon_pos <- rhdf5::h5readAttributes(file = hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["stop_codon_pos"]][1]
  return(stop_codon_pos - start_codon_pos)
}
#TEST: GetGeneLength(): returns integer numeric; TRUE
#TEST: GetGeneLength(): could possibly check gff 'width' value for gene named against value returned by GetGeneLength()?
#TEST: GetGeneLength(): returns numeric greater than 0

#' GetGeneReadsTotal(): function to get total reads per gene from .h5 file
#' 
#' @param gene gene name to get reads counts for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' 
#' @return numeric value, total number of reads per named gene
#' 
#' @examples TODO
#' 
#' @export
GetGeneReadsTotal <- function(gene, dataset, hd_file) {
  rhdf5::h5readAttributes(file = hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["reads_total"]]
}
#TEST: GetGeneReadsTotal(): returns integer numeric; TRUE
#TEST: GetGeneReadsTotal(): returns integer >= 0; TRUE

#' GetGeneReadDensity(): function to get density of reads per gene
#' 
#' @param gene gene name to get density of reads for
#' @param dataset name of dataset stored in .h5 file
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples
#' @param other_buffer numeric, default: 50. TODO: Rename sensibly. Previously this was named 'buffer' but does not relate to yaml 'buffer' parameter with default of 250. Relevant to riboviz/#83 issue 
#' 
#' @return TODO
#' 
#' @examples TODO
#' 
#' @export 
GetGeneReadDensity <- function(gene, dataset, hd_file, other_buffer = 50) {
  GetGeneReadsTotal(gene, dataset, hd_file) / (GetGeneLength(gene, dataset, hd_file) + other_buffer)
}
#TEST: GetGeneReadDensity(): TODO