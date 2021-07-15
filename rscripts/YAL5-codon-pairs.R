
# YAL5_h5 is at location $HOME/riboviz/riboviz/Mok-simYAL5/output/A/A.h5

# Given an h5 file, GFF file and .tsv file, this script creates a metafeature plot 
# for the codon/codon pair of interest. 

print('Starting process')

# source packages and functions from rscripts 
suppressMessages(source(here::here("rscripts", "read_count_functions.R")))
suppressMessages(source(here::here("rscripts", "stats_figs_block_functions.R")))


suppressMessages(library(tidyverse))
suppressMessages(library(zoo))
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))


parser <- ArgumentParser()

parser$add_argument('-i', '--input', help='Path input to h5 file')
parser$add_argument('-d', '--dataset', help='Name of the dataset being studied')
parser$add_argument('-g', '--gff', help='Path to the GFF3 file of the organism being studied')
parser$add_argument('-a', '--annotation', help='Path to codon table for organism')
parser$add_argument('--codonpair', help='Codon pair of interest')
parser$add_argument('-o', '--output', help='Path to output directory')
parser$add_argument('--expand_width', help='the desired range either side of the codon of interest', default = 5)
parser$add_argument('--startpos', help='position of the start codon', default = 1)
parser$add_argument('--startlen', help='smallest length of reads', default = 10) # please correct if wrong
parser$add_argument('--frame', help='frame to be studied', default = 0)
parser$add_argument('--minreadlen', help='minimum read length', default = 10)

args <- parser$parse_args()


hd_file <- args$input
dataset <- args$dataset
gff <- args$gff
yeast_codon_table <- args$annotation
codon_pair_of_interest <- args$codonpair
output_dir <- args$output
expand_width <- args$expand_width
startpos <- args$startpos
startlen <- args$startlen
filtering_framw <- args$frame
min_read_length <- args$minreadlen


# hd_file <- here::here("Mok-simYAL5", "output", "A", "A.h5")
# dataset <- "Mok-simYAL5"
# codon_of_interest <- 'TCC AAG'
# expand_width = 5L
# startpos <-1 
# startlen <- 10 
# filtering_frame <- 0
# min_read_length <- 10
# # gene_poscodon_codon_i200e::here("data", "yeast_codon_table.tsv")
# gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3")


### Load the necessary input files ###

# YAL5_h5 <- here::here("Mok-simYAL5", "output", "A", "A.h5")
# The h5 file from the dataset of interest 

# hd_file <- YAL5_h5
# assign YAL5_h5 file to a general name applied throughput the script 
# the script does not have to change based on different hd_files being used

gff_df <- readGFFAsDf(gff)
# The GFF file for the simulated dataset, given the general name gff_df so that
# the script does not have to change based on different gff_df files being used
gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name


# Import the .tsv file: 
yeast_codon_pos_i200 <- suppressMessages(readr::read_tsv(file = yeast_codon_table))

  # head(yeast_codon_pos_i200) 
  # 
  # # A tibble: 6 x 3
  # Gene    PosCodon Codon
  # <chr>      <dbl> <chr>
  #   1 YAL068C        1 ATG  
  # 2 YAL068C        2 GTC  
  # 3 YAL068C        3 AAA  
  # 4 YAL068C        4 TTA  
  # 5 YAL068C        5 ACT  
  # 6 YAL068C        6 TCA


# The yeast_codon_pos_i200 file is configured to show the lead and lag codon pair positions:

gene_poscodon_codon_i200 <- tibble::tibble(
  Gene = yeast_codon_pos_i200$Gene,
  CodonPos_1 = yeast_codon_pos_i200$PosCodon, 
  CodonPos_2 = dplyr::lead(yeast_codon_pos_i200$PosCodon),
  CodonPair = paste(yeast_codon_pos_i200$Codon, dplyr::lead(yeast_codon_pos_i200$Codon))
)
  
  # > gene_poscodon_codon_i200
  # # A tibble: 2,826,757 x 4
  # Gene    CodonPos_1 CodonPos_2 CodonPair
  # <chr>        <dbl>      <dbl> <chr>    
  #   1 YAL068C          1          2 ATG GTC  
  # 2 YAL068C          2          3 GTC AAA  
  # 3 YAL068C          3          4 AAA TTA  
  # 4 YAL068C          4          5 TTA ACT  
  # 5 YAL068C          5          6 ACT TCA  
  # 6 YAL068C          6          7 TCA ATC  
  # 7 YAL068C          7          8 ATC GCC  
  # 8 YAL068C          8          9 GCC GCT  
  # 9 YAL068C          9         10 GCT GGT  
  # 10 YAL068C         10         11 GGT GTC  
  # # ... with 2,826,747 more rows
  
  # Replaced yeast_codon_pos_i200 with gene_poscodon_codon_i200 throughout 
  # the script when applied to codon pairs 


# # Filter down the gene_poscodon_codon_i200 file to gene of interest
# YAL003W_pos <- dplyr::filter(gene_poscodon_codon_i200, Gene=="YAL003W")


#####


### Functions to read data from h5 file for all genes ###

# #' GetAllGeneDatamatrix(): Get matrix of read counts by length for all genes contained within the gene_names file from .h5 file.
# #'
# #' This accesses the attribute `reads/data` in .h5 file for each gene listed in gene_names. 
# #' 
# #' @param dataset name of dataset stored in .h5 file.
# #' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
# #' 
# #' @ return a list of numeric matrices of read count data for each given gene in given dataset.
# #'
# #' @examples
# #' GetAllGeneDatamatrix(dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"))
# #'
# #' @export
# # GetAllGeneDatamatrix <- function(dataset, hd_file){
# # 
# #   gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
# # 
# #   get_all_gene_datamatrix <- purrr::map(
# #       .x = gene_names,
# #       .f = GetGeneDatamatrix,
# #       dataset,
# #       hd_file
# #       )
# 
# #  #TEST: GetAllGeneDatamatrix(): returns list of matrices (where number of genes in gene_names = number of items inside list) TRUE
# gives:
# > str(get_all_gene_datamatrix)
# List of 5
#   $ : int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#   $ : int [1:41, 1:2429] 0 0 0 0 0 0 0 0 0 0 ...
#   $ : int [1:41, 1:1685] 0 0 0 0 0 0 0 0 0 0 ...
#   $ : int [1:41, 1:3509] 0 0 0 0 0 0 0 0 0 0 ...
#   $ : int [1:41, 1:2003] 0 0 0 0 0 0 0 0 0 0 ...
         
        

  
# TidyDatamatrixForAllGenes(): Convert gene data matrices for all genes into tidy dataframes
# 
# Converts each data matrix into readable tidy format with columns: `ReadLen`, `Pos`, `Counts`
# to hold read lengths, position (in transcript-centric coordinates), and number of reads for all genes
# 
# param dataset name of dataset stored in .h5 file.
# param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
# param startpos numeric value, start position along transcript-centric alignment for specific gene, default = 1.
# param startlen numeric value, value from which to begin counting lengths from (ie equivalent to '10' in read_lengths= 10:50), default = 10.
# 
# return a list of tidy format data frames (tibbles), with columns: `ReadLen`, `Pos` and `Counts`
# 
# examples
# TidyDatamatrixForAllGenes(dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), startpos = 1, startlen = 10)
# 
# export
  # # TidyDatamatrixForAllGenes <- function(dataset, hd_file, startpos = 1, startlen = 10){
  # #
  # #    gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
  # #
  # #    .x = get_all_gene_datamatrix <- GetAllGeneDatamatrix(dataset = dataset,
  # #                                                         hd_file = hd_file)
  # #
  # #    tidy_all_genes_datamatrix <- purrr::map(
  # #        .x = get_all_gene_datamatrix,
  # #        .f = TidyDatamatrix,
  # #        startpos = 1,
  # #        startlen = 10
  # #        )
  # #
  # #    return(tidy_all_genes_datamatrix)
  # #}
#TEST: TidyDatamatrixForAllGenes(): returns a list of tidy format data frames (tibbles) returns list of matrices (where number of genes in gene_names = number of items inside list) TRUE
#TEST: TidyDatamatrixForAllGenes(): for each tibble in list: number of rows of output tibble = nrow(data_mat) * ncol(data_mat)
#TEST: TidyDatamatrixForAllGenes(): column names are %in% c("ReadLen", "Pos", "Counts")
# gives:
# > str(tidy_all_genes_datamatrix)
# List of 5
# $ : Classes 'tbl_df', 'tbl' and 'data.frame':  45,961 observations of 3 variables:
#   $ ReadLen: int  10 11 12 13 14 15 16 17 18 19 ...
#   $ Pos    : int  1 1 1 1 1 1 1 1 1 1 ...
#   $ Counts : int  0 0 0 0 0 0 0 0 0 0 ...
# $ : Classes 'tbl_df', 'tbl' and 'data.frame':  99,589 observations of 3 variables:
#   $ ReadLen: int  10 11 12 13 14 15 16 17 18 19 ...
#   $ Pos    : int  1 1 1 1 1 1 1 1 1 1 ...
#   $ Counts : int  0 0 0 0 0 0 0 0 0 0 ...
# $ : Classes 'tbl_df', 'tbl' and 'data.frame':  69,085 observations of 3 variables:
#   $ ReadLen: int  10 11 12 13 14 15 16 17 18 19 ...
#   $ Pos    : int  1 1 1 1 1 1 1 1 1 1 ...
#   $ Counts : int  0 0 0 0 0 0 0 0 0 0 ...
# $ : Classes 'tbl_df', 'tbl' and 'data.frame':  143,869 observations of 3 variables:
#   $ ReadLen: int  10 11 12 13 14 15 16 17 18 19 ...
#   $ Pos    : int  1 1 1 1 1 1 1 1 1 1 ...
#   $ Counts : int  0 0 0 0 0 0 0 0 0 0 ...
# $ : Classes 'tbl_df', 'tbl' and 'data.frame':  82,123 observations of 3 variables:
#   $ ReadLen: int  10 11 12 13 14 15 16 17 18 19 ...
#   $ Pos    : int  1 1 1 1 1 1 1 1 1 1 ...
#   $ Counts : int  0 0 0 0 0 0 0 0 0 0 ...


print('Creating information tibble')


# The function CreateTranscriptInfoTibbleAllGenes applies all of the functions up until and including filtering for frame of interest to all of the genes 
# it returns one tibble, containing all of the genes and information
# to run for a single gene or to run each function separately, filter for a gene of interest and run

CreateTranscriptInfoTibbleAllGenes <- function(gene, dataset, hd_file, gff_df, filtering_frame){

          # Fetch the datamatrix for a single gene of interest , e.g YAL003W
          reads_pos_length <- GetGeneDatamatrix(gene = gene,
                                                dataset = dataset,
                                                hd_file = hd_file)

            # > str(reads_pos_length)
            # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...

          # Fetch the tidydatamatrix for a single gene of interest, e.g. YAL003W
          tidy_gene_datamatrix <- TidyDatamatrix(GetGeneDatamatrix(gene = gene,
                                                                   dataset = dataset,
                                                                   hd_file = hd_file),
                                                 startpos = 1,
                                                 startlen = 10)

            # > str(tidy_gene_datamatrix)
            # tibble [45,961 x 3] (S3: tbl_df/tbl/data.frame)
            # $ ReadLen: int [1:45961] 10 11 12 13 14 15 16 17 18 19 ...
            # $ Pos    : int [1:45961] 1 1 1 1 1 1 1 1 1 1 ...
            # $ Counts : int [1:45961] 0 0 0 0 0 0 0 0 0 0 ...
 


##### 



### Functions for A-site assignment of reads extracted from .h5 file ###
  

# #' #CalcAsiteFixedForAllGenes(): Calculate read A-site using a fixed displacement for read lengths for all genes in gene_names.
# #' 
# #' The function CalcAsiteFixed() from rscripts: read_count_functions is used within this function
# #' The assignment rules are specified in a user-supplied data frame, 'asite_displacement_length'.
# #' This function is used in TidyAsiteCountsByPositionAllGenes().
# #' 
# #' @param dataset name of dataset stored in .h5 file.
# #' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
# #' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
# #' @param asite_displacement_length data frame with columns `read_length` and `asite_displacement`
# #'  default: read_length = c(28, 29, 30), and asite_displacement = c(15, 15, 15).
# #' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
# #' 
# #' @return a list of numeric vectors if colsum_out = TRUE; matrices with a number of rows equivalent to number of rows in asite_displacement_length if colsum_out = FALSE
# #' 
# #' @examples
# #' 
# #' #### Does the function GetAllGeneDatamatrix need to be pre-defined here? Also CalcAsiteFixed?
# #' 
# #' CalcAsiteFixedForAllGenes(dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE)
# #' 
# #' @export
  # CalcAsiteFixedForAllGenes <- function(dataset, 
  #                                       hd_file, 
  #                                       min_read_length = 10, 
  #                                       asite_displacement_length = data.frame(
  #                                         read_length = c(28, 29, 30), 
  #                                         asite_displacement = c(15, 15, 15)
  #                                         ),
  #                                       colsum_out = TRUE){
  #  
  #   get_all_gene_datamatrix <- GetAllGeneDatamatrix(dataset = dataset, 
  #                                                   hd_file = hd_file)
  # 
  #   asite_counts_by_position_all_genes <- purrr::map(
  #     .x = get_all_gene_datamatrix,
  #     .f = CalcAsiteFixed,
  #     min_read_length = min_read_length,
  #     asite_displacement_length = asite_displacement_length,
  #     colsum_out = colsum_out
  #   )
  #   
  #   return(asite_counts_by_position_all_genes)
  # }
# TEST: CalcAsiteFixedForAllGenes(): if col_sum = TRUE, return list of numeric vectors returns list of matrices (where number of genes in gene_names = number of items inside list)? TRUE
# TEST: CalcAsiteFixedForAllGenes(): if col_sum = FALSE, return list of matrices? TRUE
# gives:
# > str(asite_counts_by_position_all_genes)
# List of 5
# $ : num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
# $ : num [1:2429] 0 0 0 0 0 0 0 0 0 0 ...
# $ : num [1:1685] 0 0 0 0 0 0 0 0 0 0 ...
# $ : num [1:3509] 0 0 0 0 0 0 0 0 0 0 ...
# $ : num [1:2003] 0 0 0 0 0 0 0 0 0 0 ...



# CalcAsiteFixed for single gene of interst
asite_counts_by_position <- CalcAsiteFixed(reads_pos_length,
                        min_read_length,
                        asite_displacement_length = data.frame(read_length = c(28, 29, 30),
                                                               asite_displacement = c(15, 15, 15)),
                        colsum_out = TRUE)

  # > str(asite_counts_by_position)
  # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...

          # This step extracts the asite reads from the reads_pos_length file
          # This step is incorporated in the funciton below 

      #' TidyAsiteCountsByPosition(): Aligns A-site assigned counts to the nucleotide positions for a single gene
      # Funtction with the aim of combining the asite count with the position of the nucleotides
      TidyAsiteCountsByPosition <- function(gene, dataset, hd_file, min_read_length, colsum_out = TRUE){

        asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))

        reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)

        asite_counts_by_position <- CalcAsiteFixed(reads_pos_length,
                                                   min_read_length = min_read_length,
                                                   asite_displacement_length = asite_displacement_length,
                                                   colsum_out = colsum_out)

        tidy_asite_counts <- tibble(Pos = 1:length(asite_counts_by_position), Count = asite_counts_by_position)

        return(tidy_asite_counts)

      }

      tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, dataset, hd_file, min_read_length, colsum_out = TRUE)
      
      # > str(tidy_asite_count_output)
      # tibble [1,121 x 2] (S3: tbl_df/tbl/data.frame)
      # $ Pos  : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
      # $ Count: num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...

      # The end result here is that the asite counts are aligned to the gene of interest
      # (in nucleotides, including UTRs and CDS)
      

#' TidyAsiteCountsByPositionAllGenes(): Align A-site assigned counts to the nucleotide positions for each gene contained within gene_names. 
#' 
#' CalcAsiteFixedForAllGenes() is used within this function.
#' 
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
#' 
#' @return a list of tibbles where the A-site assigned counts are aligned to the nucleotide positions of each gene of interest (UTRs and CDS)
#' 
#' @examples 
#' 
#' ### GetAllGeneDatamatrix and CalcAsiteFixedForAllGenes need to be pre-defined
#' 
#' TidyAsiteCountsByPositionAllGenes(dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE)
#' 
#' export
  # TidyAsiteCountsByPositionAllGenes <- function(dataset, hd_file, min_read_length = 10, colsum_out = TRUE){
  #   
  #   asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))
  #   
  #   asite_counts_by_position_all_genes <- CalcAsiteFixedForAllGenes(dataset = dataset, hd_file = hd_file, min_read_length = 10, colsum_out = TRUE)
  #   
  #   
  #   TidyAsiteCountsTibble <- function(asite_counts_by_position_all_genes){
  #     
  #     tidy_asite_counts_all_genes <- tibble(Pos = 1:length(asite_counts_by_position_all_genes), Count = asite_counts_by_position_all_genes)
  #     
  #   }
  #   
  #   tidy_asite_counts_all_genes <- purrr::map(
  #     .x = asite_counts_by_position_all_genes,
  #     .f = TidyAsiteCountsTibble
  #   )
  #   
  #   return(tidy_asite_counts_all_genes)
  # }
#TEST: TidyAsiteCountsByPositionAllGenes(): returns a list of tidy format data frame (tibble) where the number of items in the list = number of items in gene_names. 
#TEST: TidyAsiteCountsByPositionAllGenes(): each tibble has 2 columns.
#TEST: TidyAsiteCountsByPositionAllGenes(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df for each respective gene in gene_names.
#TEST: TidyAsiteCountsByPositionAllGenes(): the column names are %in% c("Pos", "Count")
# gives:
# > str(tidy_asite_counts_all_genes)
# List of 5
# Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 2 variables:
#   $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count: num  0 0 0 0 0 0 0 0 0 0 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   2429 observations of 2 variables:
#   $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count: num  0 0 0 0 0 0 0 0 0 0 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   1685 observations of 2 variables:
#   $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count: num  0 0 0 0 0 0 0 0 0 0 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   3509 observations of 2 variables:
#   $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count: num  0 0 0 0 0 0 0 0 0 0 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   2003 observations of 2 variables:
#   $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count: num  0 0 0 0 0 0 0 0 0 0 ...



#####


 
### Functions for assembling tibbles consisting of counts, positions and codons ###
  


#' TranscriptPosToCodonPos(): Assigns the nucleotide positions of UTRs and CDS to codon positions
#' 
#' This is a helper function for TranscriptPosToCodonPosForAllGenes().
#' 
#' @param gene from gene_names to pull information from the gff_df file
#' @param gff_df from which to extract the UTRs and CDS widths. 
#' 
#' @return Tidy data frame (tibble) containing the columns: `Gene`, `Pos` 
#' (position of nucleotides), `Pos_Codon1` (leading codon), `Pos_Codon2` (lagging codon) 
#' and `Frame` (reading frame). 
#' 
#' @example 
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' TranscriptPosToCodonPos(gene = "YAL003W", gff_df)
#' @export
TranscriptPosToCodonPos <- function(gene, gff_df){
    
    subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
    # the gff file is filtered down to the gene of interest
    
    UTR5_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR5") %>% select(width)
    # the width of the UTR5 is defined from the gff file, in this case 250 
    
    CDS_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(width)
    # the width of the CDS is defined from the gff file, in this case 621
    
    UTR3_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(width) 
    # the width of the UTR3 is defined from the gff file, in this case 250 
    
    transcript_length <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(end)
    # defines the length of the transcript including both UTRs and the CDS (in nucleotides), in this case 1121
    
    CDS_start <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start)
    # define the start position of the CDS from the gff file, in this case 251
    
    CDS_end <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end)
    # define the end position of the CDS from the gff file, in this case 871
    
    CDS_codon_positions <- rep(1:(CDS_width$width/3), each = 3)
    # In order to add the CDS codon positions the length of the CDS in nucleotides is 
    # divided by 3 (3 nucleotides = 1 codon) and repeated 3 times
      
      # > str(CDS_codon_positions)
      # int [1:621] 1 1 1 2 2 2 3 3 3 4 ...
      
    NA_UTR5_width <- rep(NA, times = UTR5_width$width)
    # NAs are repeated for the width of UTR5 
    
    NA_UTR3_width <- rep(NA, times = UTR3_width$width)
    # NAs are repeated for the width of UTR3
    
    transcript_gene_pos_poscodon_frame <- tibble(
      Gene = gene,
      Pos = 1:transcript_length$end,
      Pos_Codon1 = c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)),
      Pos_Codon2 = dplyr::lead(c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)), n = 3),
      Frame = seq(from = 2L, to = (transcript_length$end + 1L)) %% 3 # works as UTRS are 250, might not work with UTRS of others values
      # add the count column, likely something we want here 
    )
    
    # For codon pairs name of Pos_Codon became Pos_Codon1 and Pos_Codon2 was added 
    # where the line for Pos_Codon1 was copied and dplyr::lead((), n = 3) was added 
    
    return(transcript_gene_pos_poscodon_frame)
}
#TEST: TranscriptPosToCodonPos(): creates a tidy format data frame (tibble) = TRUE
#TEST: TranscriptPosToCodonPos(): the tibble contains 5 columns = TRUE
#TEST: TranscriptPosToCodonPos(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
#TEST: TranscriptPosToCodonPos(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame") 
# gives:
# > str(transcript_pos_to_codon_pos_output)
# Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 5 variables:
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...

transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene, gff_df)

# The end result is a table with the columns Gene, Pos, Pos_Codon and Frame (reading frame).
# For the Pos_Codon column the UTR positions have NA, while the CDS has the codon positions 

# #' #TranscriptPosToCodonPosForAllGenes(): Assigns the nucleotide positions of UTRs and CDS to codon positions for all genes in gene_names.
# #' 
# #' TranscriptPosToCodonPos() is used within TranscriptPosToCodonPosForAllGenes()
# #' 
# #' @param gene in gene_names to pull information from the gff_df file
# #' @param gff_df from which to extract the UTRs and CDS widths.
# #' 
# #' @return A Tidy data frame (tibble) containing the columns: `Gene`, `Pos`
# #' (position of nucleotides), `Pos_Codon1` (leading codon), `Pos_Codon2` (lagging codon)
# #' and `Frame` (reading frame) for each gene contained within gene_names.
# #' 
# #' @examples
# #' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
# #' gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
# #' 
# #' TranscriptPosToCodonPos <- function(gene, gff_df){
# #'    subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene)
# #'    UTR5_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR5") %>% select(width)
# #'    CDS_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(width)
# #'    UTR3_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(width)
# #'    transcript_length <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(end)
# #'    CDS_start <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start)
# #'    CDS_end <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end)
# #'    CDS_codon_positions <- rep(1:(CDS_width$width/3), each = 3)
# #'    NA_UTR5_width <- rep(NA, times = UTR5_width$width)
# #'    NA_UTR3_width <- rep(NA, times = UTR3_width$width)
# #'    transcript_gene_pos_poscodon_frame <- tibble(
# #'      Gene = gene,
# #'      Pos = 1:transcript_length$end,
# #'      Pos_Codon1 = c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)),
# #'      Pos_Codon2 = dplyr::lead(c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)), n = 3),
# #'      Frame = seq(from = 2L, to = (transcript_length$end + 1L)) %% 3 # works as UTRS are 250, might not work with UTRS of others values
# #'      )
# #'  return(transcript_gene_pos_poscodon_frame)
# #'  }
# #' 
# # TranscriptPosToCodonPosForAllGenes(gene_names, gff_df)

#  # export
#  # TranscriptPosToCodonPosForAllGenes <- function(.x = gene_names, gff_df){
#  #   
#  #   transcript_gene_pos_poscodon_frame_all_genes <- purrr::map(
#  #     .x = gene_names,
#  #     .f = TranscriptPosToCodonPos,
#  #     gff_df
#  #   )
#  # }
#TEST: TranscriptPosToCodonPosForAllGenes(): creates a  list of Tidy format data frames (tibbles), where the number of items in the list = number of items in gene_names = TRUE
#TEST: TranscriptPosToCodonPosForAllGenes(): each tibble contains 5 columns = TRUE
#TEST: TranscriptPosToCodonPosForAllGenes(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df for each respective gene.
#TEST: TranscriptPosToCodonPosForAllGenes(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame") 
# gives:
# List of 5
# Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 5 variables:
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   2429 observations of 5 variables:
#   $ Gene      : chr  "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   1685 observations of 5 variables:
#   $ Gene      : chr  "YAL012W" "YAL012W" "YAL012W" "YAL012W" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   3509 observations of 5 variables:
#   $ Gene      : chr  "YAL035W" "YAL035W" "YAL035W" "YAL035W" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   2003 observations of 5 variables:
#   $ Gene      : chr  "YAL038W" "YAL038W" "YAL038W" "YAL038W" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...

# transcript_gene_pos_poscodon_frame_all_genes <- TranscriptPosToCodonPosForAllGenes(gene_names, gff_df)

             

#' AddAsiteCountsToTranscriptPosToCodonPos(): merges A-site counts with transcript_pos_to_codon_pos_output
#' 
#' This is a helper function for AddAsiteCountsToTranscriptPosToCodonPosAllGenes 
#' TidyAsiteCountsByPosition() and TranscriptPosToCodonPos() are used in this function. 
#' 
#' @param gene in gene_names to pull information from the gff_df file 
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
#' @param gff_df from which to extract the UTRs and CDS widths.
#' 
#' @return a list of tibbles where the A-site assigned counts are aligned to the nucleotide positions of each gene of interest (UTRs and CDS)
#' 
#' @examples 
#' 
#' do TidyAsiteCountsByPosition() and TranscriptPosToCodonPos() need to be defined here?
#' 
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' AddAsiteCountsToTranscriptPosToCodonPos(gene = "YAL003W", dataset = "Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE, gff_df)
#' 
#' @export
AddAsiteCountsToTranscriptPosToCodonPos <- function(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
    
    tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE)
    
    transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene, gff_df)
    
    transcript_info_tibble <- dplyr::left_join(transcript_pos_to_codon_pos_output, tidy_asite_count_output, by = "Pos")
    
    return(transcript_info_tibble)
}
#TEST: AddAsiteCountsToTranscriptPosToCodonPos(): creates a tidy format data frame (tibble) = TRUE
#TEST: AddAsiteCountsToTranscriptPosToCodonPos(): the tibble contains 6 columns = TRUE
#TEST: AddAsiteCountsToTranscriptPosToCodonPos(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
#TEST: AddAsiteCountsToTranscriptPosToCodonPos(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame", "Count") 
# gives:
# > str(transcript_info_tibble)
# Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 6 variables:
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
#   $ Count     : num  0 0 0 0 0 0 0 0 0 0 ...

transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, dataset, hd_file, min_read_length, colsum_out = TRUE, gff_df)

# End result is a table which contains the asite counts, reading frame, codons and 
# codon positions for the gene of interest
  
          ##### ALL GENES #####        
  
          # AddAsiteCountsToTranscriptPosToCodonPos for all genes
              # AddAsiteCountsToTranscriptPosToCodonPosAllGenes <- function(gene_names, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
              # 
              #   transcript_info_tibble_all_genes <- purrr::map(
              #     .x = gene_names,
              #     .f = AddAsiteCountsToTranscriptPosToCodonPos,
              #     dataset, 
              #     hd_file, 
              #     min_read_length = 10, 
              #     colsum_out = TRUE,
              #     gff_df
              #   )
              #   
              # 
              # }
              # 
              # transcript_info_tibble_all_genes <- AddAsiteCountsToTranscriptPosToCodonPosAllGenes(gene_names, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
              # 
              # 
          # > str(transcript_info_tibble_all_genes)
          # List of 5
          # $ : tibble [1,121 x 6] (S3: tbl_df/tbl/data.frame)
          # ..$ Gene      : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
          # ..$ Pos       : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
          # ..$ Pos_Codon1: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Pos_Codon2: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Frame     : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...
          # ..$ Count     : num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
          # $ : tibble [2,429 x 6] (S3: tbl_df/tbl/data.frame)
          # ..$ Gene      : chr [1:2429] "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
          # ..$ Pos       : int [1:2429] 1 2 3 4 5 6 7 8 9 10 ...
          # ..$ Pos_Codon1: int [1:2429] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Pos_Codon2: int [1:2429] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Frame     : num [1:2429] 2 0 1 2 0 1 2 0 1 2 ...
          # ..$ Count     : num [1:2429] 0 0 0 0 0 0 0 0 0 0 ...
          # $ : tibble [1,685 x 6] (S3: tbl_df/tbl/data.frame)
          # ..$ Gene      : chr [1:1685] "YAL012W" "YAL012W" "YAL012W" "YAL012W" ...
          # ..$ Pos       : int [1:1685] 1 2 3 4 5 6 7 8 9 10 ...
          # ..$ Pos_Codon1: int [1:1685] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Pos_Codon2: int [1:1685] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Frame     : num [1:1685] 2 0 1 2 0 1 2 0 1 2 ...
          # ..$ Count     : num [1:1685] 0 0 0 0 0 0 0 0 0 0 ...
          # $ : tibble [3,509 x 6] (S3: tbl_df/tbl/data.frame)
          # ..$ Gene      : chr [1:3509] "YAL035W" "YAL035W" "YAL035W" "YAL035W" ...
          # ..$ Pos       : int [1:3509] 1 2 3 4 5 6 7 8 9 10 ...
          # ..$ Pos_Codon1: int [1:3509] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Pos_Codon2: int [1:3509] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Frame     : num [1:3509] 2 0 1 2 0 1 2 0 1 2 ...
          # ..$ Count     : num [1:3509] 0 0 0 0 0 0 0 0 0 0 ...
          # $ : tibble [2,003 x 6] (S3: tbl_df/tbl/data.frame)
          # ..$ Gene      : chr [1:2003] "YAL038W" "YAL038W" "YAL038W" "YAL038W" ...
          # ..$ Pos       : int [1:2003] 1 2 3 4 5 6 7 8 9 10 ...
          # ..$ Pos_Codon1: int [1:2003] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Pos_Codon2: int [1:2003] NA NA NA NA NA NA NA NA NA NA ...
          # ..$ Frame     : num [1:2003] 2 0 1 2 0 1 2 0 1 2 ...
          # ..$ Count     : num [1:2003] 0 0 0 0 0 0 0 0 0 0 ...
          
          #####
          
          
          
  AddCodonNamesToTranscriptInfoTibble <- function(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out = TRUE, gff_df){
    
    codon_table <- dplyr::filter(gene_poscodon_codon_i200, Gene == gene) %>% 
      dplyr::select("CodonPos_1", "CodonPos_2", "CodonPair")
    
    # this provides the positions of the codons and the codon names which can 
    # joined to the tibble generated by AddAsiteCountsToTranscriptPosToCodonPos 
    
    transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, dataset = dataset, hd_file = hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
    
    transcript_info_tibble <- left_join(transcript_info_tibble, codon_table, by = c("Pos_Codon1" = "CodonPos_1", "Pos_Codon2" = "CodonPos_2"), keep = FALSE)
    
    return(transcript_info_tibble)
      
  }
  
  transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
  
    # > str(transcript_info_tibble)
    # tibble [1,121 x 7] (S3: tbl_df/tbl/data.frame)
    # $ Gene      : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    # $ Pos       : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
    # $ Pos_Codon1: num [1:1121] NA NA NA NA NA NA NA NA NA NA ...
    # $ Pos_Codon2: num [1:1121] NA NA NA NA NA NA NA NA NA NA ...
    # $ Frame     : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...
    # $ Count     : num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
    # $ CodonPair : chr [1:1121] NA NA NA NA ...
  

            # ##### ALL GENES #####
            # 
            # AddCodonNamesToTranscriptInfoTibbleForAllGenes <- function(gene_poscodon_codon_i200, gene = gene_names, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
            #   
            #   codon_names_info_tibble_all_genes <- purrr::map(
            #     .x = gene_names,
            #     .f = AddCodonNamesToTranscriptInfoTibble,
            #     gene_poscodon_codon_i200,
            #     dataset, 
            #     hd_file, 
            #     min_read_length = 10, 
            #     colsum_out = TRUE,
            #     gff_df
            #   )
            #   
            # }
            # 
            # codon_names_info_tibble_all_genes <- AddCodonNamesToTranscriptInfoTibbleForAllGenes(gene_poscodon_codon_i200, gene = gene_names, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
            # 
            # #####

  #Filter for reading frame of interest, default is 0
  FilterForFrame <- function(transcript_info_tibble, filtering_frame = 0){
    
    transcript_info_tibble <- dplyr::filter(transcript_info_tibble, Frame == filtering_frame)
    
    return(transcript_info_tibble)
    
  }
  
  transcript_info_tibble <- FilterForFrame(transcript_info_tibble, filtering_frame)
  
    # > str(transcript_info_tibble)
    # tibble [374 x 7] (S3: tbl_df/tbl/data.frame)
    # $ Gene      : chr [1:374] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    # $ Pos       : int [1:374] 2 5 8 11 14 17 20 23 26 29 ...
    # $ Pos_Codon1: num [1:374] NA NA NA NA NA NA NA NA NA NA ...
    # $ Pos_Codon2: num [1:374] NA NA NA NA NA NA NA NA NA NA ...
    # $ Frame     : num [1:374] 0 0 0 0 0 0 0 0 0 0 ...
    # $ Count     : num [1:374] 0 0 0 0 0 0 0 0 0 0 ...
    # $ CodonPair : chr [1:374] NA NA NA NA ...

}

# CreateTranscriptInfoTibbleAllGenes (reads_pos_length, gene, dataset, hd_file, gff_df)
transcript_info_tibble <- suppressMessages(purrr::map_dfr(.x = gene_names, 
                                                          .f = CreateTranscriptInfoTibbleAllGenes, 
                                                          dataset, 
                                                          hd_file, 
                                                          gff_df, 
                                                          filtering_frame = 0))

# this takes all of the functions defined above and applies them to all of the genes in the sample. 

# For each gene, the first and last 250 nt may have NA in the codon column and Pos_Codon
# This is due to the 250 nt utr buffer, present in many riboviz example datasets 

# > str(transcript_info_tibble)
# tibble [1,870 x 7] (S3: tbl_df/tbl/data.frame)
# $ Gene      : chr [1:1870] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
# $ Pos       : int [1:1870] 2 5 8 11 14 17 20 23 26 29 ...
# $ Pos_Codon1: num [1:1870] NA NA NA NA NA NA NA NA NA NA ...
# $ Pos_Codon2: num [1:1870] NA NA NA NA NA NA NA NA NA NA ...
# $ Frame     : num [1:1870] 0 0 0 0 0 0 0 0 0 0 ...
# $ Count     : num [1:1870] 0 0 0 0 0 0 0 0 0 0 ...
# $ CodonPair : chr [1:1870] NA NA NA NA ..


# CalculateCountsPerCodonPos <- function(transcript_info_tibble){
#   
#   counts_per_codon_pos <- dplyr::group_by(transcript_info_tibble, Pos_Codon) %>% 
#     summarise(Count = sum(Count))
#   
#   
#   return(counts_per_codon_pos)
#   
#   
# }
#
# Because we have run the FilterForFrame function the CalculateCountsPerCodonPos 
# function gives you the counts per codon which are already present in Counts in 
# the filtered transcript_info_tibble

# test_function_2 <- function(gene, dataset, hd_file, gff_df, filtering_frame, codon_pair_of_interest){
  
  #transcript_info_tibble <- purrr::map_dfr(.x = gene_names, .f = test_function, dataset, hd_file, gff_df, filtering_frame)
  
  ###

print('Filtering for codon of interest')

  FilterForCodonPairOfInterestPositions <- function(transcript_info_tibble, codon_pair_of_interest, gene_name = gene_name){
  
      interesting_codon_table <- dplyr::filter(transcript_info_tibble, CodonPair == codon_pair_of_interest & Gene == gene_name)
      # Changed Codon == codon_of_interest to CodonPair = codon_of_interest
  
      interesting_first_codon_positions <- interesting_codon_table$Pos_Codon1
  
      # return(interesting_codon_positions)
      return(interesting_first_codon_positions)
  }
  
# interesting_first_codon_positions <- FilterForCodonPairOfInterestPositions(transcript_info_tibble, codon_pair_of_interest, gene_name)

  # > str(interesting_first_codon_positions)
  # num [1:8] 7 57 383 508 535 321 90 412
  
# interesting_first_codon_positions <- test_function_2(gene, dataset, hd_file, gff_df, filtering_frame = 0, codon_pair_of_interest = "TCC AAG")

# # Returns: 
# # > interesting_first_codon_position
# # [1]  7 57
# # 
# # Which refers to the position of the first codon in the codon pair as the second codon
# # in the pair is the first codon position + 1

  
### Expand frame around codon of interest 
  
ExpandCodonPairRegion <- function(.x = interesting_first_codon_positions, 
                                  transcript_info_tibble, 
                                  gene, 
                                  dataset, 
                                  hd_file, 
                                  expand_width = 5L, 
                                  remove_overhang = TRUE) {

  gene_poscodon_count <- na.omit(tibble(
        Gene = transcript_info_tibble$Gene,
        Pos_Codon1 = transcript_info_tibble$Pos_Codon1,
        Pos_Codon2 = transcript_info_tibble$Pos_Codon2,
        Rel_Count = transcript_info_tibble$Count
      ))

      # added na.omit() as the slice function goes by the row number instead of Pos_codon
      # na.omit() removes the UTRs so only CDS remains
  
  gene_length <- GetGeneLength(gene, dataset, hd_file)

  #if (remove_overhang) {
  # return an empty tibble if the desired region hangs over the edge of the coding region

  if (.x <= expand_width  |.x + expand_width > gene_length/3) {
    return()
  } else {
    output_codonpair_info <- tibble(
      dplyr::slice(gene_poscodon_count, (.x - expand_width):(.x + expand_width), each = FALSE),
      Rel_Pos =  seq(- expand_width, expand_width)
    )
    
    if(dim(output_codonpair_info)[1] == (2*expand_width + 1)){
      
    return(output_codonpair_info)
  }else{
    return()
    }
  }
}

# The if statement ensures that codon positions that are less/more than the 
# expand_width value are discarded 



# output_codonpair_info <- ExpandCodonPairRegion(.x = interesting_first_codon_positions, transcript_info_tibble, gene = gene_names, dataset, hd_file, expand_width = 5L, remove_overhang = TRUE)

# > output_codonpair_info
# # A tibble: 11 x 5
# Gene    Pos_Codon1 Pos_Codon2 Rel_Count Rel_Pos
# <chr>        <dbl>      <dbl>     <dbl>   <int>
#   1 YAL003W          2          3       429      -5
# 2 YAL003W          3          4       488      -4
# 3 YAL003W          4          5       102      -3
# 4 YAL003W          5          6       994      -2
# 5 YAL003W          6          7       146      -1
# 6 YAL003W          7          8       173       0
# 7 YAL003W          8          9       762       1
# 8 YAL003W          9         10        13       2
# 9 YAL003W         10         11       176       3
# 10 YAL003W         11         12        98       4
# 11 YAL003W         12         13       123       5




# Apply the ExpandCodonRegion function to the codons of interest to generate expanded tibbles for each position 

ExpandCodonPairRegionForList <- function(transcript_info_tibble, 
                                         gene, 
                                         dataset, 
                                         hd_file, 
                                         startpos, 
                                         startlen, 
                                         codon_pair_of_interest,
                                         gff_df, 
                                         expand_width){
  
  expand_codon_pair_region <- list()
  
  
  # A loop is used here so each unique gene can be processed separately. there were issues with using map and lapply 
  # on the filter function of FilterForCodonOfInterestPositions, and filtering using only one gene would produce NULLs
  
  for(gene_name in unique(gene)){
  
    interesting_codon_positions <- FilterForCodonPairOfInterestPositions(transcript_info_tibble, codon_pair_of_interest, gene_name)  
  
    transcript_info_tibble_gene <- dplyr::filter(transcript_info_tibble, Gene == gene_name)
    
    tmp_expand_codon_pair_region <- purrr::map(
      .x = interesting_codon_positions,
      .f = ExpandCodonPairRegion,
      transcript_info_tibble,
      gene = gene_name,
      dataset,
      hd_file,
      gff_df,
      expand_width = expand_width
    )
   
    # remove any tibbles from the list that are null, so are within one expand_width of the UTRs
    tmp_expand_codon_pair_region <- tmp_expand_codon_pair_region[!sapply(tmp_expand_codon_pair_region,is.null)]
    
    
    expand_codon_pair_region <- append(expand_codon_pair_region, tmp_expand_codon_pair_region)
    
  }
    
    return(expand_codon_pair_region)
   
}

expand_codon_pair_region <- ExpandCodonPairRegionForList(transcript_info_tibble, 
                                                         gene = unique(transcript_info_tibble$Gene), 
                                                         dataset, 
                                                         hd_file, 
                                                         startpos, 
                                                         startlen, 
                                                         codon_pair_of_interest, 
                                                         gff_df, 
                                                         expand_width)

# three fuctions called in one.returns a list containing a tibble for each occurrence of the codon of interest
# TEST:: expand_codon_region should be a list
# type(expand_codon_region)
# [1] "list"

# the dimensions of each item in the list shoud be [(2*expand_width+1)X4] as there are 4 rows; Gene, Pos_Codon, Rel_Count, Rel_Pos
# and one column for each position being included ie (-5 to 5) relative to the stop codon
# TEST:: Where an expand_width of =5 is given,the following would be observed:
# > dim(expand_codon_region[[1]])
# [1] 11  4

  # > str(expand_codon_pair_region)
  # List of 2
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene      : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon1: num [1:11] 2 3 4 5 6 7 8 9 10 11 ...
  # ..$ Pos_Codon2: num [1:11] 3 4 5 6 7 8 9 10 11 12 ...
  # ..$ Rel_Count : num [1:11] 429 488 102 994 146 173 762 13 176 98 ...
  # ..$ Rel_Pos   : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene      : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon1: num [1:11] 52 53 54 55 56 57 58 59 60 61 ...
  # ..$ Pos_Codon2: num [1:11] 53 54 55 56 57 58 59 60 61 62 ...
  # ..$ Rel_Count : num [1:11] 42 53 648 293 121 92 519 79 765 196 ...
  # ..$ Rel_Pos   : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...



### Normalisation

# Normalization carried out within each expanded frame so that they are comparable 
# Normalizes the ExpandCodonRegion list generating a RelCount column with the normalization values

print('Normalising data')

ExpandedCodonRegionNormalization <- function(.x, expand_width){
  
  # dplyr::mutate(.x, RelCount = PerCodonCounts / sum(PerCodonCounts) * (2 * expand_width + 1))
  dplyr::mutate(.x, RelCount = Rel_Count / sum(Rel_Count) * (2 * expand_width + 1))
  
}


# normalized_expanded_codonpair_region <- ExpandedCodonRegionNormalization(expand_codon_pair_region[[1]], expand_width = 5L)

  # > str(normalized_expanded_codonpair_region)
  # tibble [11 x 6] (S3: tbl_df/tbl/data.frame)
  # $ Gene      : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # $ Pos_Codon1: num [1:11] 2 3 4 5 6 7 8 9 10 11 ...
  # $ Pos_Codon2: num [1:11] 3 4 5 6 7 8 9 10 11 12 ...
  # $ Rel_Count : num [1:11] 429 488 102 994 146 173 762 13 176 98 ...
  # $ Rel_Pos   : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ RelCount  : num [1:11] 1.347 1.532 0.32 3.12 0.458 ...



# Normalization carried out for all the tibbles within ExpandList 
normalized_expand_list <- purrr::map(
  .x = expand_codon_pair_region,
  .f = ExpandedCodonRegionNormalization,
  expand_width = 5L
)

# TEST:: normalised_expand_list should be of type list
# type(normalized_expand_list)
# [1] "list"

# TEST:: Normalized_expand_list should be the same length as expand_codon_region
# > length(normalized_expand_list)==length(expand_codon_region)
# [1] TRUE

# TEST:: the dimensions of each item in the list shoud be [(2*expand_width+1) X 5] as there are now 5 rows; Gene, Pos_Codon, Rel_Count, Rel_Pos, RelCount
# > dim(expand_codon_region[[1]])
# [1] 11  5

# TEST:: At each position, the sum of RelCount should be equal to (2*expand_width+1)
# ie if the expand width was 5:
# sum(normalized_expand_list[[1]]$RelCount)
# [1] 11


### Overlaying the normalized expanded tibbles 

print('Calculating average')

# Function to overlay graphs into a single graph. Need to generate a single tibble 
# from NormalizedExpandList. Need to join by Rel_Pos, in RelCount need the mean for 
# each Rel_Pos (sum row(x) / number of row(x))

OverlayedTable <- function(normalized_expand_list, expand_width){
  
  number_of_objects <- length(normalized_expand_list)
  # The number of objects inside normalized_expand_list
  
  result = lapply(normalized_expand_list, "[", c("Rel_Pos", "RelCount"))
  # Reduces normalized_expand_list to the columns Rel_Pos and RelCount
  
  joined_result = result %>% purrr::reduce(full_join, by = c("Rel_Pos"), sum("RelCount"))
  
  joined_rows = joined_result %>% 
    mutate(SumRows = rowSums(select(joined_result, -"Rel_Pos")) / number_of_objects)
  
  overlayed_tibbles <- tibble::tibble(
    Rel_Pos = seq(- expand_width, expand_width),
    RelCount = joined_rows$SumRows
  )
}


Over <- OverlayedTable(normalized_expand_list, expand_width) 

  # > str(overlayed_tibbles)
  # tibble [11 x 2] (S3: tbl_df/tbl/data.frame)
  # $ Rel_Pos : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ RelCount: num [1:11] 0.751 0.863 1.351 2.099 0.452 ...

print('Creating plot')


overlayed_plot <- ggplot(Over, mapping = aes(x = Rel_Pos, y = RelCount)) + geom_line()

save_plot_pdf <- function(overlayed_plot, output_dir){
  overlayed_plot %>%
    ggsave(
      filename = file.path(output_dir,"inhibitory_codon_plot.pdf"),
      width = 6, height = 5
    )
}

save_plot_pdf(overlayed_plot, output_dir)

print('Done')




