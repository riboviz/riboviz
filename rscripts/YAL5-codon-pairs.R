
# YAL5_h5 is at location $HOME/riboviz/riboviz/Mok-simYAL5/output/A/A.h5

# source packages and functions from rscripts 
source(here::here("rscripts", "read_count_functions.R"))
source(here::here("rscripts", "stats_figs_block_functions.R"))


# other packages I loaded that I used 
library(tidyverse)
# install.packages("zoo") - new package required for function rollapply
library(zoo)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("plotly")
library(plotly)
library(purrr)
library(dplyr)

# Given an h5 file, GFF file and .tsv file, this script creates a metafeature plot 
# for the codon/codon pair of interest. 


# FOR CODON PAIRS: USE gene_poscodon_codon_i200 THROUGHOUT (derived from yeast_codon_pos_i200)


# Load the necessary input files 

YAL5_h5 <- here::here("Mok-simYAL5", "output", "A", "A.h5")
# The h5 file from the dataset of interest 
hd_file <- YAL5_h5
# assign YAL5_h5 file to a general name applied throughput the script 
# the script does not have to change based on different hd_files being used

gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
# The GFF file for the simulated dataset, given the general name gff_df so that
# the script does not have to change based on different gff_df files being used

dataset = "Mok-simYAL5"

gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
# names of the genes used in the Mok-simYAL5 dataset


# Import the .tsv file: 
yeast_codon_pos_i200 <- readr::read_tsv(file = here::here("data", "yeast_codon_table.tsv"))

  # Confirm that the .tsv file has been imported and check structure:
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

#' GetAllGeneDatamatrix(): Get matrix of read counts by length for all genes contained within the gene_names file from .h5 file.
#'
#' This accesses the attribute `reads/data` in .h5 file for each gene listed in gene_names. 
#' 
#' @param list of each gene name to pull out read counts for.
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' 
#' @ return a list of numeric matrices of read count data for each given gene in given dataset.
#'
#' @examples
#' GetAllGeneDatamatrix(dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"))
#'
#' @export
GetAllGeneDatamatrix <- function(dataset, hd_file){
          
  gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
          
  get_all_gene_datamatrix <- purrr::map(
      .x = gene_names,
      .f = GetGeneDatamatrix,
      dataset,
      hd_file
      )
          
    return(get_all_gene_datamatrix)
}
#TEST: GetAllGeneDatamatrix(): returns list of matrices (where number of genes in gene_names = number of items inside list) TRUE
# gives:
# > str(get_all_gene_datamatrix)
# List of 5
#   $ : int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#   $ : int [1:41, 1:2429] 0 0 0 0 0 0 0 0 0 0 ...
#   $ : int [1:41, 1:1685] 0 0 0 0 0 0 0 0 0 0 ...
#   $ : int [1:41, 1:3509] 0 0 0 0 0 0 0 0 0 0 ...
#   $ : int [1:41, 1:2003] 0 0 0 0 0 0 0 0 0 0 ...


get_all_gene_datamatrix <- GetAllGeneDatamatrix(dataset = "Mok-simYAL5", hd_file = YAL5_h5)          
        

  
#' TidyDatamatrixForAllGenes(): Convert gene data matrices for all genes into tidy dataframes
#' 
#' Converts each data matrix into readable tidy format with columns: `ReadLen`, `Pos`, `Counts`
#' to hold read lengths, position (in transcript-centric coordinates), and number of reads for all genes
#' 
#' @param list of each gene name to pull out read counts for.
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param startpos numeric value, start position along transcript-centric alignment for specific gene, default = 1.
#' @param startlen numeric value, value from which to begin counting lengths from (ie equivalent to '10' in read_lengths= 10:50), default = 10.
#' 
#' @return a list of tidy format data frames (tibbles), with columns: `ReadLen`, `Pos` and `Counts`
#' 
#' @examples 
#' TidyDatamatrixForAllGenes(dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), startpos = 1, startlen = 10)
#' 
#' @export
TidyDatamatrixForAllGenes <- function(dataset, hd_file, startpos = 1, startlen = 10){
          
    gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
          
    .x = get_all_gene_datamatrix <- GetAllGeneDatamatrix(dataset = dataset, 
                                                         hd_file = hd_file)
          
    tidy_all_genes_datamatrix <- purrr::map(
        .x = get_all_gene_datamatrix,
        .f = TidyDatamatrix,
        startpos = 1,
        startlen = 10
        )
    
    return(tidy_all_genes_datamatrix)
}
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
          

tidy_all_genes_datamatrix <- TidyDatamatrixForAllGenes(hd_file, dataset, startpos = 1, startlen = 10)          


          
          
# # test_function <- function(gene, dataset, hd_file, gff_df, filtering_frame){
# 
#   # Fetch the datamatrix for a single gene of interest , e.g YAL003W
#   reads_pos_length <- GetGeneDatamatrix(gene = "YAL003W", 
#                                         dataset  = dataset, 
#                                         hd_file)
#   
#     # > str(reads_pos_length)
#     # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
#   
#   # Fetch the tidydatamatrix for a single gene of interest, e.g. YAL003W 
#   tidy_gene_datamatrix <- TidyDatamatrix(GetGeneDatamatrix(gene = "YAL003W", 
#                                                            dataset = dataset, 
#                                                            hd_file), 
#                                          startpos = 1, 
#                                          startlen = 10)
#   
#     # > str(tidy_gene_datamatrix)
#     # tibble [45,961 x 3] (S3: tbl_df/tbl/data.frame)
#     # $ ReadLen: int [1:45961] 10 11 12 13 14 15 16 17 18 19 ...
#     # $ Pos    : int [1:45961] 1 1 1 1 1 1 1 1 1 1 ...
#     # $ Counts : int [1:45961] 0 0 0 0 0 0 0 0 0 0 ...
 


##### 

### Functions for A-site assignment of reads extracted from .h5 file ###
  

#' CalcAsiteFixedForAllGenes(): Calculate read A-site using a fixed displacement for read lengths for all genes in gene_names.
#' 
#' The assignment rules are specified in a user-supplied data frame, 'asite_displacement_length'.
#' 
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param asite_displacement_length data frame with columns `read_length` and `asite_displacement`
#'  default: read_length = c(28, 29, 30), and asite_displacement = c(15, 15, 15).
#' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
#' 
#' @return a list of numeric vectors if colsum_out = TRUE; matrices with a number of rows equivalent to number of rows in asite_displacement_length if colsum_out = FALSE
#' 
#' @examples 
#' CalcAsiteFixedForAllGenes(dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE)
#' 
#' @export
CalcAsiteFixedForAllGenes <- function(dataset, 
                                      hd_file, 
                                      min_read_length = 10, 
                                      asite_displacement_length = data.frame(
                                        read_length = c(28, 29, 30), 
                                        asite_displacement = c(15, 15, 15)
                                        ),
                                      colsum_out = TRUE){
 
  get_all_gene_datamatrix <- GetAllGeneDatamatrix(dataset = dataset, 
                                                  hd_file = hd_file)

  asite_counts_by_position_all_genes <- purrr::map(
    .x = get_all_gene_datamatrix,
    .f = CalcAsiteFixed,
    min_read_length = min_read_length,
    asite_displacement_length = asite_displacement_length,
    colsum_out = colsum_out
  )
  
  return(asite_counts_by_position_all_genes)
}
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


asite_counts_by_position_all_genes <- CalcAsiteFixedForAllGenes(dataset = "Mok-simYAL5", 
                                                                hd_file, 
                                                                min_read_length = 10, 
                                                                colsum_out = TRUE)



      # # CalcAsiteFixed for single gene of interst 
      # asite_counts_by_position <- CalcAsiteFixed(reads_pos_length, 
      #                         min_read_length = 10, 
      #                         asite_displacement_length = data.frame(read_length = c(28, 29, 30), 
      #                                                                asite_displacement = c(15, 15, 15)), 
      #                         colsum_out = TRUE)
      # 
      #   # > str(asite_counts_by_position)
      #   # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...


  
  
  # Funtction with the aim of combining the asite count with the position of the nucleotides           
  TidyAsiteCountsByPosition <- function(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE){
    
    asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))
    
    reads_pos_length <- GetGeneDatamatrix(gene, dataset = dataset, hd_file)
    
    asite_counts_by_position <- CalcAsiteFixed(reads_pos_length, 
                                               min_read_length = min_read_length, 
                                               asite_displacement_length = asite_displacement_length,
                                               colsum_out = colsum_out)
    
    tidy_asite_counts <- tibble(Pos = 1:length(asite_counts_by_position), Count = asite_counts_by_position)
    
    return(tidy_asite_counts)
    
  }
  
  tidy_asite_count_output <- TidyAsiteCountsByPosition(gene = "YAL003W", dataset, hd_file, min_read_length = 10, colsum_out = TRUE)
  
    # > str(tidy_asite_count_output)
    # tibble [1,121 x 2] (S3: tbl_df/tbl/data.frame)
    # $ Pos  : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
    # $ Count: num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
  
  
            ##### ALL GENES #####
            # TidyAsiteCountsByPosition function applied to all genes of interest
            TidyAsiteCountsByPositionAllGenes <- function(dataset, hd_file, min_read_length = 10, colsum_out = TRUE){

              asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))

              asite_counts_by_position_all_genes <- CalcAsiteFixedForAllGenes(dataset = dataset, hd_file = hd_file, min_read_length = 10, colsum_out = TRUE)


              TidyAsiteCountsTibble <- function(asite_counts_by_position_all_genes){

                tidy_asite_counts_all_genes <- tibble(Pos = 1:length(asite_counts_by_position_all_genes), Count = asite_counts_by_position_all_genes)

              }


              tidy_asite_counts_all_genes <- purrr::map(
                .x = asite_counts_by_position_all_genes,
                .f = TidyAsiteCountsTibble
              )

              return(tidy_asite_counts_all_genes)

            }

            tidy_asite_counts_all_genes <- TidyAsiteCountsByPositionAllGenes(dataset, hd_file, min_read_length = 10, colsum_out = TRUE)

            # > str(tidy_asite_counts_all_genes)
            # List of 5
            # $ : tibble [1,121 x 2] (S3: tbl_df/tbl/data.frame)
            # ..$ Pos  : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
            # ..$ Count: num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
            # $ : tibble [2,429 x 2] (S3: tbl_df/tbl/data.frame)
            # ..$ Pos  : int [1:2429] 1 2 3 4 5 6 7 8 9 10 ...
            # ..$ Count: num [1:2429] 0 0 0 0 0 0 0 0 0 0 ...
            # $ : tibble [1,685 x 2] (S3: tbl_df/tbl/data.frame)
            # ..$ Pos  : int [1:1685] 1 2 3 4 5 6 7 8 9 10 ...
            # ..$ Count: num [1:1685] 0 0 0 0 0 0 0 0 0 0 ...
            # $ : tibble [3,509 x 2] (S3: tbl_df/tbl/data.frame)
            # ..$ Pos  : int [1:3509] 1 2 3 4 5 6 7 8 9 10 ...
            # ..$ Count: num [1:3509] 0 0 0 0 0 0 0 0 0 0 ...
            # $ : tibble [2,003 x 2] (S3: tbl_df/tbl/data.frame)
            # ..$ Pos  : int [1:2003] 1 2 3 4 5 6 7 8 9 10 ...
            # ..$ Count: num [1:2003] 0 0 0 0 0 0 0 0 0 0 ...

            #####

  
  # TranscriptPosToCodonPos generates a table with the columns Gene, Pos (position 
  # of the nuceotide), Pos_Codon1 (leading codon), Pos_Codon2 (lagging codon) and
  # frame (reading frame).
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
  
  
  transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene = "YAL003W", gff_df)
    
    # > str(transcript_pos_to_codon_pos_output)
    # tibble [1,121 x 5] (S3: tbl_df/tbl/data.frame)
    # $ Gene      : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    # $ Pos       : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
    # $ Pos_Codon1: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
    # $ Pos_Codon2: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
    # $ Frame     : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...
  
  
              ##### ALL GENES #####
              # TranscriptPosToCodonPos applied to all genes, generates a list of tibbles
              # with the same structure as described above
              TranscriptPosToCodonPosForAllGenes <- function(.x = gene_names, gff_df){

                transcript_gene_pos_poscodon_frame_all_genes <- purrr::map(
                  .x = gene_names,
                  .f = TranscriptPosToCodonPos,
                  gff_df
                )

              }

              transcript_gene_pos_poscodon_frame_all_genes <- TranscriptPosToCodonPosForAllGenes(gene_names, gff_df)

              # List of 5
              # $ : tibble [1,121 x 5] (S3: tbl_df/tbl/data.frame)
              # ..$ Gene      : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
              # ..$ Pos       : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
              # ..$ Pos_Codon1: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Pos_Codon2: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Frame     : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...
              # $ : tibble [2,429 x 5] (S3: tbl_df/tbl/data.frame)
              # ..$ Gene      : chr [1:2429] "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
              # ..$ Pos       : int [1:2429] 1 2 3 4 5 6 7 8 9 10 ...
              # ..$ Pos_Codon1: int [1:2429] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Pos_Codon2: int [1:2429] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Frame     : num [1:2429] 2 0 1 2 0 1 2 0 1 2 ...
              # $ : tibble [1,685 x 5] (S3: tbl_df/tbl/data.frame)
              # ..$ Gene      : chr [1:1685] "YAL012W" "YAL012W" "YAL012W" "YAL012W" ...
              # ..$ Pos       : int [1:1685] 1 2 3 4 5 6 7 8 9 10 ...
              # ..$ Pos_Codon1: int [1:1685] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Pos_Codon2: int [1:1685] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Frame     : num [1:1685] 2 0 1 2 0 1 2 0 1 2 ...
              # $ : tibble [3,509 x 5] (S3: tbl_df/tbl/data.frame)
              # ..$ Gene      : chr [1:3509] "YAL035W" "YAL035W" "YAL035W" "YAL035W" ...
              # ..$ Pos       : int [1:3509] 1 2 3 4 5 6 7 8 9 10 ...
              # ..$ Pos_Codon1: int [1:3509] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Pos_Codon2: int [1:3509] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Frame     : num [1:3509] 2 0 1 2 0 1 2 0 1 2 ...
              # $ : tibble [2,003 x 5] (S3: tbl_df/tbl/data.frame)
              # ..$ Gene      : chr [1:2003] "YAL038W" "YAL038W" "YAL038W" "YAL038W" ...
              # ..$ Pos       : int [1:2003] 1 2 3 4 5 6 7 8 9 10 ...
              # ..$ Pos_Codon1: int [1:2003] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Pos_Codon2: int [1:2003] NA NA NA NA NA NA NA NA NA NA ...
              # ..$ Frame     : num [1:2003] 2 0 1 2 0 1 2 0 1 2 ...

              #####

              
  # create a codon_per_codons_count_table - merge asite counts with transcript_info_tibble  
  AddAsiteCountsToTranscriptPosToCodonPos <- function(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
    
    tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE)
    
    transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene, gff_df)
    
    transcript_info_tibble <- dplyr::left_join(transcript_pos_to_codon_pos_output, tidy_asite_count_output, by = "Pos")
    
    return(transcript_info_tibble)
    
  }
   
  transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene = "YAL003W", dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
  
    # > str(transcript_info_tibble)
    # tibble [1,121 x 6] (S3: tbl_df/tbl/data.frame)
    # $ Gene      : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    # $ Pos       : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
    # $ Pos_Codon1: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
    # $ Pos_Codon2: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
    # $ Frame     : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...
    # $ Count     : num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
  
  
          ##### ALL GENES #####        
  
          # AddAsiteCountsToTranscriptPosToCodonPos for all genes
          AddAsiteCountsToTranscriptPosToCodonPosAllGenes <- function(gene_names, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){

            transcript_info_tibble_all_genes <- purrr::map(
              .x = gene_names,
              .f = AddAsiteCountsToTranscriptPosToCodonPos,
              dataset, 
              hd_file, 
              min_read_length = 10, 
              colsum_out = TRUE,
              gff_df
            )
            

          }

          transcript_info_tibble_all_genes <- AddAsiteCountsToTranscriptPosToCodonPosAllGenes(gene_names, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)

  
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
          
          
          
  AddCodonNamesToTranscriptInfoTibble <- function(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
    
    codon_table <- dplyr::filter(gene_poscodon_codon_i200, Gene == gene) %>% 
      dplyr::select("CodonPos_1", "CodonPos_2", "CodonPair")
    
    transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, dataset = dataset, hd_file = hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
    
    transcript_info_tibble_join <- left_join(transcript_info_tibble, codon_table, by = c("Pos_Codon1" = "CodonPos_1", "Pos_Codon2" = "CodonPos_2"), keep = FALSE)
    
    return(transcript_info_tibble_join)
      
  }
  
  transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(gene_poscodon_codon_i200, gene = "YAL003W", dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
  
    # > str(transcript_info_tibble)
    # tibble [1,121 x 7] (S3: tbl_df/tbl/data.frame)
    # $ Gene      : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    # $ Pos       : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
    # $ Pos_Codon1: num [1:1121] NA NA NA NA NA NA NA NA NA NA ...
    # $ Pos_Codon2: num [1:1121] NA NA NA NA NA NA NA NA NA NA ...
    # $ Frame     : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...
    # $ Count     : num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
    # $ CodonPair : chr [1:1121] NA NA NA NA ...
  

          ##### ALL GENES #####
          
          AddCodonNamesToTranscriptInfoTibbleForAllGenes <- function(gene_poscodon_codon_i200, gene = gene_names, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
            
            codon_names_info_tibble_all_genes <- purrr::map(
              .x = gene_names,
              .f = AddCodonNamesToTranscriptInfoTibble,
              gene_poscodon_codon_i200,
              dataset, 
              hd_file, 
              min_read_length = 10, 
              colsum_out = TRUE,
              gff_df
            )
            
          }
  
          codon_names_info_tibble_all_genes <- AddCodonNamesToTranscriptInfoTibbleForAllGenes(gene_poscodon_codon_i200, gene = gene_names, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
  
          #####

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

# }

# test_function(reads_pos_length, gene, dataset, hd_file, gff_df)
# transcript_info_tibble <- purrr::map_dfr(.x = gene_names, .f = test_function, dataset, hd_file, gff_df, filtering_frame = 0)

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

  FilterForCodonPairOfInterestPositions <- function(transcript_info_tibble, codon_pair_of_interest){
  
      interesting_codon_table <- dplyr::filter(transcript_info_tibble, CodonPair == "TCC AAG")
      # Changed Codon == codon_of_interest to CodonPair = codon_of_interest
  
      interesting_first_codon_positions <- interesting_codon_table$Pos_Codon1
  
      # return(interesting_codon_positions)
      return(interesting_first_codon_positions)
  }
  
  interesting_first_codon_positions <- FilterForCodonPairOfInterestPositions(transcript_info_tibble, codon_pair_of_interest)
# }

# interesting_first_codon_positions <- test_function_2(gene, dataset, hd_file, gff_df, filtering_frame = 0, codon_pair_of_interest = "TCC AAG")

# # Returns: 
# # > interesting_first_codon_position
# # [1]  7 57
# # 
# # Which refers to the position of the first codon in the codon pair as the second codon
# # in the pair is the first codon position + 1


ExpandCodonPairRegion <- function(.x = interesting_first_codon_positions, transcript_info_tibble, gene, dataset, hd_file, startpos = 1, startlen = 10, gff_df, expand_width = 5L, remove_overhang = TRUE) {

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

  if (.x < expand_width | .x + expand_width > gene_length ) {
    return( tibble() )
  } else {
    output_codonpair_info <- tibble(
      dplyr::slice(gene_poscodon_count, (.x - expand_width):(.x + expand_width), each = FALSE),
      Rel_Pos =  seq(- expand_width, expand_width)
    )
    return(output_codonpair_info)
  }
  #  } # if(remove_overhang)
}


output_codonpair_info <- ExpandCodonPairRegion(.x = interesting_first_codon_positions, transcript_info_tibble, gene = gene_names, dataset, hd_file, expand_width = 5L, remove_overhang = TRUE)

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

ExpandCodonPairRegionForList <- function(transcript_info_tibble, gene, dataset, hd_file, startpos = 1, startlen = 10, codon_pair_of_interest, gff_df, expand_width = 5L){
  
  interesting_codon_positions <- FilterForCodonPairOfInterestPositions(transcript_info_tibble, codon_pair_of_interest)  
  
  expand_codon_pair_region <- purrr::map(
    .x = interesting_codon_positions,
    .f = ExpandCodonPairRegion,
    transcript_info_tibble,
    gene,
    dataset,
    hd_file,
    startpos = 1,
    startlen = 10,
    gff_df,
    expand_width = 5L
  )
 
  return(expand_codon_pair_region)
   
}

expand_codon_pair_region <- ExpandCodonPairRegionForList(transcript_info_tibble, gene = gene_names, dataset, hd_file, startpos = 1, startlen = 10, codon_pair_of_interest = "TCC AAG", gff_df, expand_width = 5L)

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


# Normalization carried out within each expanded frame so that they are comparable 
# Normalizes the ExpandCodonRegion list generating a RelCount column with the normalization values
ExpandedCodonRegionNormalization <- function(.x, expand_width = 5L){
  
  # dplyr::mutate(.x, RelCount = PerCodonCounts / sum(PerCodonCounts) * (2 * expand_width + 1))
  dplyr::mutate(.x, RelCount = Rel_Count / sum(Rel_Count) * (2 * expand_width + 1))
  
}


normalized_expanded_codonpair_region <- ExpandedCodonRegionNormalization(expand_codon_pair_region[[1]], expand_width = 5L)

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



# Function to overlay tibbles to generate a single graph
OverlayedTable <- function(normalized_expand_list, expand_width = 5L){
  number_of_objects <- length(normalized_expand_list)
  
  result = lapply(normalized_expand_list, "[", c("Rel_Pos", "RelCount"))
  
  joined_result = result %>% purrr::reduce(full_join, by = c("Rel_Pos"), sum("RelCount"))
  
  joined_rows = joined_result %>% 
    mutate(SumRows = rowSums(select(joined_result, -"Rel_Pos")) / number_of_objects)
  
  overlayed_tibbles <- tibble::tibble(
    Rel_Pos = seq(- expand_width, expand_width),
    RelCount = joined_rows$SumRows
  )
}


overlayed_tibbles <- OverlayedTable(normalized_expand_list, expand_width = 5L)  

  # > str(overlayed_tibbles)
  # tibble [11 x 2] (S3: tbl_df/tbl/data.frame)
  # $ Rel_Pos : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ RelCount: num [1:11] 0.751 0.863 1.351 2.099 0.452 ...


overlayed_plot <- ggplot(Over, mapping = aes(x = Rel_Pos, y = RelCount)) + geom_line()



# how to apply to all genes?




