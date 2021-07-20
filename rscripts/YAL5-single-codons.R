# YAL5_h5 is at location $HOME/riboviz/riboviz/Mok-simYAL5/output/A/A.h5

# Given an h5 file, GFF file and .tsv file, this script creates a metafeature plot 
# for the codon of interest

print('Starting process')

# source packages and functions from rscripts 
source(here::here("rscripts", "read_count_functions.R"))
source(here::here("rscripts", "stats_figs_block_functions.R"))

suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument('-i', '--input', help='Path input to h5 file')
parser$add_argument('-d', '--dataset', help='Name of the dataset being studied')
parser$add_argument('-g', '--gff', help='Path t the GFF3 file of the organism being studied')
parser$add_argument('-a', '--annotation', help='Path to codon table for organism')
parser$add_argument('--codon', help='Codon of interest')
parser$add_argument('-o', '--output', help='Path to output directory')
parser$add_argument('--expand_width', help='the desired range either side of the codon of interest', default = 5)
parser$add_argument('--startpos', help='position of the start codon', default = 1)
parser$add_argument('--startlen', help='smallest length of reads', default = 10) # please correct if wrong
parser$add_argument('--frame', help='frame to be studied', default = 0)
parser$add_argument('--minreadlen', help='minimum read length', default = 10)
parser$add_argument('--colsum_out', help='logical', default = TRUE)

args <- parser$parse_args()

hd_file <- here::here("Mok-simYAL5", "output", "A", "A.h5")
dataset <- "Mok-simYAL5"
codon_of_interest <- 'TCC'
expand_width = 5L
startpos <-1
startlen <- 10
filtering_frame <- 0
min_read_length <- 10
yeast_codon_table <- here::here("data", "yeast_codon_table.tsv")
gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3")
colsum_out <- TRUE

hd_file <- args$input
dataset <- args$dataset
gff <- args$gff
yeast_codon_table <- args$annotation
codon_of_interest <- args$codon
output_dir <- args$output
expand_width <- args$expand_width
startpos <- args$startpos
startlen <- args$startlen
filtering_frame <- args$frame
min_read_length <- args$minreadlen
colsum_out <- args$colsum_out

# 

gff_df <- readGFFAsDf(gff)
gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name

# The GFF file for the simulated dataset, given the general name gff_df so that 
# the script does not have to change based on different gff_df files being used 

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




# Filter down the yeast_codon_pos_i200 file to the gene that you are 
# working with, in this case YAL003W
# YAL003W_pos <- dplyr::filter(yeast_codon_pos_i200, Gene=="YAL003W")


print('Creating information tibble')


# The function CreateTranscriptInfoTibbleAllGenes applies all of the functions up until filtering for frame of interest to all of the genes 
# it returns one tibble, containing all of the genes and information
# to run for a single gene or to run each function separately, filter for a gene of interest and run

# YAL003W_pos <- dplyr::filter(gene_poscodon_codon_i200, Gene=="YAL003W")
# gene <- 'YAL003W'

CreateTranscriptInfoTibbleAllGenes<- function(gene, dataset, hd_file, gff_df, colsum_out, startpos, startlen){
    
    
    ### Functions for A-site assignment of reads extracted from .h5 file ###
    
    
    #' TidyAsiteCountsByPosition(): Align A-site assigned counts to the nucleotide positions for a single gene
    #' 
    #' The function ReadAsiteDisplacementLengthFromFile() from rscripts: is used within this function.
    #' The function GetGeneDatamatrix() from rscripts: is used within this function. 
    #' The function CalcAsiteFixed() from rscripts: read_count_functions is used within this function.
    #' 
    #' @param dataset name of dataset stored in .h5 file.
    #' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
    #' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
    #' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
    #' 
    #' @return a tidy format data frame (tibble) with the columns "Pos" and "Count"
    #' 
    #' @examples
    #' 
    #' TidyAsiteCountsByPosition(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE)
    #' 
    #' @export      
    TidyAsiteCountsByPosition <- function(gene, dataset, hd_file, min_read_length, colsum_out = TRUE){
      
      asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", 
                                                                                  "yeast_standard_asite_disp_length.txt"))
      
      reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
      
      asite_counts_by_position <- CalcAsiteFixed(reads_pos_length, 
                                                 min_read_length = min_read_length, 
                                                 asite_displacement_length = asite_displacement_length,
                                                 colsum_out = colsum_out)
      
      tidy_asite_counts <- tibble(Pos = 1:length(asite_counts_by_position), 
                                  Count = asite_counts_by_position)
      
      return(tidy_asite_counts)
      
    }
    #TEST: TidyAsiteCountsByPosition(): returns a tidy format data frame (tibble). 
    #TEST: TidyAsiteCountsByPosition(): the tibble has 2 columns.
    #TEST: TidyAsiteCountsByPosition(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df for the gene.
    #TEST: TidyAsiteCountsByPosition(): the column names are %in% c("Pos", "Count")
    # gives:          
    # > str(tidy_asite_count_output)
    # Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 2 variables:
    # $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
    # $ Count: num  0 0 0 0 0 0 0 0 0 0 ...
    
    tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, 
                                                         dataset, 
                                                         hd_file, 
                                                         min_read_length, 
                                                         colsum_out)
    # The end result here is that the A-site assigned counts are aligned to the gene of interest
    # (in nucleotides, including UTRs and UTRs)
    
   
    
     #####
    
    
    
    ### Functions for assembling tibbles consisting of counts, positions and codons ###
    

    
    #' TranscriptPosToCodonPos(): Assigns the nucleotide positions of UTRs and CDS to codon positions
    #' 
    #' Generates a tidy format data frame (tibble) with the columns Gene, Pos (position 
    #' of the nucleotide), Pos_Codon and frame (reading frame).
    #' 
    #' @param gene from gene_names to pull information from the gff_df file
    #' @param gff_df from which to extract the UTRs and CDS widths. 
    #' 
    #' @return Tidy data frame (tibble) containing the columns: "Gene", "Pos" 
    #' (position of nucleotides), "Pos_Codon" (leading codon), "Pos_Codon2" (lagging codon) 
    #' and "Frame" (reading frame). 
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
        Pos_Codon = c(rep(NA, times = UTR5_width$width), 
                      CDS_codon_positions, 
                      rep(NA, times = UTR3_width$width)),
        Frame = seq(from = 2L, to = (transcript_length$end + 1L)) %% 3 # works as UTRS are 250, might not work with UTRS of others values
        # add the count column, likely something we want here 
      )
      
      return(transcript_gene_pos_poscodon_frame)
    }
    #TEST: TranscriptPosToCodonPos(): creates a tibble = TRUE
    #TEST: TranscriptPosToCodonPos(): the tibble contains 4 columns = TRUE
    #TEST: TranscriptPosToCodonPos(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
    #TEST: TranscriptPosToCodonPos(): the column names are %in% c("Gene", "Pos", "Pos_Codon", "Frame" (reading frame)) 
    #TEST: TranscriptPosToCodonPos(): The first and last 250 rows of the column "Pos_codon" for the gene contains NA = TRUE. 
    # gives:
    # > str(transcript_pos_to_codon_pos_output)
    # Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 4 variables:
    #   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    #   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
    #   $ Pos_Codon : int  NA NA NA NA NA NA NA NA NA NA ...
    #   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
    
    transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene, gff_df)
    
    
    
    
    #' AddAsiteCountsToTranscriptPosToCodonPos(): merges A-site counts with transcript_pos_to_codon_pos_output
    #' 
    #' This is a helper function for AddAsiteCountsToTranscriptPosToCodonPosAllGenes 
    #' TidyAsiteCountsByPosition() and TranscriptPosToCodonPos() are used in this function.
    #' 
    #' The function joins the tables that are output by TidyAsiteCountsByPosition and 
    #' TranscriptPosToCodonPos functions to generate transcript_info_tibble   
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
    #' AddAsiteCountsToTranscriptPosToCodonPos(gene = "YAL003W", dataset = Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE, gff_df)
    #' 
    #' @export
    AddAsiteCountsToTranscriptPosToCodonPos <- function(gene, 
                                                        dataset, 
                                                        hd_file, 
                                                        min_read_length, 
                                                        colsum_out, 
                                                        gff_df){
      
      tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, 
                                                           dataset, 
                                                           hd_file, 
                                                           min_read_length, 
                                                           colsum_out)
      
      transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene, 
                                                                    gff_df)
      
      transcript_info_tibble <- dplyr::left_join(transcript_pos_to_codon_pos_output, 
                                                 tidy_asite_count_output, 
                                                 by = "Pos")
      
      return(transcript_info_tibble)
      
    }
    #TEST: AddAsiteCountsToTranscriptPosToCodonPos(): creates a tidy format data frame (tibble) = TRUE
    #TEST: AddAsiteCountsToTranscriptPosToCodonPos(): the tibble contains 5 columns = TRUE
    #TEST: AddAsiteCountsToTranscriptPosToCodonPos(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
    #TEST: AddAsiteCountsToTranscriptPosToCodonPos(): the column names are %in% c("Gene", "Pos", "Pos_Codon", "Frame", "Count") 
    # gives:
    # > str(transcript_info_tibble)
    # Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 5 variables:
    #   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    #   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
    #   $ Pos_Codon : int  NA NA NA NA NA NA NA NA NA NA ...
    #   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
    #   $ Count     : num  0 0 0 0 0 0 0 0 0 0 ...
    
    transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      colsum_out, 
                                                                      gff_df)
    
      # End result is a table which contains the asite counts, reading frame, codons and 
      # codon positions for the gene of interest
    
    
    
    #' AddCodonNamesToTranscriptInfoTibble(): Provides codon names from the .tsv file to the 
    #' tidy format data frame (tibble) generated from AddAsiteCountsToTranscriptPosToCodonPos   
    #' 
    #' AddAsiteCountsToTranscriptPosToCodonPos() is used in this function
    #' 
    #' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
    #' @param gene from gene_names   
    #' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
    #' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
    #' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
    #' @param gff_df from which to extract the UTRs and CDS widths.
    #' 
    #' @example 
    #' 
    #' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
    #' 
    #' AddCodonNamesToTranscriptInfoTibble(yeast_codon_pos_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df)
    #' 
    #' @export   
    AddCodonNamesToTranscriptInfoTibble <- function(yeast_codon_pos_i200, 
                                                    gene, 
                                                    dataset, 
                                                    hd_file, 
                                                    min_read_length, 
                                                    colsum_out, 
                                                    gff_df){
      
      codon_table <- dplyr::filter(yeast_codon_pos_i200, Gene == gene) %>%
        dplyr::select("PosCodon", "Codon")
      # this provides the positions of the codons and the codon names which can 
      # joined to the tibble generated by AddAsiteCountsToTranscriptPosToCodonPos 
      
      transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, 
                                                                        dataset, 
                                                                        hd_file, 
                                                                        min_read_length, 
                                                                        colsum_out, 
                                                                        gff_df)
    
      transcript_info_tibble <- left_join(transcript_info_tibble, 
                                          codon_table, 
                                          by = c("Pos_Codon" = "PosCodon"), 
                                          keep = FALSE)
      
      return(transcript_info_tibble)
      
    }
    #TEST: AddCodonNamesToTranscriptInfoTibble(): creates a tidy format data frame (tibble) = TRUE
    #TEST: AddCodonNamesToTranscriptInfoTibble(): the tibble contains 6 columns = TRUE
    #TEST: AddCodonNamesToTranscriptInfoTibble(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
    #TEST: AddCodonNamesToTranscriptInfoTibble(): the column names are %in% c("Gene", "Pos", "Pos_Codon", "Frame", "Count", "Codon") 
    #TEST: AddCodonNamesToTranscriptInfoTibble(): the first and last 250 rows of the columns "Pos_Codon" and "Codon" for each gene contain NA. 
    # gives:
    # > str(transcript_info_tibble)
    # Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 6 variables:
    #   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    #   $ Pos      : int  1 2 3 4 5 6 7 8 9 10 ...
    #   $ Pos_Codon: num  NA NA NA NA NA NA NA NA NA NA ...
    #   $ Frame    : num  2 0 1 2 0 1 2 0 1 2 ...
    #   $ Count    : num  0 0 0 0 0 0 0 0 0 0 ...
    #   $ Codon    : chr  NA NA NA NA ...
    
    transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(yeast_codon_pos_i200, 
                                                                  gene, 
                                                                  dataset, 
                                                                  hd_file, 
                                                                  min_read_length, 
                                                                  colsum_out, 
                                                                  gff_df)

    
    
}

# CreateTranscriptInfoTibbleAllGenes(reads_pos_length, gene, dataset, hd_file, gff_df)
# this takes all of the functions defined above and applies them to all of the genes in the sample. 

transcript_info_tibble <- suppressMessages(purrr::map_dfr(.x = gene_names, 
                                                          .f = CreateTranscriptInfoTibbleAllGenes, 
                                                          dataset, 
                                                          hd_file, 
                                                          gff_df, 
                                                          colsum_out,
                                                          startpos, 
                                                          startlen))

    # > str(transcript_info_tibble)
    # tibble [10,747 x 6] (S3: tbl_df/tbl/data.frame)
    # $ Gene     : chr [1:10747] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
    # $ Pos      : int [1:10747] 1 2 3 4 5 6 7 8 9 10 ...
    # $ Pos_Codon: num [1:10747] NA NA NA NA NA NA NA NA NA NA ...
    # $ Frame    : num [1:10747] 2 0 1 2 0 1 2 0 1 2 ...
    # $ Count    : num [1:10747] 0 0 0 0 0 0 0 0 0 0 ...
    # $ Codon    : chr [1:10747] NA NA NA NA ...

# For each gene, the first and last 250 nt may have NA in the codon column and Pos_Codon
# This is due to the 250 nt utr buffer, present in many riboviz example datasets 


# TEST:: Transcript_info_tibble should be a tibble
# class(transcript_info_tibble)
# [1] "tbl_df"     "tbl"        "data.frame"

# TEST:: Should contain the same number of genes as listed in the h5 file
# length(unique(transcript_info_tibble$Gene)) == length(gene_names)
# [1] TRUE

# TEST:: The genes included should be the same as those listed in the h5 file
# all(unique(transcript_info_tibble$Gene) %in% gene_names)
# [1] TRUE

# TEST:: should contain 6 columns 
# ncol(transcript_info_tibble)
# [1] 6

# TEST:: Should have 3 different frames listed 
# length(unique(transcript_info_tibble$Frame)) 
# [1] 3



#' FilterForFrame(): Filter for reading frame of interest
#' 
#' The default reading frame is 0
#' FIXME: Don't want to have the output from a function as the input for another function
#' 
#' @param transcript_info_tibble from AddCodonNamesToTranscriptInfoTibble 
#' @param filtering_frame which sets the reading frame that you filter for 
#' 
#' @example 
#' 
#' FilterForFrame(transcript_info_tibble, filtering_frame = 0)
#' 
#' @export
FilterForFrame <- function(transcript_info_tibble, filtering_frame){
  
  transcript_info_tibble <- dplyr::filter(transcript_info_tibble, 
                                          Frame == filtering_frame)
  
  return(transcript_info_tibble)
  
}
#TEST: FilterForFrame(): creates a tidy format data frame (tibble) = TRUE
#TEST: FilterForFrame(): the tibble contains 6 columns = TRUE
#TEST: FilterForFrame(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df divided by 3, for YAL003W = 1121/3 = ~374
#TEST: FilterForFrame(): The "Frame" column only contains the parameter you filtered for, e.g. 0.
#TEST: FilterForFrame(): the column names are %in% c("Gene", "Pos", "Pos_Codon", "Frame", "Count", "Codon") 
#TEST: FilterForFrame(): the first and last ~83 rows of the columns "Pos_Codon" and "Codon" for each gene contain NA. 
# gives: 
# > str(transcript_info_tibble) for YAL003W
# Classes 'tbl_df', 'tbl' and 'data.frame':   374 observations of 6 variables
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos       : int  2 5 8 11 14 17 20 23 26 29 ...
#   $ Pos_Codon : num  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  0 0 0 0 0 0 0 0 0 0 ...
#   $ Count     : num  0 0 0 0 0 0 0 0 0 0 ...
#   $ Codon     : chr  NA NA NA NA ...
#
# > str(transcript_info_tibble) for sim-YAL5
# Classes 'tbl_df', 'tbl' and 'data.frame':   3,584 observations of 6 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos      : int  2 5 8 11 14 17 20 23 26 29 ...
#   $ Pos_Codon: num  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame    : num  0 0 0 0 0 0 0 0 0 0 ...
#   $ Count    : num  0 0 0 0 0 0 0 0 0 0 ...
#   $ Codon    : chr  NA NA NA NA ...

transcript_info_tibble <- FilterForFrame(transcript_info_tibble, 
                                         filtering_frame)


print('Filtering for codon of interest')



# Function to filter for codons of interest, generates a list of the positions of the codon of interest

#' FilterForCodonOfInterestPositions(): Filters for the codon of interest 
#' 
#' FIXME: takes transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes() as its input
#' 
#' @param gene in gene_names to pull information from the gff_df file 
#' @param transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes()
#' @param codon_of_interest character, each incidence of the codon pair will be extracted from transcript_info_tibble
#' 
#' @example 
#' 
#' FilterForCodonOfInterestPositions(transcript_info_tibble, codon_of_interest = "CGA", gene_name)
#' 
#' @export
FilterForCodonOfInterestPositions <- function(transcript_info_tibble, codon_of_interest, gene_name){

  # transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(yeast_codon_pos_i200,
  #                                                               gene,
  #                                                               dataset,
  #                                                               hd_file,
  #                                                               min_read_length,
  #                                                               colsum_out,
  #                                                               gff_df)

  # generating the transcript_info_tibble inside of the function seems to result
  # in the list that is generated being triplicated. Skipping this step seems to
  # resolve that issue

  interesting_codon_table <- dplyr::filter(transcript_info_tibble, Codon == codon_of_interest, Gene == gene_name)
  interesting_codon_positions <- interesting_codon_table$Pos_Codon

  return(interesting_codon_positions)
}
#TEST: FilterForCodonPairOfInterestPositions(): returns a list of numeric values = TRUE
#gives:
# > str(interesting_codon_positions) for CCT in YAL003W
# num [1:7] 3 7 49 57 86 121 180

# interesting_codon_table <- dplyr::filter(transcript_info_tibble, Codon == codon_of_interest)
# interesting_codon_positions <- FilterForCodonOfInterestPositions(transcript_info_tibble, codon_of_interest, gene_name = "YAL003W")



### Expand frame around codon of interest ###



#' ExpandCodonPairRegion(): Expands the window around the codon position of interest
#' 
#' This function uses codon position values, not nucleotide position values
#' Helperfunction for ExpandCodonPairRegionForList()
#' 
#' This function is not looped and will therefore only output one expanded region when used on its own
#' 
#' FIXME: takes transcript_info_tibble as its input, do not want outputs from other functions as input
#' 
#' The expanded windows surrounding the positions of interest allow you to compare the rate of 
#' translation for the feature of interest to its surrounding translational environment 
#' 
#' @param gene in gene_names to pull information from the gff_df file 
#' @param transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes()
#' @param interesting_codon_positions list of numeric values describing position of the feature of interest
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' @param remove_overhang default = TRUE, removes features of interest who's positions do not allow for complete windows to be generated
#' 
#' @example 
#' 
#' ExpandCodonRegion(.x = interesting_codon_positions, 
#'                       transcript_info_tibble, 
#'                       gene = "YAL003W", 
#'                       dataset = "Mok-simYAL5", 
#'                       hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), 
#'                       expand_width = 5L, 
#'                       remove_overhang = TRUE)
#'                       
#' @export     
ExpandCodonRegion <- function(.x = interesting_codon_positions, transcript_info_tibble_gene, gene, dataset, hd_file, expand_width, remove_overhang = TRUE) {

  gene_poscodon_count <- na.omit(tibble(
    Gene = transcript_info_tibble_gene$Gene,
    Pos_Codon = transcript_info_tibble_gene$Pos_Codon,
    Rel_Count = transcript_info_tibble_gene$Count
  ))

  # added na.omit() as the slice function goes by the row number instead of Pos_codon
  # na.omit() removes the UTRs so only CDS remains

  gene_length <- GetGeneLength(gene, dataset, hd_file) # giving me 618??

  if (.x <= expand_width  |.x + expand_width > gene_length/3) {
    return()
    }else{
    output_codon_info <- tibble(
      dplyr::slice(gene_poscodon_count, (.x - expand_width):(.x + expand_width), each = FALSE),
      Rel_Pos =  seq(- expand_width, + expand_width)
    )
    
    if(dim(output_codon_info)[1] == (2*expand_width + 1)){
      
      return(output_codon_info)
    }else{
      return()
    }
    
  }
  
  # The if statement ensures that codon positions that are less/more than the 
  # expand_width value are discarded 
}
#TEST: ExpandCodonRegion(): creates a tidy format data frame (tibble) = TRUE
#TEST: ExpandCodonRegion(): the tibble contains 4 columns = TRUE
#TEST: ExpandCodonRegion(): the column names are %in% c("Gene", "Pos_Codon", "Rel_Count", "Rel_Pos")
#TEST: ExpandCodonRegion(): number of observations in the output tibble = "expand_width" * 2 + 1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandCodonRegion(): the position from "interesting_codon_positions" has "Rel_Pos" value 0 = TRUE
#TEST: ExpandCodonRegion(): the column "Rel_Pos" goes from -"expand_width to +"expand_width"
#gives:
# > str(output_codon_info) for "CGA" in YAL035W
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 4 variables
#   $ Gene     : chr  "YAL005C" "YAL005C" "YAL005C" "YAL005C" "YAL005C"
#   $ Pos_Codon: num  18 19 20 21 22 23 24 25 26 27 28
#   $ Rel_Count: num  32 33 26 21 41 44 152 57 1 2 26
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 5


    # output_codon_info <- ExpandCodonRegion(.x = interesting_codon_positions,
    #                                        transcript_info_tibble,
    #                                        gene = "YAL035W",
    #                                        dataset,
    #                                        hd_file,
    #                                        expand_width = 5L,
    #                                        remove_overhang = TRUE)



#' ExpandCodonPairRegionForList(): Applies the ExpandCodonRegion function to all of the codons of interest. 
#' 
#' This function generates expanded tibbles for each position contained within the codon_of_interest list.
#' Each occurence in codon_of_interest list generates one tibble 
#' 
#' FIXME: Takes codon_of_interest and transcript_info_tibble as inputs 
#' 
#' @param transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes()
#' @param gene in gene_names to pull information from the gff_df file 
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param startpos position of the start codon, default = 1
#' @param startlen, smallest length of reads', default = 10
#' @param codon_of_interest character, each incidence of the codon pair will be extracted from transcript_info_tibble
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' 
#' @example 
#' 
#' ExpandCodonPairRegionForList(transcript_info_tibble, 
#'                              gene = unique(transcript_info_tibble$Gene), 
#'                              dataset = "Mok-simYAL5", 
#'                              hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5") 
#'                              startpos = 1, 
#'                              startlen = 10, 
#'                              codon_of_interest = "TCC", 
#'                              gff_df, 
#'                              expand_width = 5L)
#' 
#' @export
ExpandCodonRegionForList <- function(transcript_info_tibble, gene, dataset, hd_file, startpos, startlen, codon_of_interest, gff_df, expand_width){
  
  expand_codon_region <- list()
  
  # A loop is used here so each unique gene can be processed separately. there were issues with using map and lapply 
  # on the filter function of FilterForCodonOfInterestPositions, and filtering using only one gene would produce NULLs
  
  for(gene_name in unique(gene)){
    
    interesting_codon_positions <- FilterForCodonOfInterestPositions(transcript_info_tibble, codon_of_interest, gene_name)
    
    transcript_info_tibble_gene <- dplyr::filter(transcript_info_tibble, Gene == gene_name)
    
    tmp_expand_codon_region <- purrr::map(
      .x = interesting_codon_positions, # vector of codon positions to iterate over for single codons
      .f = ExpandCodonRegion,   # function to use at each codon position of interest for single codons
      transcript_info_tibble_gene,
      gene = gene_name,
      dataset,
      hd_file,
      gff_df,
      expand_width = expand_width
    ) 
    
    # remove any tibbles from the list that are null, so are within one expand_width of the UTRs
    tmp_expand_codon_region <- tmp_expand_codon_region[!sapply(tmp_expand_codon_region,is.null)]
    
    
    expand_codon_region <- append(expand_codon_region, tmp_expand_codon_region)
  }
  return(expand_codon_region)
}
#TEST: ExpandCodonRegionForList(): output is a list of tidy format data frames (tibbles) = TRUE. type(expand_codon_region) = "list"
#TEST: ExpandCodonRegionForList(): number of tibbles in list matches the number of occurrences in "codon_of_interest" list = TRUE
#TEST: ExpandCodonRegionForList(): each tibble contains 4 columns = TRUE
#TEST: ExpandCodonRegionForList(): the column names are %in% c("Gene", "Pos_Codon", "Rel_Count", "Rel_Pos")
#TEST: ExpandCodonRegionForList(): number of observations in each output tibble = "expand_width" * 2 + 1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandCodonRegionForList(): the position from "interesting_codon_positions" has "Rel_Pos" value 0 = TRUE
#TEST: ExpandCodonRegionForList(): the column "Rel_Pos" goes from -"expand_width to +"expand_width"
#gives:
# > str(expand_codon_region)
# List of 68 (shows 3 of 68 occurrences)
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 4 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos_Codon: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Rel_Count: num  429 488 102 994 146 173 762 13 176 98 ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 4 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos_Codon: num  44 45 46 47 48 49 50 51 52 53 ...
#   $ Rel_Count: num  128 279 417 389 237 ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 4 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos_Codon: num  52 53 54 55 56 57 58 59 60 61 ...
#   $ Rel_Count: num  42 53 648 293 121 92 519 79 765 196 ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...

expand_codon_region <- ExpandCodonRegionForList(transcript_info_tibble, 
                                                gene = unique(transcript_info_tibble$Gene), 
                                                dataset, 
                                                hd_file, 
                                                startpos, 
                                                startlen, 
                                                codon_of_interest, 
                                                gff_df, 
                                                expand_width)


# the dimensions of each item in the list shoud be [(2*expand_width+1)X4] as there are 4 rows; Gene, Pos_Codon, Rel_Count, Rel_Pos
# and one column for each position being included ie (-5 to 5) relative to the stop codon
# TEST:: Where an expand_width of =5 is given,the following would be observed:
# > dim(expand_codon_region[[1]])
# [1] 11  4



### Normalization ###

# Normalization carried out within each expanded frame so that they are comparable 
# Normalizes the expand_codon_region list generating a RelCount column with the normalization values

print('Normalising data')

#' ExpandedRegionNormalization(): carries out normalization within each expanded frame so that they are comparable
#' 
#' Normalizes the ExpandCodonRegion list generating a RelCount column with the normalization values.
#' As this function is not looped it will only generate one normalized tibble for one occurence 
#' 
#' @param .x which is the list of tidy format data frames (tibbles) generated by the function ExpandCodonRegion
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' 
#' @example 
#' 
#' ExpandedCodonRegionNormalization(expand_codon_region, expand_width = 5L)
#' 
ExpandedRegionNormalization <- function(.x, expand_width){
  
  # dplyr::mutate(.x, RelCount = PerCodonCounts / sum(PerCodonCounts) * (2 * expand_width + 1))
  dplyr::mutate(.x, RelCount = Rel_Count / sum(Rel_Count) * (2 * expand_width + 1))
  
}
#TEST: ExpandedRegionNormalization(): creates a tidy format data frame (tibble) = TRUE
#TEST: ExpandedRegionNormalization(): the tibble contains 5 columns = TRUE
#TEST: ExpandedRegionNormalization(): the column names are %in% c("Gene", "Pos_Codon", "Rel_Count", "Rel_Pos", "RelCount")
#TEST: ExpandedRegionNormalization(): number of observations in the output tibble = "expand_width" * 2 + 1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandedRegionNormalization(): the position from "interesting_codon_positions" has "Rel_Pos" value 0 = TRUE
#TEST: ExpandedRegionNormalization(): the column "Rel_Pos" goes from -"expand_width to +"expand_width" 
#gives:
# > str(normalized_expanded_codon_region)
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 5 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos_Codon: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Rel_Count: num  429 488 102 994 146 173 762 13 176 98 ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
#   $ RelCount : num  1.347 1.532 0.32 3.12 0.458 ...

# normalized_expanded_codon_region <- ExpandedRegionNormalization(expand_codon_region[[1]], expand_width)



# Normalization carried out for all the tibbles within ExpandList 
normalized_expand_list <- purrr::map(
  .x = expand_codon_region,
  .f = ExpandedCodonRegionNormalization,
  expand_width
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


  # > str(normalized_expand_list)
  # List of 6
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 13 14 15 16 17 18 19 20 21 22 ...
  # ..$ Rel_Count: num [1:11] 490 561 51 567 67 65 158 71 323 288 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 1.367 1.565 0.142 1.582 0.187 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 26 27 28 29 30 31 32 33 34 35 ...
  # ..$ Rel_Count: num [1:11] 224 429 201 94 0 54 118 485 242 50 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 1.067 2.043 0.957 0.448 0 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 38 39 40 41 42 43 44 45 46 47 ...
  # ..$ Rel_Count: num [1:11] 64 251 128 191 232 89 128 279 417 389 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 0.293 1.148 0.585 0.874 1.061 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 59 60 61 62 63 64 65 66 67 68 ...
  # ..$ Rel_Count: num [1:11] 79 765 196 184 626 84 62 395 104 72 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 0.312 3.026 0.775 0.728 2.476 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 64 65 66 67 68 69 70 71 72 73 ...
  # ..$ Rel_Count: num [1:11] 84 62 395 104 72 214 99 163 128 69 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 0.435 0.321 2.044 0.538 0.373 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 192 193 194 195 196 197 198 199 200 201 ...
  # ..$ Rel_Count: num [1:11] 364 960 162 52 115 58 89 235 7 86 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 1.686 4.446 0.75 0.241 0.533 ...



### Overlaying the normalized expanded tibbles 

print('Calculating average')

# Function to overlay graphs into a single graph. Need to generate a single tibble 
# from NormalizedExpandList. Need to join by Rel_Pos, in RelCount need the mean for 
# each Rel_Pos (sum row(x) / number of row(x))

OverlayedTable <- function(normalized_expand_list, expand_width){
  number_of_objects <- length(normalized_expand_list)
  # The number of objects inside normalized_expand_list
  
  result <- lapply(normalized_expand_list, "[", c("Rel_Pos", "RelCount"))
  # Reduces normalized_expand_list to the columns Rel_Pos and RelCount
  
  joined_result <- result %>% purrr::reduce(full_join, by = c("Rel_Pos"), sum("RelCount"))
  
  joined_rows = joined_result %>% 
    mutate(SumRows = rowSums(select(joined_result, -"Rel_Pos")) / number_of_objects)
  
  overlayed_tibbles <- tibble::tibble(
    Rel_Pos = seq(- expand_width, expand_width),
    RelCount = joined_rows$SumRows
  )
}

# something is wrong within the overlaying function, need to figure out why this isnt working anymore  

Over <- OverlayedTable(normalized_expand_list, expand_width) 

  # output from OverlayedTable function for single codons:
  # # > Over
  # # # A tibble: 11 x 2
  # # Rel_Pos RelCount
  # # <int>    <dbl>
  # #   1      -5    0.860
  # # 2      -4    2.09
  # # 3      -3    0.876
  # # 4      -2    0.735
  # # 5      -1    0.772
  # # 6       0    0.426
  # # 7       1    0.460
  # # 8       2    1.21
  # # 9       3    0.845
  # # 10       4    0.644
  # # 11       5    2.08

# At each position, the sum of RelCount should be equal to (2*expand_width+1)
# ie if the expand width was 5:
# sum(Over$RelCount)
# [1] 11

# TEST:: max Rel_Pos should equal expand_width and min Rel_Pos should equal -expand_width
# max(Over$Rel_Pos)
# [1] 5
# min(Over$Rel_Pos)
# [1] -5

# TEST:: Over should be a tibble
# class(Over)
# [1] "tbl_df"     "tbl"        "data.frame"

print('Creating plot')

overlayed_plot <- ggplot(Over, mapping = aes(x = Rel_Pos, y = RelCount)) + 
  geom_line() +
  theme_bw()+
  theme(text=element_text(size=14),
        axis.title=element_text(size=14, face='bold'),
        title = element_text(size = 14, face='bold'))+
  labs(title = paste0('Relative read counts around codon ', codon_of_interest),
       x = 'Position relative to codon of interest',
       y = 'Relative number of reads', size = 2)

save_plot_pdf <- function(overlayed_plot, output_dir){
  overlayed_plot %>%
    ggsave(
      filename = file.path(output_dir,"inhibitory_codon_plot.pdf"),
      width = 6, height = 5
    )
}

save_plot_pdf(overlayed_plot, output_dir)

print('Done')
