
# YAL5_h5 is at location $HOME/riboviz/riboviz/Mok-simYAL5/output/A/A.h5

# Given an h5 file, GFF file and .tsv file, this script creates a metafeature plot 
# for the codon/codon pair of interest. 

print('Starting process')

# source packages and functions from rscripts 
suppressMessages(source(here::here("rscripts", "read_count_functions.R")))
suppressMessages(source(here::here("rscripts", "stats_figs_block_functions.R")))


suppressMessages(library(tidyverse))
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
parser$add_argument('--colsum_out', help='logical', default = TRUE)

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
filtering_frame <- args$frame
min_read_length <- args$minreadlen
colsum_out <- args$colsum_out


hd_file <- here::here("Mok-simYAL5", "output", "A", "A.h5")
dataset <- "Mok-simYAL5"
codon_pair_of_interest <- 'TCC AAG'
expand_width = 5L
startpos <-1
startlen <- 10
filtering_frame <- 0
min_read_length <- 10
yeast_codon_table <- here::here("data", "yeast_codon_table.tsv")
gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3")
colsum_out <- TRUE

#

gff_df <- readGFFAsDf(gff)
gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name


# Import the .tsv file: 
yeast_codon_pos_i200 <- suppressMessages(readr::read_tsv(file = yeast_codon_table))

    # > str(yeast_codon_pos_i200)
    # spec_tbl_df [2,826,757 x 3] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
    # $ Gene    : chr [1:2826757] "YAL068C" "YAL068C" "YAL068C" "YAL068C" ...
    # $ PosCodon: num [1:2826757] 1 2 3 4 5 6 7 8 9 10 ...
    # $ Codon   : chr [1:2826757] "ATG" "GTC" "AAA" "TTA" ...


# The yeast_codon_pos_i200 file is configured to show the lead and lag codon pair positions:
gene_poscodon_codon_i200 <- tibble::tibble(
  Gene = yeast_codon_pos_i200$Gene,
  CodonPos_1 = yeast_codon_pos_i200$PosCodon, 
  CodonPos_2 = dplyr::lead(yeast_codon_pos_i200$PosCodon),
  CodonPair = paste(yeast_codon_pos_i200$Codon, dplyr::lead(yeast_codon_pos_i200$Codon))
)
  
    # > str(gene_poscodon_codon_i200)
    # tibble [2,826,757 x 4] (S3: tbl_df/tbl/data.frame)
    # $ Gene      : chr [1:2826757] "YAL068C" "YAL068C" "YAL068C" "YAL068C" ...
    # $ CodonPos_1: num [1:2826757] 1 2 3 4 5 6 7 8 9 10 ...
    # $ CodonPos_2: num [1:2826757] 2 3 4 5 6 7 8 9 10 11 ...
    # $ CodonPair : chr [1:2826757] "ATG GTC" "GTC AAA" "AAA TTA" "TTA ACT" ...


#####


print('Creating information tibble')


# The function CreateTranscriptInfoTibbleAllGenes applies all of the functions up until and including filtering for frame of interest to all of the genes 
# it returns one tibble, containing all of the genes and information
# to run for a single gene or to run each function separately, filter for a gene of interest and run

# Filter down the gene_poscodon_codon_i200 file to gene of interest, e.g. YAL003W
YAL003W_pos <- dplyr::filter(gene_poscodon_codon_i200, Gene=="YAL003W")
gene <- 'YAL003W'


CreateTranscriptInfoTibbleAllGenes <- function(gene, dataset, hd_file, gff_df, startpos, startlen, colsum_out){

          
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
TidyAsiteCountsByPosition <- function(gene, dataset, hd_file, min_read_length, colsum_out){
  
  asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))
  #This step extracts the asite reads from the reads_pos_length file  
  
    # > str(asite_displacement_length)
    # spec_tbl_df [3 x 2] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
    # $ read_length       : num [1:3] 28 29 30
    # $ asite_displacement: num [1:3] 15 15 15
            
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  
    # > str(reads_pos_length) (for YAL003W)
    # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
              
  asite_counts_by_position <- CalcAsiteFixed(reads_pos_length,
                                             min_read_length = min_read_length,
                                             asite_displacement_length = asite_displacement_length,
                                             colsum_out = colsum_out
                                             )
    # > str(asite_counts_by_position)(for YAL003W)
    # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
            
  tidy_asite_counts <- tibble(Pos = 1:length(asite_counts_by_position), Count = asite_counts_by_position)
            
  return(tidy_asite_counts)
            
}
#TEST: TidyAsiteCountsByPosition(): returns a a tidy format data frame (tibble). 
#TEST: TidyAsiteCountsByPosition(): the tibble has 2 columns.
#TEST: TidyAsiteCountsByPosition(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df for the gene.
#TEST: TidyAsiteCountsByPosition(): the column names are %in% c("Pos", "Count")
# gives:          
# > str(tidy_asite_count_output)
# Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 2 variables:
# $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
# $ Count: num  0 0 0 0 0 0 0 0 0 0 ...
          
tidy_asite_counts <- TidyAsiteCountsByPosition(gene, dataset, hd_file, min_read_length, colsum_out)
          
      # The end result here is that the A-site assigned counts are aligned to the gene of interest
      # (in nucleotides, including UTRs and UTRs)
          

#####


 
### Functions for assembling tibbles consisting of counts, positions and codons ###
  


#' TranscriptPosToCodonPos(): Assigns the nucleotide positions of UTRs and CDS to codon positions
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
      Pos_Codon1 = c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)),
      Pos_Codon2 = dplyr::lead(c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)), n = 3),
      Frame = seq(from = 2L, to = (transcript_length$end + 1L)) %% 3 # works as UTRS are 250, might not work with UTRS of others values
      # add the count column, likely something we want here 
    )
    
    # For codon pairs name of Pos_Codon became Pos_Codon1 and Pos_Codon2 was added 
    # where the line for Pos_Codon1 was copied and dplyr::lead((), n = 3) was added 
    
    return(transcript_gene_pos_poscodon_frame)
}
#TEST: TranscriptPosToCodonPos(): creates a tibble = TRUE
#TEST: TranscriptPosToCodonPos(): the tibble contains 5 columns = TRUE
#TEST: TranscriptPosToCodonPos(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
#TEST: TranscriptPosToCodonPos(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame") 
#TEST: TranscriptPosToCodonPos(): The first and last 250 rows of the columns "Pos_codon1" and "Pos_Codon2" for each gene contain NA. 
# gives:
# > str(transcript_pos_to_codon_pos_output)
# Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 5 variables:
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...


transcript_gene_pos_poscodon_frame <- TranscriptPosToCodonPos(gene, gff_df)

  # The end result is a table with the columns Gene, Pos, Pos_Codon and Frame (reading frame).
  # For the Pos_Codon column the UTR positions have NA, while the CDS has the codon positions 
             

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
#' AddAsiteCountsToTranscriptPosToCodonPos(gene = "YAL003W", dataset = Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE, gff_df)
#' 
#' @export
AddAsiteCountsToTranscriptPosToCodonPos <- function(gene, dataset, hd_file, min_read_length, colsum_out, gff_df){
    
    tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, dataset, hd_file, min_read_length, colsum_out = colsum_out)
    
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

    # End result is a table which contains the asite counts, reading frame, codons and 
    # codon positions for the gene of interest
  
transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)



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
#' AddCodonNamesToTranscriptInfoTibble(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df)
#' 
#' @export                                                                     
AddCodonNamesToTranscriptInfoTibble <- function(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df){
    
    codon_table <- dplyr::filter(gene_poscodon_codon_i200, Gene == gene) %>% 
      dplyr::select("CodonPos_1", "CodonPos_2", "CodonPair")
    # this provides the positions of the codons and the codon names 
    
    transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, dataset = dataset, hd_file = hd_file, min_read_length = min_read_length, colsum_out = colsum_out, gff_df)
    # this provides the tibble with the columns Gene, Pos, Pos_Codon and Frame (reading frame)
    
    transcript_info_tibble <- left_join(transcript_info_tibble, codon_table, by = c("Pos_Codon1" = "CodonPos_1", "Pos_Codon2" = "CodonPos_2"), keep = FALSE)

    return(transcript_info_tibble)
}
#TEST: AddCodonNamesToTranscriptInfoTibble(): creates a tidy format data frame (tibble) = TRUE
#TEST: AddCodonNamesToTranscriptInfoTibble(): the tibble contains 7 columns = TRUE
#TEST: AddCodonNamesToTranscriptInfoTibble(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
#TEST: AddCodonNamesToTranscriptInfoTibble(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame", "Count", "CodonPair") 
#TEST: AddCodonNamesToTranscriptInfoTibble(): the first and last 250 rows of the columns "Pos_codon1", "Pos_Codon2" and "CodonPair" for each gene contain NA. 
# gives:
# > str(transcript_info_tibble)
# Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 7 variables:
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Pos_Codon1: num  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: num  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
#   $ Count     : num  0 0 0 0 0 0 0 0 0 0 ...
#   $ CodonPair : chr  NA NA NA NA ...

transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df)
 

 
}

# single_gene_info_tibble <- CreateTranscriptInfoTibbleAllGenes(gene, dataset, hd_file, gff_df, filtering_frame, startpos, startlen)

# CreateTranscriptInfoTibbleAllGenes (reads_pos_length, gene, dataset, hd_file, gff_df)
transcript_info_tibble <- suppressMessages(purrr::map_dfr(.x = gene_names, 
                                                          .f = CreateTranscriptInfoTibbleAllGenes, 
                                                          dataset, 
                                                          hd_file, 
                                                          gff_df, 
                                                          colsum_out,
                                                          startpos, 
                                                          startlen))

# this takes all of the functions defined above and applies them to all of the genes in the sample. 

# For each gene, the first and last 250 nt may have NA in the codon column and Pos_Codon
# This is due to the 250 nt utr buffer, present in many riboviz example datasets 

# > str(transcript_info_tibble)
# tibble [3,584 x 7] (S3: tbl_df/tbl/data.frame)
# $ Gene      : chr [1:3584] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
# $ Pos       : int [1:3584] 2 5 8 11 14 17 20 23 26 29 ...
# $ Pos_Codon1: num [1:3584] NA NA NA NA NA NA NA NA NA NA ...
# $ Pos_Codon2: num [1:3584] NA NA NA NA NA NA NA NA NA NA ...
# $ Frame     : num [1:3584] 0 0 0 0 0 0 0 0 0 0 ...
# $ Count     : num [1:3584] 0 0 0 0 0 0 0 0 0 0 ...
# $ CodonPair : chr [1:3584] NA NA NA NA ...



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
#TEST: FilterForFrame(): the tibble contains 7 columns = TRUE
#TEST: FilterForFrame(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df divided by 3, for YAL003W = 1121/3 = ~374
#TEST: FilterForFrame(): The "Frame" column only contains the paramter you filtered for, e.g. 0.
#TEST: AddCodonNamesToTranscriptInfoTibble(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame", "Count", "CodonPair") 
#TEST: AddCodonNamesToTranscriptInfoTibble(): the first and last ~83 rows of the columns "Pos_codon1", "Pos_Codon2" and "CodonPair" for each gene contain NA. 
# gives: 
# > str(transcript_info_tibble)
# Classes 'tbl_df', 'tbl' and 'data.frame':   374 observations of 7 variables
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos       : int  2 5 8 11 14 17 20 23 26 29 ...
#   $ Pos_Codon1: num  NA NA NA NA NA NA NA NA NA NA ...
#   $ Pos_Codon2: num  NA NA NA NA NA NA NA NA NA NA ...
#   $ Frame     : num  0 0 0 0 0 0 0 0 0 0 ...
#   $ Count     : num  0 0 0 0 0 0 0 0 0 0 ...
#   $ CodonPair : chr  NA NA NA NA ...

transcript_info_tibble <- FilterForFrame(transcript_info_tibble, filtering_frame)

  
  ###



print('Filtering for codon of interest')


#' FilterForCodonPairOfInterestPositions(): Filters for the codon pairs of interest 
#' 
#' FIXME: takes transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes() as its input
#' 
#' @param gene in gene_names to pull information from the gff_df file 
#' @param transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes()
#' @param codon_pair_of_interest character, each incidence of the codon pair will be extracted from transcript_info_tibble
#' 
#' @example 
#' 
#' FilterForCodonPairOfInterestPositions(transcript_info_tibble, codon_pair_of_interest = "TCC AAG", gene_name)
#' 
#' @export
FilterForCodonPairOfInterestPositions <- function(transcript_info_tibble, codon_pair_of_interest, gene_name){
  
      interesting_codon_table <- dplyr::filter(transcript_info_tibble, CodonPair == codon_pair_of_interest & Gene == gene_name)
      # Changed Codon == codon_of_interest to CodonPair = codon_of_interest
  
      interesting_first_codon_positions <- interesting_codon_table$Pos_Codon1
  
      # return(interesting_codon_positions)
      return(interesting_first_codon_positions)
}
#TEST: FilterForCodonPairOfInterestPositions(): returns a list of numeric values = TRUE
#gives:
  # > str(interesting_first_codon_positions)
  # num [1:8] 7 57 383 508 535 321 90 412
  
    # # Which refers to the position of the first codon in the codon pair as the second codon
    # # in the pair is the first codon position + 1

# interesting_first_codon_positions <- FilterForCodonPairOfInterestPositions(transcript_info_tibble, codon_pair_of_interest, gene_name)
  


### Expand frame around codon of interest 


#' ExpandCodonPairRegion(): Expands the window around the codon_pair position of interest
#' 
#' Helperfunction for ExpandCodonPairRegionForList()
#' 
#' FIXME: takes transcript_info_tibble as its input, do not want outputs from other functions as input
#' 
#' The expanded windows surrounding the positions of interest allow you to compare the rate of 
#' translation for the feature of interest to its surrounding translational environment 
#' 
#' @param gene in gene_names to pull information from the gff_df file 
#' @param transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes()
#' @param interesting_first_codon_positions list of numeric values describing position of the feature of interest
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' @param remove_overhang default = TRUE, removes features of interest who's positions do not allow for complete windows to be generated
#' 
#' @example 
#' 
#' ExpandCodonPairRegion(.x = interesting_first_codon_positions, 
#'                       transcript_info_tibble, 
#'                       gene = "YAL003W", 
#'                       dataset = "Mok-simYAL5", 
#'                       hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), 
#'                       expand_width = 5L, 
#'                       remove_overhang = TRUE)
#'                       
#' @export                       
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
#TEST:
#gives: 


output_codonpair_info <- ExpandCodonPairRegion(.x = interesting_first_codon_positions, transcript_info_tibble, gene, dataset, hd_file, expand_width = 5L, remove_overhang = TRUE)

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

#' ExpandCodonPairRegionForList(): Applies the ExpandCodonRegion function to all of the codons of interest. 
#' 
#' This function generates expanded tibbles for each position contained within the codon_pair_of_interest list.
#' 
#' FIXME: Takes codon_pair_of_interest and transcript_info_tibble as inputs 
#' 
#' @param transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes()
#' @param gene in gene_names to pull information from the gff_df file 
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param startpos position of the start codon, default = 1
#' @param startlen, smallest length of reads', default = 10
#' @param codon_pair_of_interest character, each incidence of the codon pair will be extracted from transcript_info_tibble
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
#'                              codon_pair_of_interest = "TCC AAG", 
#'                              gff_df, 
#'                              expand_width = 5L)
#' 
#' @export
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
  
    interesting_codon_positions <- FilterForCodonPairOfInterestPositions(transcript_info_tibble, codon_pair_of_interest = codon_pair_of_interest, gene_name = gene_name)  
  
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
#TEST:
#gives:

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



### Normalization


print('Normalising data')

#' ExpandedCodonRegionNormalization(): carries out normalization within each expanded frame so that they are comparable
#' 
#' Normalizes the ExpandCodonRegion list generating a RelCount column with the normalization values 
#' 
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


overlayed_plot <- ggplot(Over, mapping = aes(x = Rel_Pos, y = RelCount)) + 
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        title = element_text(size = 12, face = 'bold'))+
  labs(title = paste0('Meta-feature plot of codon pair ', codon_pair_of_interest),
       x = 'Distance from codon pair (3nt)',
       y = 'Normalized ribosomal occupancy', size = 2)
save_plot_pdf <- function(overlayed_plot, output_dir){
  overlayed_plot %>%
    ggsave(
      filename = file.path(output_dir,"inhibitory_codon_plot.pdf"),
      width = 6, height = 5
    )
}

save_plot_pdf(overlayed_plot, output_dir)

print('Done')




