# Given an h5 file, GFF file and .tsv file, this script creates a metafeature plot 
# for the feature of interest

## TEST: run on TinySim Dataset 

print('Starting process')

# source packages and functions from rscripts 
source(here::here("rscripts", "read_count_functions.R"))
source(here::here("rscripts", "stats_figs_block_functions.R"))

suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))

# parser <- ArgumentParser()
# 
# parser$add_argument('-i', '--input', help='Path input to h5 file')
# parser$add_argument('-d', '--dataset', help='Name of the dataset being studied')
# parser$add_argument('-g', '--gff', help='Path t the GFF3 file of the organism being studied')
# parser$add_argument('-a', '--annotation', help='Path to codon table for organism')
# parser$add_argument('--feature', help='Feature of interest, e.g. codon')
# parser$add_argument('-o', '--output', help='Path to output directory')
# parser$add_argument('--expand_width', help='the desired range either side of the feature of interest', default = 5)
# parser$add_argument('--startpos', help='position of the start codon', default = 1)
# parser$add_argument('--startlen', help='smallest length of reads', default = 10) # please correct if wrong
# parser$add_argument('--frame', help='reading frame to be studied', default = 0)
# parser$add_argument('--minreadlen', help='minimum read length', default = 10)
# parser$add_argument('--colsum_out', help='logical', default = TRUE)
# parser$add_argument('--filter_for_frame', help='Filter to include only the reads from the first nucleotide of a codon', default = FALSE)
# parser$add_argument('--snapdisp', help='frame to filer to when using SnapToCodon', default = 0L)
# 
# args <- parser$parse_args()

hd_file <- here::here("Mok-tinysim", "output", "A", "A.h5")
dataset <- "Mok-tinysim"
feature_of_interest <- 'GCC'
expand_width = 1L
startpos <-1
startlen <- 10
filtering_frame <- 0
min_read_length <- 10
yeast_codon_table <- here::here("Mok-tinysim", "tinysim_codon_pos.tsv")
gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "tiny_2genes_20utrs.gff3")
colsum_out <- TRUE
output_dir <- '.'
filter_for_frame <- FALSE
snapdisp <- 0L


# hd_file <- args$input
# dataset <- args$dataset
# gff <- args$gff
# yeast_codon_table <- args$annotation
# feature_of_interest <- args$feature
# output_dir <- args$output
# expand_width <- args$expand_width
# startpos <- args$startpos
# startlen <- args$startlen
# filtering_frame <- args$frame
# min_read_length <- args$minreadlen
# colsum_out <- args$olsum_out
# filter_for_frame <- args$filter_for_frame
# snapdisp <- args$snapdisp


gff_df <- readGFFAsDf(gff)

## TEST: When running on tidy sim data set, the gff appears as follows
# There are 2 genes, MAT and MIKE, which have CDSs flanked by 20nt 

# > gff_df
# # A tibble: 6 x 10
# seqnames start   end width strand source   type  score phase Name 
# <fct>    <int> <int> <int> <fct>  <fct>    <fct> <dbl> <int> <chr>
#   1 MAT          1    20    20 +      ewallace UTR5     NA    NA MAT  
# 2 MAT         21    32    12 +      ewallace CDS      NA    NA MAT  
# 3 MAT         33    52    20 +      ewallace UTR3     NA    NA MAT  
# 4 MIKE         1    20    20 +      ewallace UTR5     NA    NA MIKE 
# 5 MIKE        21    35    15 +      ewallace CDS      NA    NA MIKE 
# 6 MIKE        36    55    20 +      ewallace UTR3     NA    NA MIKE 


gene_names <- unique(gff_df$Name)


# Import the .tsv file: 

yeast_codon_pos_i200 <- suppressMessages(readr::read_tsv(file = yeast_codon_table))


# > str(yeast_codon_pos_i200)
# spec_tbl_df [2,826,757 x 3] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
# $ Gene    : chr [1:2826757] "YAL068C" "YAL068C" "YAL068C" "YAL068C" ...
# $ PosCodon: num [1:2826757] 1 2 3 4 5 6 7 8 9 10 ...
# $ Codon   : chr [1:2826757] "ATG" "GTC" "AAA" "TTA" ...

print('Get codon positions and read counts')



##### GetAllCodonPosCounts #####

GetGeneCodonPosReads1dsnap <- function(gene, dataset, hd_file, left, right,
                                       min_read_length,
                                       asite_displacement_length = data.frame(
                                         read_length = c(28, 29, 30),
                                         asite_displacement = c(15, 15, 15)
                                       ),
                                       snapdisp) {
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
  reads_asitepos <- CalcAsiteFixed(
    reads_pos_length,
    min_read_length,
    asite_displacement_length
  )
  SnapToCodon(reads_asitepos, left, right, snapdisp)
}

## TEST:: CalAsiteFixed:
#
# mat_datamatrix <- GetGeneDatamatrix('MAT', dataset, hd_file)
# mike_datamatrix <- GetGeneDatamatrix('MIKE', dataset, hd_file) 
#
# asite_mat <- CalcAsiteFixed(mat_datamatrix, min_read_length = 10, asite_displacement_length)
# mat_tibble <- tibble(counts =asite_mat,
#                      pos = 1:length(asite_mat))
#
# mat_tibble <- tibble(pos = 1:length(asite_mat)
#                      counts =asite_mat)
#
# 
# asite_mike <- CalcAsiteFixed(mike_datamatrix, min_read_length = 10, asite_displacement_length)
# mike_tibble <- tibble(counts = asite_mike,
#                       pos = 1:length(asite_mike))
# 
# mike_tibble <- tibble(pos = 1:length(asite_mike)
#                       counts =asite_mike)
#
# Expected contents of mat_tibble 
#
# pos   counts
# <int>  <dbl>
# 19    1
# 25    1
# 26    1
# 27    2
# 
# Expected contents of mike_tibble
#   
# pos   counts
# <int>  <dbl>
# 20    1
# 25    1


FilterForFrameFunction <- function(gene, dataset, hd_file, asite_displacement_length ,reads_asitepos, left, right){
  
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
  
  reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length, asite_displacement_length)
  
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
  # where gene = YAL003W
  
  left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start))
  # 251
  
  right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end))
  # 871
  
  cds <- reads_asitepos[left:right]
  # num [1:621] 811 460 2978 429 251 ...
  
  cds_length <- length(cds)/3
  # 207
  
  cds_frames <- tibble(Count = cds,
                       Frame = rep(c(0, 1, 2), times = cds_length)
  )
  # > str(cds_frames)
  # tibble [621 x 2] (S3: tbl_df/tbl/data.frame)
  # $ Count: num [1:621] 811 460 2978 429 251 ...
  # $ Frame: num [1:621] 0 1 2 0 1 2 0 1 2 0 ...
  
  filtered <- dplyr::filter(cds_frames, Frame == 0L)
  # > str(cds_frames)
  # tibble [621 x 2] (S3: tbl_df/tbl/data.frame)
  # $ Count: num [1:621] 811 460 2978 429 251 ...
  # $ Frame: num [1:621] 0 1 2 0 1 2 0 1 2 0 ...
  
  filtered_counts <- filtered$Count
  # > str(filtered_counts)
  # num [1:207] 811 429 488 102 994 146 173 762 13 176 ...
  
}

#' GetAllCodonPosCounts(): extracts A-site assigned counts for a list of genes 
#' 
#' Applies the GetGeneCodonPosReads1dsnap() function to a list of genes and generates 
#' a tidy data frame (tibble) which contains the counts for all genes 
#' 
#' @param gene from gene_names to get read lengths for
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' 
#' @return a tibble which contains the columns "Gene", "PosCodon" and "Count" for a list of genes
#' 


GetAllCodonPosCounts <- function(gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame){
  
  gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
  
  GetAllCodonPosCounts1Gene <- function(gene, dataset, hd_file, min_read_length, asite_displacement_length, snapdisp, filter_for_frame){
    
    subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
    
    left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start))
    
    right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end))
    
    asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", 
                                                                                "yeast_standard_asite_disp_length.txt"))
      
    if(filter_for_frame == FALSE){
      
      codon_counts_1_gene <- GetGeneCodonPosReads1dsnap(gene, dataset, hd_file, 
                                                        left, right, min_read_length, 
                                                        asite_displacement_length, snapdisp)
      
      # > str(codon_counts_1_gene)
      # num [1:207] 4249 825 1017 1176 1116 ...
      
    } else {
      
      codon_counts_1_gene <- FilterForFrameFunction(gene,dataset, hd_file, asite_displacement_length, left, right)
      
    }
    
    
    codon_pos_counts <- tibble(Gene = gene,
                               PosCodon = 1:length(codon_counts_1_gene),
                               Count = codon_counts_1_gene)
    
    as.data.frame(codon_pos_counts, row.names = NULL, optional = FALSE)
    
    return(codon_pos_counts)
    
  }
  
  total_codon_pos_counts <- purrr::map_dfr(.x = gene_names,
                                           .f = GetAllCodonPosCounts1Gene,
                                           dataset,
                                           hd_file,
                                           min_read_length,
                                           snapdisp = snapdisp,
                                           filter_for_frame = filter_for_frame
  )
  
  return (total_codon_pos_counts)
}

total_codon_pos_counts <- suppressMessages(GetAllCodonPosCounts(gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame))

#TEST: GetAllCodonPosCounts(): returns a tibble. 
#TEST: GetAllCodonPosCounts(): the tibble has 3 columns.
#TEST: GetAllCodonPosCounts(): the column names are %in% c("Gene", "PosCodon" and "Count").
#TEST: GetAllCodonPosCounts(): number of observations in the output tibble = sum of CDS (codon co-ordinates) for all genes in gene_names.
#TEST: GetAllCodonPosCounts(): the unique gene names in column "Gene" match the genes in gene_names (unique(total_codon_pos_counts$Gene) = gene_names) = TRUE.
#gives: 
# > str(total_codon_pos_counts)
# Classes 'tbl_df', 'tbl' and 'data.frame':   2749 observations of 3 variables:
#   $ Gene    : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ PosCodon: int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count   : num  4249 825 1017 1176 1116 ...
# 
# 
# #TEST:: GetAllCodonPosCounts example with tinysim, filter_for_frame = FALSE
# 
# 
# Gene      PosCodon  Count
# 1 MAT       1         0
# 2 MAT       2         2
# 3 MAT       3         2
# 4 MAT       4         0
# 5 MIKE      1         0
# 6 MIKE      2         1
# 7 MIKE      3         0
# 8 MIKE      4         0
# 9 MIKE      5         0
#
## TEST:: GetAllCodonPosCounts example with tinysim, filter_for_frame = TRUE
# 
# A tibble: 9 x 3
# Gene  PosCodon Count
# <chr>    <int> <dbl>
#   1 MAT          1     0
# 2 MAT          2     0
# 3 MAT          3     2
# 4 MAT          4     0
# 5 MIKE         1     0
# 6 MIKE         2     0
# 7 MIKE         3     0
# 8 MIKE         4     0
# 9 MIKE         5     0
# 
# this is expected as only reads mapping to the first nucleotide of the codon are retained
    
AddCodonNamesToCodonPosCounts <- function(yeast_codon_pos_i200, gene_names, dataset, hd_file, min_read_length, colsum_out, gff_df, filter_for_frame, snapdisp){
  
  total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame)
  
  transcript_tibbles <- left_join(total_codon_pos_counts, yeast_codon_pos_i200, by = c("PosCodon", "Gene"), keep = FALSE, copy = TRUE)
  
  transcript_gene_pos_poscodon_frame <- tibble(
    Gene = transcript_tibbles$Gene,
    PosCodon = transcript_tibbles$PosCodon,
    Count = transcript_tibbles$Count,
    Codon = transcript_tibbles$Codon
  )
  
  return(transcript_gene_pos_poscodon_frame)
}   


transcript_gene_pos_poscodon_frame <- suppressMessages(AddCodonNamesToCodonPosCounts(yeast_codon_pos_i200, gene_names, dataset, hd_file, min_read_length, colsum_out, gff_df, filter_for_frame, snapdisp))


##TEST: Expect to produce a tibble with each position in the CDS having the correct codon beside it.

#TEST: AddCodonNamesToCodonPosCounts(): returns a tibble. 
#TEST: AddCodonNamesToCodonPosCounts(): the tibble has 4 columns.
#TEST: AddCodonNamesToCodonPosCounts(): the column names are %in% c("Gene", "PosCodon", "Count", "Codon").
#TEST: AddCodonNamesToCodonPosCounts(): number of observations in the output tibble = sum of CDS (codon co-ordinates) for all genes in gene_names.
#TEST: AddCodonNamesToCodonPosCounts(): the unique gene names in column "Gene" match the genes in gene_names (unique(total_codon_pos_counts$Gene) = gene_names) = TRUE.
#gives: 
#Example: using Mok-tinysim data, the following tibble is returned 
# # A tibble: 9 x 4
# Gene  PosCodon Count Codon
# <chr>    <dbl> <dbl> <chr>
#   1 MAT          1     0 ATG  
# 2 MAT          2     2 GCC  
# 3 MAT          3     2 ACA  
# 4 MAT          4     0 TGA  
# 5 MIKE         1     0 ATG  
# 6 MIKE         2     1 ATC  
# 7 MIKE         3     0 AAG  
# 8 MIKE         4     0 GAG  
# 9 MIKE         5     0 TAA  

        
FilterForFeatureOfInterestPositions <- function(yeast_codon_pos_i200, gene, gene_names, dataset, hd_file, min_read_length, colsum_out, gff_df, feature_of_interest, filter_for_frame, snapdisp ){
  
  transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(yeast_codon_pos_i200, 
                                                                      gene_names, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      colsum_out, 
                                                                      gff_df,
                                                                      filter_for_frame,
                                                                      snapdisp)
  
  interesting_feature_table <- dplyr::filter(transcript_gene_pos_poscodon_frame, Codon == feature_of_interest & Gene == gene)
  
  return(interesting_feature_table)
}

# gene <- 'MAT'

# interesting_feature_table <- FilterForFeatureOfInterestPositions(yeast_codon_pos_i200, gene, gene_names, dataset, hd_file, min_read_length, colsum_out, gff_df, feature_of_interest, filter_for_frame, snapdisp)

# When run individually FilterForFeatureOfInterestPositions will only run on an individual gene, so needs a gene name as an input.
#TEST: FilterForFeatureOfInterestPositions(): returns a tibble if there is at least one occurrence of the feature_of_interest.
#TEST: FilterForFeatureOfInterestPositions(): returns an empty tibble if there are no occurrences of the feature_of_interest 
#TEST: FilterForFeatureOfInterestPositions(): the tibble has 4 columns.
#TEST: FilterForFeatureOfInterestPositions(): the column names are %in% c("Gene", "PosCodon", "Count", "Codon").
#TEST: FilterForFeatureOfInterestPositions(): number of observations in the output tibble = number of occurances of the feature of interest in the gene being studied .
#TEST: FilterForFeatureOfInterestPositions(): the unique gene names in column "Gene" match the genes in gene_names (unique(total_codon_pos_counts$Gene) = gene_names) = TRUE.
#gives: 
#Example: using Mok-tinysim data, when feature_of_interest is GCC and gene is MAT, the following tibble is returned 
# # A tibble: 1 x 4
# Gene  PosCodon Count Codon
# <chr>    <dbl> <dbl> <chr>
#   1 MAT          2     2 GCC


## Slice out interesting features

# ExpandFeatureRegionAllGenes takes a slice out of transcript_gene_pos_poscodon_frame, centered on the position of the feature of interest
# It then sets the position of the feature of interest to 0, and changes the positions of adjacent codons to be relative to the feature of interest
# This is done for all occurrences of the feature of interest on all genes in the sample provided.

ExpandFeatureRegionAllGenes <- function(yeast_codon_pos_i200, 
                                gene_names, dataset, hd_file, 
                                min_read_length, colsum_out, 
                                gff_df, feature_of_interest,
                                expand_width, 
                                remove_overhang,
                                filter_for_frame,
                                snapdisp) {
  
  transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(yeast_codon_pos_i200, 
                                                                      gene_names, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      colsum_out, 
                                                                      gff_df,
                                                                      filter_for_frame,
                                                                      snapdisp)

 # take as inputs and select for positions on separate genes
  AllGeneInterestingFeatures <- function(yeast_codon_pos_i200, 
                                         gene, gene_names, dataset, hd_file, 
                                         min_read_length, colsum_out, 
                                         gff_df, feature_of_interest, transcript_gene_pos_poscodon_frame,
                                         filter_for_frame){ 
    print(paste('Checking',gene))
    
    TranscriptForOneGene <- function(yeast_codon_pos_i200, 
                                     gene, gene_names, dataset, hd_file, 
                                     min_read_length, colsum_out, 
                                     gff_df, feature_of_interest){
      
      interesting_feature_tibble <- dplyr::filter(transcript_gene_pos_poscodon_frame, Codon == feature_of_interest)
      
      transcript_for_one_gene <- dplyr::filter(interesting_feature_tibble, Gene == gene)
      
      return(transcript_for_one_gene)
      
    }
    
    transcript_for_one_gene <- TranscriptForOneGene(yeast_codon_pos_i200, 
                                                 gene, gene_names, dataset, hd_file, 
                                                 min_read_length, colsum_out, 
                                                 gff_df, feature_of_interest)
    
    
    #if (remove_overhang) {
    # return an empty tibble if the desired region hangs over the edge of the coding region
    
    ExpandRegions <- function(transcript_for_one_gene, transcript_gene_pos_poscodon_frame, gene, dataset, hd_file, expand_width, remove_overhang = TRUE){
      
      interesting_features <- transcript_for_one_gene
      
      transcript_gene_pos_poscodon_gene_interest <- dplyr::filter(transcript_gene_pos_poscodon_frame, Gene == gene)
      
      gene_length <- filter(gff_df, gff_df$type == 'CDS' & gff_df$Name == gene)$width
      
      if (interesting_features <= expand_width  |interesting_features + expand_width > gene_length/3) {
        return()
      } else {
        output_feature_info <- tibble(
          dplyr::slice(transcript_gene_pos_poscodon_gene_interest, (interesting_features - expand_width):(interesting_features + expand_width), each = FALSE),
          Rel_Pos =  seq(- expand_width, expand_width)
        )
        
        if(dim(output_feature_info)[1] == (2*expand_width + 1)){
          
          return(output_feature_info)
        }else{
          return()
        }
      }
    }
    # The if statement ensures that feature positions that are less/more than the 
    # expand_width value are discarded 
    output_feature_info <- purrr::map(.x = transcript_for_one_gene$PosCodon, .f = ExpandRegions, 
                                      transcript_gene_pos_poscodon_frame,
                                      gene,
                                      dataset,
                                      hd_file,
                                      expand_width,
                                      remove_overhang = TRUE)
    
    
    return(output_feature_info)
  }
  
  output_feature_info <- purrr::map(.x = gene_names, .f = AllGeneInterestingFeatures,
                                    gene_names = gene_names,
                                    yeast_codon_pos_i200 = yeast_codon_pos_i200, 
                                    dataset,
                                    hd_file, 
                                    min_read_length,
                                    colsum_out, 
                                    gff_df,
                                    feature_of_interest,
                                    transcript_gene_pos_poscodon_frame)
  
  # produces a list for each gene, containing a list for each occurrence of the feature of interest
  # Unlist to produce one list, containing each occurrence of the feature of interest
  
  output_feature_info <- unlist(output_feature_info, recursive = F)
  
  # remove NULLS, which represent features of interest occuring within one expand width of the UTRs
  output_feature_info <- output_feature_info[!sapply(output_feature_info, is.null)]
  
  return(output_feature_info)
}


output_feature_info <- suppressMessages(ExpandFeatureRegionAllGenes(yeast_codon_pos_i200 = yeast_codon_pos_i200, 
                                           gene_names = gene_names, dataset, hd_file, 
                                           min_read_length, colsum_out, 
                                           gff_df, feature_of_interest, 
                                           expand_width, remove_overhang, filter_for_frame, snapdisp))


#TEST: ExpandFeatureRegionAllGenes(): creates an object of type "list" 
#TEST: ExpandFeatureRegionAllGenes(): Returns an empty list if there are no occurrences of the feature_of_interest 
#TEST: ExpandFeatureRegionAllGenes(): length(list) == nrow(interesting_feature_table)
#TEST: ExpandFeatureRegionAllGenes(): the tibble contains 5 columns = TRUE
#TEST: ExpandFeatureRegionAllGenes(): the column names are %in% c("Gene", "PosCodon", "Count", "Codon, "Rel_Pos") 
#TEST: ExpandFeatureRegionAllGenes(): number of observations in the output tibble = "expand_width" * 2 + 1, so if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandFeatureRegionAllGenes(): the position from "interesting_feature_positions" has "Rel_Pos" value 0 = TRUE
#TEST: ExpandFeatureRegionAllGenes(): the column "Rel_Pos" goes from -"expand_width to +"expand_width"
## Example, tidysim with feature_of_interest == 'GCC' and expand_width == '1L'
# [[1]]
# # A tibble: 3 x 5
# Gene  PosCodon Count Codon Rel_Pos
# <chr>    <dbl> <dbl> <chr>   <int>
#   1 MAT          1     0 ATG        -1
# 2 MAT          2     2 GCC         0
# 3 MAT          3     2 ACA         1




if(length(output_feature_info) == 0){
  print('No occurrances of the feature of interest')
  
  if(expand_width>1){
    
    print('Try script with an expand_width of 1L to check for occurances near to start or stop codon')
  }
  
  print('Done')
  stop()
}



### Normalization ###

# Normalization carried out within each expanded frame so that they are comparable 
# Normalizes the expand_feature_region list generating a RelCount column with the normalization values

print('Normalising data')

#' ExpandedRegionNormalization(): carries out normalization within each expanded frame 
#' 
#' Normalizes the ExpandFeatureRegion list generating a RelCount column with the normalization values.
#' As this function is not looped it will only generate one normalized tibble for each occurrence of the feature of interest 
#' 
#' @param .x which is the list of tidy format data frames (tibbles) generated by the function ExpandFeatureRegion
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' 
#' @return the list of tibbles which contain the normalized counts within the window so that the feature is comparable despite overall varying levels of expression between genes 
#' 

 
 

ExpandedRegionNormalization <- function(.x, expand_width){
  
  # dplyr::mutate(.x, RelCount = PerCodonCounts / sum(PerCodonCounts) * (2 * expand_width + 1))
  normalized_expand_tibble <- dplyr::mutate(.x, RelCount = Count / sum(Count) * (2 * expand_width + 1))
  
  CheckForNaN <- function(normalized_expand_tibble){
    
    Relcount_values <- unlist(normalized_expand_tibble$RelCount)
    
    SetNaNToZero <- function(Relcount_values){
      
      if(is.nan(Relcount_values)){
        Relcount_values <- 0
      }else{
        Relcount_values <- Relcount_values 
      }
    }
    Relcount_values <- unlist(purrr::map(.x = Relcount_values, .f = SetNaNToZero))
    
    dplyr::mutate(.x, RelCount = Relcount_values)
  }
  CheckForNaN(normalized_expand_tibble)
}

#TEST: ExpandedRegionNormalization(): creates a tidy format data frame (tibble) = TRUE
#TEST: ExpandedRegionNormalization(): the tibble contains 5 columns = TRUE
#TEST: ExpandedRegionNormalization(): the column names are %in% c("Gene", "Pos_Codon", "Rel_Count", "Rel_Pos", "RelCount")
#TEST: ExpandedRegionNormalization(): number of observations in the output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandedRegionNormalization(): the column "Rel_Pos" goes from -"expand_width to +"expand_width" 
#TEST: ExpandedRegionNormalization(): sum(normalized_expand_list[[1]]$RelCount)/nrow(normalized_expand_list[[1]]) == 1 
#TEST: ExpandedRegionNormalizetion(): None of the RelCount columns should contain NaN. These should all be set to 0
#
# Example: using tinysim
# A tibble: 3 x 6
# Gene  PosCodon Count Codon Rel_Pos RelCount
# <chr>    <dbl> <dbl> <chr>   <int>    <dbl>
#   1 MAT          1     0 ATG        -1      0  
# 2 MAT          2     2 GCC         0      1.5
# 3 MAT          3     2 ACA         1      1.5
#



# Normalization carried out for all the tibbles within ExpandList 
normalized_expand_list <- purrr::map(
  .x = output_feature_info,
  .f = ExpandedRegionNormalization,
  expand_width
)


# TEST:: normalised_expand_list should be of type list
# type(normalized_expand_list)
# [1] "list"

# TEST:: Normalized_expand_list should be the same length as expand_feature_region
# > length(normalized_expand_list)==length(expand_feature_region)
# [1] TRUE

# TEST:: the dimensions of each item in the list shoud be [(2*expand_width+1) X 5] as there are now 5 rows; Gene, Pos_Codon, Rel_Count, Rel_Pos, RelCount
# > dim(expand_feature_region[[1]])
# [1] 11  5

# TEST:: At each position, the sum of RelCount should be equal to (2*expand_width+1)
# ie if the expand width was 5:
# sum(normalized_expand_list[[1]]$RelCount)
# [1] 11


### Overlaying the normalized expanded tibbles ###

print('Calculating average')

# Function to overlay graphs into a single graph. Need to generate a single tibble 
# from NormalizedExpandList. Join by Rel_Pos, in RelCount need the mean for 
# each Rel_Pos (sum row(x) / number of row(x))


#' OverlayedTable(): overlays tidy format data frames (tibbles) to create a single overlayed tibble.
#' 
#' FIXME: Takes normalized_expand_list as its input.
#' 
#' @param normalized_expand_list the output from the looped function ExpandedRegionNormalization()
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' 
#' @return a tibble which contains the mean counts for each position from the normalized tibbles. 
#' 
#' @example 
#' 
#' normalized_expand_list <- purrr::map(.x = expand_feature_region,
#'                                      .f = ExpandedRegionNormalization,
#'                                      expand_width)
#' 
#' OverlayedTable(normalized_expand_list, expand_width)
#' 
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
#TEST: OverlayedTable(): creates a tidy format data frame (tibble) = TRUE
#TEST: OverlayedTable(): the tibble contains 2 columns = TRUE
#TEST: OverlayedTable(): the column names are %in% c("Rel_Pos", "RelCount")
#TEST: OverlayedTable(): number of observations in the output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: OverlayedTable(): the column "Rel_Pos" goes from -"expand_width to +"expand_width" 
#TEST: OverlayedTable(): RelCount is a numeric 
#gives:
# > str(overlayed_tibbles)
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 2 variables
#   $ Rel_Pos : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
#   $ RelCount: num  0.893 1.125 0.992 0.998 0.779 ...

overlayed_tibbles <- OverlayedTable(normalized_expand_list, expand_width) 

print('Creating plot')

overlayed_plot <- ggplot(overlayed_tibbles, mapping = aes(x = Rel_Pos, y = RelCount)) + 
  geom_line() +
  theme_bw()+
  theme(text=element_text(size=14),
        axis.title=element_text(size=14, face='bold'),
        title = element_text(size = 14, face='bold'))+
  labs(title = paste0('Relative read counts around feature ', feature_of_interest),
       x = 'Position relative to feature of interest',
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
