rm(list=ls())

# Given an h5 file, GFF file and .tsv file, this script creates a metafeature plot 
# for the feature of interest. 

### TEST ###

# run the script on the tinysim dataset. Documentation of tests in issue ticket 
# #402 "Writing tests for visualisation of feature scripts" 
# (https://github.com/riboviz/riboviz/issues/402)

print('Starting process')

# source packages and functions from rscripts 
suppressMessages(source(here::here("rscripts", "read_count_functions.R")))
suppressMessages(source(here::here("rscripts", "stats_figs_block_functions.R")))

# read in dependent packages
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))


# Set optparse arguments

option_list <- list(make_option(c('-i', '--input'), type = "character", help='Path input to h5 file'),
                    make_option(c('-d', '--dataset'), type = "character", help='Name of dataset being studied'),
                    make_option(c('-g', '--gff'), type = "character", help='Path to the GFF3 file of the organism being studied'),
                    make_option(c('-a', '--annotation'), type = "character", help='Path to codon positions table for organism being studied'),
                    make_option(c('--feature'), type = "character", help='Feature of interest, e.g. codon pair'),
                    make_option(c('-o', '--output'), type = "character", help='Path to output directory'),
                    make_option(c('--expand_width'), type = "integer", help='the desired range either side of the feature of interest', default = 5),
                    make_option(c('--frame'), type = "integer", help='Reading frame to filter counts for, either 0, 1 or 2', default = 0),
                    make_option(c('--minreadlen'), type = "integer", help='Minimum read length', default = 10),
                    make_option(c('--filter_for_frame'), type = "logical", help='FALSE: keep counts for all reading frames, TRUE: filter for a reading frame', default = TRUE),
                    make_option(c('--snapdisp'), type = "integer", help='Frame to filter for wjem using SnaptToCodon', default = 0L))

opt <- optparse::parse_args(OptionParser(option_list = option_list))

# parser$add_argument('--colsum_out', help='logical', default = TRUE)


hd_file <- opt$input
dataset <- opt$dataset
gff <- opt$gff
yeast_codon_table <- opt$annotation
feature_of_interest <- opt$feature
output_dir <- opt$output
expand_width <- opt$expand_width
filter_for_frame <- opt$frame
min_read_length <- opt$minreadlen
# colsum_out <- opt$colsum_out
snapdisp <- opt$snapdisp


# hd_file <- here::here("Mok-simYAL5", "output", "A", "A.h5")
# dataset <- "Mok-simYAL5"
# feature_of_interest <- 'CGA TAG'
# expand_width = 5L
# filter_for_frame <- TRUE
# min_read_length <- 10
yeast_codon_table <- here::here("data", "yeast_codon_table.tsv")
# gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3")
# colsum_out <- TRUE
# snapdisp <- 0L
# output_dir <- here::here(".")


# hd_file <- here::here("Mok-tinysim", "A", "A.h5")
# dataset <- "Mok-tinysim"
# gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "tiny_2genes_20utrs.gff3")
# gene <- "MAT"
# gene_names <- unique(gff_df$Name)
# yeast_codon_table <- here::here("data", "tinysim_codon_table.tsv")
# 
# 
# 
# hd_file <- here::here("J-Sc_2014", "sec63_1.h5")
# dataset <- "J-Sc_2014"
# gff <- here::here("..", "example-datasets", "fungi", "saccharomyces", "annotation", "Saccharomyces_cerevisiae_yeast_CDS_w_250utrs.gff3")
# gene_names <- c("YOL012C", "YOL039W", "YOL086C")

gff_df <- readGFFAsDf(gff)

gene_names <- rhdf5::h5ls(hd_file, 
                          recursive = 1)$name

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
  CodonPos1 = yeast_codon_pos_i200$PosCodon, 
  CodonPos2 = dplyr::lead(yeast_codon_pos_i200$PosCodon),
  CodonPair = paste(yeast_codon_pos_i200$Codon, (dplyr::lead(yeast_codon_pos_i200$Codon) 
                                                 %>% str_replace_all("ATG", "NA")))
)

    # > str(gene_poscodon_codon_i200)
    # tibble [2,826,757 x 4] (S3: tbl_df/tbl/data.frame)
    # $ Gene      : chr [1:2826757] "YAL068C" "YAL068C" "YAL068C" "YAL068C" ...
    # $ CodonPos_1: num [1:2826757] 1 2 3 4 5 6 7 8 9 10 ...
    # $ CodonPos_2: num [1:2826757] 2 3 4 5 6 7 8 9 10 11 ...
    # $ CodonPair : chr [1:2826757] "ATG GTC" "GTC AAA" "AAA TTA" "TTA ACT" ...


#####



#' FilterForFrame(): fetches A-site assigned counts, filters for a reading frame
#' 
#' Counts are fetched with GetGeneDatamatrix() and A-site assignment is carried out with CalcAsiteFixed() (in nt positions)
#' the reading frame is then filtered for and the counts are returned (in codon positions)
#' 
#' @param gene from gene_names to get read lengths for
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping, default = 0L
#' 
#' @return a list of numeric values (counts) for a single reading frame (0, 1 or 2)
#' 
#' @example 
#' 
#' FilterForFrame(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file, min_read_length = 10, snapdisp = 0L)
#' 
#' @export
FilterForFrame <- function(gene, dataset, hd_file, min_read_length, snapdisp){
  
  asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", 
                                                                              "yeast_standard_asite_disp_length.txt"))
  
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
  
  reads_asitepos <- CalcAsiteFixed(reads_pos_length, 
                                   min_read_length, 
                                   asite_displacement_length)
  
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
    # where gene = YAL003W
  
  left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") 
                     %>%  select(start))
    # 251
  
  right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") 
                      %>%  select(end))
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
  
  filtered <- dplyr::filter(cds_frames, Frame == snapdisp)
    # > str(cds_frames)
    # tibble [621 x 2] (S3: tbl_df/tbl/data.frame)
    # $ Count: num [1:621] 811 460 2978 429 251 ...
    # $ Frame: num [1:621] 0 1 2 0 1 2 0 1 2 0 ...
  
  filtered_counts <- filtered$Count
    # > str(filtered_counts)
    # num [1:207] 811 429 488 102 994 146 173 762 13 176 ...
  
}
# TEST: FilterForFrame(): returns a list of numeric values = TRUE
# TEST: FilterForFrame(): length(filtered_counts) = length(CDS) for gene of interest = TRUE
# gives:
# > str(filtered_counts)
# num [1:207] 811 429 488 102 994 146 173 762 13 176 ...



### Fetch and format counts from the h5 file ###



#' GetAllCodonPosCounts(): extracts A-site assigned counts for a list of genes 
#' 
#' This function extracts the A-site assigned counts and generates a tidy data frame
#' (tibble) which contains the counts for all genes in the list of genes.
#' 
#' If filter_for_frame = FALSE the GetGeneCodonPosReads1dsnap() function is applied to a 
#' list of genes and generates a tidy data frame (tibble) which contains the counts from all reading frames.
#' 
#' If filter_for_frame = TRUE the FilterForFrame() function is applied instead, this generates
#' a tibble which contains the read counts for a reading frame (snapdisp = 0L or 1L or 2L). 
#' 
#' @param gene from gene_names to get read lengths for
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping, default = 0L
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if keeping and grouping all reading frames for each codon
#' 
#' @return a tibble which contains the columns "Gene", "PosCodon" and "Count" for a list of genes
#' 
#' @example 
#' 
#' GetAllCodonPosCounts(gene_names, 
#'                      dataset = "Mok-simYAL5", 
#'                      hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), 
#'                      min_read_length = 10, 
#'                      snapdisp = 0L, 
#'                      filter_for_frame = TRUE)
#' 
#' @export 
GetAllCodonPosCounts <- function(gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame){
  
  gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
  
  GetAllCodonPosCounts1Gene <- function(gene, dataset, hd_file, min_read_length, 
                                        asite_displacement_length, snapdisp, 
                                        filter_for_frame){
    
    subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
    
    left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") 
                       %>%  select(start))
    
    right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") 
                        %>%  select(end))
    
    asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", 
                                                                                "yeast_standard_asite_disp_length.txt"))
    
    reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
    
    reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length, 
                                     asite_displacement_length)
    
    
    if(filter_for_frame == FALSE){
      
      codon_counts_1_gene <- GetGeneCodonPosReads1dsnap(gene, dataset, hd_file, 
                                                        left, right, min_read_length, 
                                                        asite_displacement_length, 
                                                        snapdisp = snapdisp)
      
        # > str(codon_counts_1_gene)
        # num [1:207] 4249 825 1017 1176 1116 ...
      
    } else {
      
      codon_counts_1_gene <- FilterForFrame(gene, dataset, hd_file, 
                                                    min_read_length, snapdisp)
      
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
#TEST: GetAllCodonPosCounts(): returns a tibble. 
#TEST: GetAllCodonPosCounts(): the tibble has 3 columns.
#TEST: GetAllCodonPosCounts(): the column names are %in% c("Gene", "PosCodon" and "Count").
#TEST: GetAllCodonPosCounts(): number of observations in the output tibble = sum of CDS (codon co-ordinates) for all genes in gene_names.
#TEST: GetAllCodonPosCounts(): the unique gene names in column "Gene" match the genes in gene_names (unique(total_codon_pos_counts$Gene) = gene_names) = TRUE.
#gives: 
# for filter_for_frame = FALSE
# > str(total_codon_pos_counts)
# Classes 'tbl_df', 'tbl' and 'data.frame':   2749 observations of 3 variables:
#   $ Gene    : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ PosCodon: int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count   : num  4249 825 1017 1176 1116 ...
# for filter_for_frame = TRUE, snapdisp = 0L
# > str(total_codon_pos_counts)
# Classes 'tbl_df', 'tbl' and 'data.frame':   2749 observations of 3 variables:
#   $ Gene    : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ PosCodon: int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count   : num  811 429 488 102 994 146 173 762 13 176 ...


 

#' AddCodonNamesToCodonPosCounts(): takes codon names from the annotation .tsv file, joined to the codon counts table 
#' 
#' Uses the function GetAllCodonPosCounts().
#' 
#' @return a  tidy format data frame (tibble) which contains the genes in gene_names, codon positions, counts and the codon pair. 
#' 
#' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
#' @param gene from gene_names  
#' @param dataset name of dataset stored in .h5 file. 
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if keeping and grouping all reading frames for each codon
#' 
#' @example 
#' 
#' AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200,
#'                               gene_names,
#'                               dataset = "Mok-simYAL5",
#'                               hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"),
#'                               min_read_length = 10,
#'                               snapdisp = 0L,
#'                               filter_for_frame = TRUE)
#' 
#' @export
AddCodonNamesToCodonPosCounts <- function(gene_poscodon_codon_i200, gene_names, 
                                          dataset, hd_file, min_read_length, 
                                          snapdisp, filter_for_frame){
  
  total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, dataset, hd_file, 
                                                 min_read_length, snapdisp, 
                                                 filter_for_frame = filter_for_frame)
  
  transcript_tibbles <- total_codon_pos_counts %>% left_join(gene_poscodon_codon_i200, 
                                                             by = c("PosCodon" = "CodonPos1", 
                                                                    "Gene" = "Gene"), 
                                                             keep = FALSE, 
                                                             copy = TRUE)
  
  transcript_gene_pos_poscodon_frame <- tibble(
    Gene = transcript_tibbles$Gene,
    CodonPos1 = transcript_tibbles$PosCodon,
    CodonPos2 = transcript_tibbles$CodonPos2,
    Count = transcript_tibbles$Count,
    CodonPair = transcript_tibbles$CodonPair
  )
  
  return(transcript_gene_pos_poscodon_frame)
}
#TEST: AddCodonNamesToCodonPosCounts(): creates a tibble = TRUE
#TEST: AddCodonNamesToCodonPosCounts(): the tibble contains columns = TRUE
#TEST: AddCodonNamesToCodonPosCounts(): number of observations in the output tibble = sum of CDS (codon co-ordinates) for all genes in gene_names.
#TEST: AddCodonNamesToCodonPosCounts(): the column names are %in% c("Gene", "CodonPos1", "CodonPos2", "Count", "CodonPair") 
#TEST: AddCodonNamesToCodonPosCounts(): the unique gene names in column "Gene" match the genes in gene_names (unique(total_codon_pos_counts$Gene) = gene_names) = TRUE.
# gives:
# > str(transcript_gene_pos_poscodon_frame)
# Classes 'tbl_df', 'tbl' and 'data.frame':   2,749 observations of 5 variables:
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  1 2 3 4 5 6 7 8 9 10 ...
#   $ CodonPos2: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Count    : num  4249 825 1017 1176 1116 ...
#   $ CodonPair: chr  "ATG GCA" "GCA TCC" "TCC ACC" "ACC GAT" ...



### Filtering for feature of interest ###



#' FilterForFeatureOfInterestPositions(): Filters for the feature of interest 
#' 
#' Uses the AddCodonNamesToCodonPosCounts to get counts 
#' 
#' @return a tidy format data frame (tibble) which contains the gene, codon positions and counts for the feature of interest 
#' 
#' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
#' @param gene from gene_names  
#' @param dataset name of dataset stored in .h5 file. 
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if keeping and grouping all reading frames for each codon
#' @param feature_of_interest character, each incidence of the feature will be extracted from transcript_info_tibble
#' 
#' @example 
#' 
#' FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, 
#'                                     gene, 
#'                                     dataset = "Mok-simYAL5", 
#'                                     hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), 
#'                                     min_read_length = 10, 
#'                                     snapdisp = 0L, 
#'                                     filter_for_frame = TRUE, 
#'                                     feature_of_interest = 'TCC AAG')
#' 
#' @export
FilterForFeatureOfInterestPositions <- function(gene_poscodon_codon_i200, gene_names, 
                                                dataset, hd_file, min_read_length, 
                                                snapdisp, filter_for_frame, feature_of_interest){
  
  transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, 
                                                                      gene_names, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      snapdisp, 
                                                                      filter_for_frame)
  
  interesting_feature_table <- dplyr::filter(transcript_gene_pos_poscodon_frame, 
                                             CodonPair == feature_of_interest)
  
  return(interesting_feature_table)
}
#TEST: FilterForFeatureOfInterestPositions(): returns a tibble = TRUE
#TEST: FilterForFeatureOfInterestPositions(): the tibble contains 5 columns = TRUE
#TEST: FilterForFeatureOfInterestPositions(): the column names are %in% c("Gene", "CodonPos1", "CodonPos2", "Count", "CodonPair") 
#TEST: FilterForFeatureOfInterestPositions(): the unique feature of interest in column "CodonPair" (unique(interesting_feature_positions$CodonPair)) = 1
#gives:
# > str(interesting_feature_positions)
# Classes 'tbl_df', 'tbl' and 'data.frame':   8 observations of 5 variables:
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL005C" "YAL005C" ...
#   $ CodonPos1: num  7 57 383 508 535 321 90 412
#   $ CodonPos2: num  8 58 384 509 536 322 91 413
#   $ Count    : num  1553 1147 358 317 498 ...
#   $ CodonPair: chr  "TCC AAG" "TCC AAG" "TCC AAG" "TCC AAG" ...



### Expand frame around feature of interest ###



#' ExpandFeatureRegion(): Generates a window around each feature of interest across all genes in gene_names
#' 
#' @return a list of tibbles, a tibble for each feature of interest with an expanded window around the feature
#' 
#' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
#' @param gene from gene_names  
#' @param dataset name of dataset stored in .h5 file. 
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param feature_of_interest character, each incidence of the feature will be extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' @param remove_overhang logical, default = TRUE, removes features of interest who's positions fall within expand_width and therefore do not allow for complete windows to be generated
#' @param snapdisp integer any additional displacement in the snapping
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if keeping and grouping all reading frames for each codon
#' 
#' @example 
#' 
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' ExpandFeatureRegion(gene_poscodon_codon_i200, 
#'                     gene_names, 
#'                     dataset = "Mok-simYAL5", 
#'                     hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), 
#'                     min_read_length = 10, 
#'                     gff_df, 
#'                     feature_of_interest = 'TCC AAG', 
#'                     expand_width = 5L,
#'                     remove_overhang = TRUE, 
#'                     snapdisp = 0L, 
#'                     filter_for_frame = TRUE)
#' 
#' @export                                                           
ExpandFeatureRegion <- function(gene_poscodon_codon_i200, gene_names, dataset, 
                                hd_file, min_read_length, gff_df, 
                                feature_of_interest, expand_width, remove_overhang = TRUE,
                                snapdisp, filter_for_frame) {
  
  transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, 
                                                                      gene_names, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      snapdisp, 
                                                                      filter_for_frame)
  
  # take as inputs and select for positions on separate genes
  AllGeneInterestingFeatures <- function(gene_poscodon_codon_i200, 
                                         gene, gene_names, dataset, hd_file, 
                                         min_read_length, 
                                         gff_df, feature_of_interest, 
                                         transcript_gene_pos_poscodon_frame,
                                         filter_for_frame){ 
    print(gene)
    
    TranscriptForOneGene <- function(gene_poscodon_codon_i200, 
                                     gene, gene_names, dataset, hd_file, 
                                     min_read_length, 
                                     gff_df, feature_of_interest){
      
      interesting_feature_tibble <- dplyr::filter(transcript_gene_pos_poscodon_frame, 
                                                  CodonPair == feature_of_interest)
      
      transcript_for_one_gene <- dplyr::filter(interesting_feature_tibble, 
                                               Gene == gene)
      
      return(transcript_for_one_gene)
      
    }
    
    transcript_for_one_gene <- TranscriptForOneGene(gene_poscodon_codon_i200, 
                                                    gene, gene_names, dataset, hd_file, 
                                                    min_read_length, 
                                                    gff_df, feature_of_interest)
    
    
    #if (remove_overhang) {
    # return an empty tibble if the desired region hangs over the edge of the coding region
    
    ExpandRegions <- function(transcript_for_one_gene, transcript_gene_pos_poscodon_frame, 
                              gene, dataset, hd_file, expand_width, remove_overhang = TRUE, gff_df){
      
      interesting_features <- transcript_for_one_gene
      
      transcript_gene_pos_poscodon_gene_interest <- dplyr::filter(transcript_gene_pos_poscodon_frame, 
                                                                  Gene == gene)
      
      subset_gff_df_by_gene <- dplyr::filter(gff_df, seqnames == gene)
      
      gene_length <- dplyr::filter(subset_gff_df_by_gene, type == "CDS") %>% select(width) 
      
      if (interesting_features <= expand_width  |interesting_features + expand_width > gene_length/3) {
        return()
      } else {
        expand_feature_region <- tibble(
          dplyr::slice(transcript_gene_pos_poscodon_gene_interest, 
                       (interesting_features - expand_width):(interesting_features + expand_width), 
                       each = FALSE),
          RelPos =  seq(- expand_width, expand_width)
        )
        
        if(dim(expand_feature_region)[1] == (2*expand_width + 1)){
          
          return(expand_feature_region)
        }else{
          return()
        }
      }
    }
    # The if statement ensures that feature positions that are less/more than the 
    # expand_width value are discarded 
    expand_feature_region <- purrr::map(.x = transcript_for_one_gene$CodonPos1, 
                                        .f = ExpandRegions, 
                                      transcript_gene_pos_poscodon_frame,
                                      gene,
                                      dataset,
                                      hd_file,
                                      expand_width,
                                      remove_overhang = TRUE, 
                                      gff_df)
    
    
    return(expand_feature_region)
  }
  
  expand_feature_region <- purrr::map(.x = gene_names, 
                                      .f = AllGeneInterestingFeatures,
                                      gene_names = gene_names,
                                      gene_poscodon_codon_i200 = gene_poscodon_codon_i200, 
                                      dataset, hd_file, min_read_length, 
                                      gff_df, feature_of_interest,
                                      transcript_gene_pos_poscodon_frame)
  
  expand_feature_region <- unlist(expand_feature_region, 
                                  recursive = F)
  
  # remove NULLS, which represent features of interest occuring within one expand width of the UTRs
  expand_feature_region <- expand_feature_region[!sapply(expand_feature_region, 
                                                         is.null)]
  
  return(expand_feature_region)
}
#' #TEST: ExpandFeatureRegion(): output is a list of tidy format data frames (tibbles) = TRUE. type(expand_feature_region) = "list"
#' #TEST: ExpandFeatureRegion(): number of tibbles in list matches the number of occurrences in "feature_of_interest" list = TRUE
#' #TEST: ExpandFeatureRegion(): each tibble contains 6 columns = TRUE
#' #TEST: ExpandFeatureRegion(): the column names are %in% c("Gene", "Pos_Codon1", "Pos_Codon2", "Count", "CodonPair", "RelPos")
#' #TEST: ExpandFeatureRegion(): number of observations in each output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#' #TEST: ExpandFeatureRegion(): the position from "interesting_feature_positions" has "RelPos" value 0 = TRUE
#' #TEST: ExpandFeatureRegion(): the column "RelPos" goes from -"expand_width to +"expand_width"
#gives:
# > str(expand_feature_region)
# List of 8
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ CodonPos2: num  3 4 5 6 7 8 9 10 11 12 ...
#   $ Count    : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ CodonPair: chr  "GCA TCC" "TCC ACC" "ACC GAT" "GAT TTC" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  52 53 54 55 56 57 58 59 60 61 ...
#   $ CodonPos2: num  53 54 55 56 57 58 59 60 61 62 ...
#   $ Count    : num  42 53 648 293 121 92 519 79 765 196 ...
#   $ CodonPair: chr  "TTC AAC" "AAC CAC" "CAC ATC" "ATC GCT" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
#   $ CodonPos1: num  378 379 380 381 382 383 384 385 386 387 ...
#   $ CodonPos2: num  379 380 381 382 383 384 385 386 387 388 ...
#   $ Count    : num  28 74 56 201 4 166 46 15 34 206 ...
#   $ CodonPair: chr  "ACT GGT" "GGT GAC" "GAC GAA" "GAA TCT" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
#   $ CodonPos1: num  503 504 505 506 507 508 509 510 511 512 ...
#   $ CodonPos2: num  504 505 506 507 508 509 510 511 512 513 ...
#   $ Count    : num  266 468 224 550 42 19 133 31 115 31 ...
#   $ CodonPair: chr  "GAC AAG" "AAG GGT" "GGT AGA" "AGA TTG" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
#   $ CodonPos1: num  530 531 532 533 534 535 536 537 538 539 ...
#   $ CodonPos2: num  531 532 533 534 535 536 537 538 539 540 ...
#   $ Count    : num  21 32 47 3 7 85 123 40 124 39 ...
#   $ CodonPair: chr  "TCT CAA" "CAA AGA" "AGA ATT" "ATT GCT" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL012W" "YAL012W" "YAL012W" "YAL012W" ...
#   $ CodonPos1: num  316 317 318 319 320 321 322 323 324 325 ...
#   $ CodonPos2: num  317 318 319 320 321 322 323 324 325 326 ...
#   $ Count    : num  179 47 260 63 24 29 176 20 33 99 ...
#   $ CodonPair: chr  "GGT GCT" "GCT GAA" "GAA GCT" "GCT GCT" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL035W" "YAL035W" "YAL035W" "YAL035W" ...
#   $ CodonPos1: num  85 86 87 88 89 90 91 92 93 94 ...
#   $ CodonPos2: num  86 87 88 89 90 91 92 93 94 95 ...
#   $ Count    : num  20 57 81 19 39 17 133 144 39 49 ...
#   $ CodonPair: chr  "AAG CCT" "CCT ATA" "ATA CTA" "CTA AAG" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL038W" "YAL038W" "YAL038W" "YAL038W" ...
#   $ CodonPos1: num  407 408 409 410 411 412 413 414 415 416 ...
#   $ CodonPos2: num  408 409 410 411 412 413 414 415 416 417 ...
#   $ Count    : num  185 446 307 37 0 286 271 171 124 161 ...
#   $ CodonPair: chr  "ACC CCA" "CCA AGA" "AGA TTG" "TTG GTT" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...


### Normalisation ###


#' ExpandedRegionNormalisation(): carries out normalization within each expanded frame so that they are comparable
#' 
#' Normalizes the ExpandRegion list generating a RelCount column with the normalization values 
#' As this function is not looped it will only generate one normalized tibble for one occurence 
#' 
#' @return a list of normalised tibbles 
#' 
#' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
#' @param gene from gene_names  
#' @param dataset name of dataset stored in .h5 file. 
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param feature_of_interest character, each incidence of the feature will be extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' @param remove_overhang default = TRUE, removes features of interest who's positions do not allow for complete windows to be generated
#' @param snapdisp integer any additional displacement in the snapping
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if keeping and grouping all reading frames for each codon
#' 
#' @example 
#' 
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' ExpandedRegionNormalisation(gene_poscodon_codon_i200,
#'                             gene_names, 
#'                             dataset = "Mok-simYAL5", 
#'                             hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"),
#'                             min_read_length = 10, 
#'                             gff_df, 
#'                             feature_of_interest = 'TCC AAG',
#'                             expand_width = 5L,
#'                             remove_overhang = TRUE,
#'                             snapdisp = 0L, 
#'                             filter_for_frame = TRUE)
#' 
#' @export 
ExpandedRegionNormalisation <- function(gene_poscodon_codon_i200, 
                                        gene_names, dataset, hd_file, 
                                        min_read_length, colsum_out, 
                                        gff_df, feature_of_interest,
                                        expand_width, 
                                        remove_overhang = TRUE,
                                        snapdisp, filter_for_frame){
  
  expand_feature_region <- ExpandFeatureRegion(gene_poscodon_codon_i200, 
                                               gene_names, dataset, hd_file, 
                                               min_read_length, colsum_out, 
                                               gff_df, feature_of_interest,
                                               expand_width, 
                                               remove_overhang = TRUE,
                                               snapdisp, filter_for_frame)
  
  Normalization <- function(.x, expand_width){
    dplyr::mutate(.x, RelCount = Count / sum(Count) * (2 * expand_width + 1))
 }
  
  
  normalized_expand_list <- purrr::map(
    .x = expand_feature_region,
    .f = Normalization,
    expand_width
  )
  
  CheckForNaN <- function(.x){
    
    normalized_expand_list <- .x 
    
    Relcount_values <- unlist(normalized_expand_list$RelCount)
    
    print(Relcount_values)
    
    SetNaNToZero <- function(Relcount_values){
      
      if(is.nan(Relcount_values)){
        print('NaN present')
        Relcount_values <- 0
      }else{
        Relcount_values <- Relcount_values 
      }
    }
    Relcount_values <- unlist(purrr::map(.x = Relcount_values, 
                                         .f = SetNaNToZero))
    
    print(Relcount_values)
    
    dplyr::mutate(.x, RelCount = Relcount_values)
  }
  
  normalized_expand_list <- purrr::map(.x = normalized_expand_list, 
                                       .f = CheckForNaN)
  
}
#TEST: ExpandedRegionNormalisation(): creates a list of tidy format data frame (tibble) (type(normalized_expand_list)) = TRUE
#TEST: ExpandedRegionNormalisation(): the tibble contains 6 columns = TRUE
#TEST: ExpandedRegionNormalisation(): the column names are %in% c("Gene", "CodonPos1", "CodonPos2", "Count", "RelPos", "RelCount")
#TEST: ExpandedRegionNormalisation(): number of observations in the output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandedRegionNormalisation(): the column "RelPos" goes from -"expand_width to +"expand_width" 
#gives:
# > str(normalized_expanded_feature_region)
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1 : num  2 3 4 5 6 7 8 9 10 11 ...
#   $ CodonPos2 : num  3 4 5 6 7 8 9 10 11 12 ...
#   $ Count     : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ RelPos    : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
#   $ RelCount  : num  1.347 1.532 0.32 3.12 0.458 ...



### Overlaying the normalized expanded tibbles ###



#' OverlayedTibble: overlays the expanded tibbles generated by ExpandedRegionNormalisation
#'
#' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
#' @param gene from gene_names  
#' @param dataset name of dataset stored in .h5 file. 
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param feature_of_interest character, each incidence of the feature will be extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' @param remove_overhang default = TRUE, removes features of interest who's positions do not allow for complete windows to be generated
#' @param snapdisp integer any additional displacement in the snapping
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if keeping and grouping all reading frames for each codon
#' 
#' @return  OverlayedTibble(): overlays tidy format data frames (tibbles) to create a single overlayed tibble.
#' 
#' @example 
#' 
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' OverlayedTibble(gene_poscodon_codon_i200, 
#'                  gene_names, 
#'                  dataset = "Mok-simYAL5",
#'                  hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"),
#'                  min_read_length = 10,
#'                  gff_df,
#'                  feature_of_interest = 'TCC AAG',
#'                  expand_width = 5L,
#'                  remove_overhang = TRUE,
#'                  snapdisp = 0L,
#'                  filter_for_frame = TRUE)
#' 
#' @export
OverlayedTibble <- function(gene_poscodon_codon_i200, 
                           gene_names, dataset, hd_file, 
                           min_read_length, colsum_out, 
                           gff_df, feature_of_interest,
                           expand_width, 
                           remove_overhang = TRUE,
                           snapdisp, filter_for_frame){
  
  normalized_expand_list <- ExpandedRegionNormalisation(gene_poscodon_codon_i200,
                                                        gene_names, dataset, hd_file,
                                                        min_read_length, colsum_out,
                                                        gff_df, feature_of_interest,
                                                        expand_width,
                                                        remove_overhang = TRUE,
                                                        snapdisp, filter_for_frame)
  
  number_of_objects <- length(normalized_expand_list)
  # The number of objects inside normalized_expand_list
  
  result <- lapply(normalized_expand_list, "[", c("RelPos", "RelCount")) 
  # Reduces normalized_expand_list to the columns RelPos and RelCount
  
  joined_result <- purrr::reduce(.x = result, .f =`+`)
  
  joined_rows <- mutate(joined_result / number_of_objects)
  
  overlayed_tibbles <- tibble::tibble(
    RelPos = seq(- expand_width, expand_width),
    RelCount = joined_rows$RelCount
  )
}
#TEST: OverlayedTibble(): creates a tidy format data frame (tibble) = TRUE
#TEST: OverlayedTibble(): the tibble contains 2 columns = TRUE
#TEST: OverlayedTibble(): the column names are %in% c("RelPos", "RelCount")
#TEST: OverlayedTibble(): number of observations in the output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: OverlayedTibble(): the column "RelPos" goes from -"expand_width to +"expand_width" 
#TEST: OverlayedTibble(): RelCount is a numeric 
#gives:
# > str(overlayed_tibbles)
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 2 variables
# $ RelPos : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
# $ RelCount: num [1:11] 1.029 1.07 0.987 1.45 1.151 ...



### Generate plot around feature of interest ###



#' GeneratePlot(): Generates a metafeature plot around the feature of interest
#' 
#' @return A plot which shows the ribosomal occupancy around a feature of interest
#'   
#' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
#' @param gene from gene_names  
#' @param dataset name of dataset stored in .h5 file. 
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param feature_of_interest character, each incidence of the feature will be extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' @param remove_overhang default = TRUE, removes features of interest who's positions do not allow for complete windows to be generated
#' @param snapdisp integer any additional displacement in the snapping
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if keeping and grouping all reading frames for each codon
#' @param size the size of the text on the metafeature plot
#' 
#' @example 
#' 
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' GeneratePlot(gene_poscodon_codon_i200, 
#'              gene_names, 
#'              dataset = "Mok-simYAL5", 
#'              hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), 
#'              min_read_length = 10,
#'              colsum_out, 
#'              gff_df, 
#'              feature_of_interest = 'TCC AAG', 
#'              expand_width = 5L,
#'              remove_overhang = TRUE, 
#'              snapdisp = 0L, 
#'              filter_for_frame = TRUE, 
#'              size = 12)
#'             
#' @export              
GeneratePlot <- function(gene_poscodon_codon_i200, gene_names, dataset, hd_file, 
                         min_read_length, colsum_out, gff_df, feature_of_interest, 
                         expand_width, remove_overhang = TRUE, snapdisp, 
                         filter_for_frame, size = 12){
  
  overlayed_tibbles <- OverlayedTibble(gene_poscodon_codon_i200, gene_names, 
                                       dataset, hd_file, min_read_length, colsum_out, 
                                       gff_df, feature_of_interest, expand_width, 
                                       remove_overhang = TRUE, snapdisp, filter_for_frame)
  
  overlayed_plot <- ggplot(overlayed_tibbles, 
                           mapping = aes(x = RelPos, 
                                         y = RelCount)) + 
    geom_line() +
    theme_bw() +
    theme(text = element_text(size = size),
          axis.title = element_text(size = size),
          title = element_text(size = size, face = 'bold'))+
    labs(title = paste0('Meta-feature plot of codon pair ', 
                        feature_of_interest),
         x = 'Distance from codon pair (3nt)',
         y = 'Normalised ribosomal occupancy', size = 2)
  
  # save_plot_pdf <- function(overlayed_plot, output_dir){
  #   overlayed_plot %>%
  #     ggsave(
  #       filename = file.path(output_dir,"inhibitory_codon_plot.pdf"),
  #       width = 6, height = 5
  #     )
  # }
  # 
  # save_plot_pdf(overlayed_plot, output_dir)
}
#TEST: GeneratePlot(): produces a single plot =TRUE
#TEST: GeneratePlot(): title is "Meta-feature plot of codon pair <feature_of_interest>
#Test: GeneratePlot(): the x-axis is "Distance from codon pair (3nt) and y-axis is "Normalised ribosomal occupancy"
#gives: 


# print('Fetching read counts')
# 
# total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, dataset, hd_file, 
#                                                min_read_length, snapdisp, 
#                                                filter_for_frame)
# 
# print('Creating information tibble')
# 
# transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, 
#                                                                     gene_names, dataset, hd_file, 
#                                                                     min_read_length, snapdisp, 
#                                                                     filter_for_frame)
# 
# print('Filtering for feature of interest')
# 
# interesting_feature_positions <- FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, 
#                                                                      gene_names, 
#                                                                      dataset, 
#                                                                      hd_file, 
#                                                                      min_read_length, 
#                                                                      snapdisp, 
#                                                                      filter_for_frame, 
#                                                                      feature_of_interest)
# 
# print('Expanding regions around features of interest')
# 
# expand_feature_region <- ExpandFeatureRegion(gene_poscodon_codon_i200, 
#                                              gene_names, dataset, hd_file, 
#                                              min_read_length, colsum_out, 
#                                              gff_df, feature_of_interest,
#                                              expand_width, 
#                                              remove_overhang = TRUE,
#                                              snapdisp, filter_for_frame) 
# print('Normalising data')
# 
# normalized_expand_list <- ExpandedRegionNormalisation(gene_poscodon_codon_i200, 
#                                                       gene_names, dataset, hd_file, 
#                                                       min_read_length, colsum_out, 
#                                                       gff_df, feature_of_interest,
#                                                       expand_width, 
#                                                       remove_overhang = TRUE,
#                                                       snapdisp, filter_for_frame)
# print('Calculating average')
# 
# overlayed_tibbles <- OverlayedTibble(gene_poscodon_codon_i200, gene_names, dataset,
#                                      hd_file, min_read_length, colsum_out, gff_df,
#                                      feature_of_interest, expand_width,
#                                      remove_overhang = TRUE, snapdisp, filter_for_frame)
# 
# print('Creating plot')
# 
# plot <- GeneratePlot(gene_poscodon_codon_i200, gene_names, dataset, hd_file,
#                      min_read_length, colsum_out, gff_df, feature_of_interest,
#                      expand_width, remove_overhang = TRUE, snapdisp,
#                      filter_for_frame, size = 12)
# 
# print('Done')
# 








if(length(feature_of_interest) == 1){
  
  # Run ExpandFeatureRegionAllGenes to get a list of occurrances of the codon of interest 
  
  print(paste0('Finding occurances of ', feature_of_interest))
  
  output_feature_info <- suppressMessages(ExpandFeatureRegion(gene_poscodon_codon_i200, 
                                                              gene_names, dataset, hd_file, 
                                                              min_read_length, colsum_out, 
                                                              gff_df, feature_of_interest,
                                                              expand_width, 
                                                              remove_overhang = TRUE,
                                                              snapdisp, filter_for_frame))
  
  # Check for the presence of the feature of interest. Output_feature_info being empty will cause problems with normalization
  
  if(length(output_feature_info) == 0){
    
    print('No occurrances of the feature of interest')
    
    if(expand_width > 1){
      
      print('Try script with an expand_width of 1L to check for occurances near to start or stop codon')
    }
    
    print('Done')
    stop()
  }
  
  # Run ExpandedRegionNormalisation to calculate the relative number of reads mapping to each position arounf the feature of interest
  
  print('Normalising read counts')
  
  normalized_expand_list <- ExpandedRegionNormalisation(gene_poscodon_codon_i200, 
                                                        gene_names, dataset, hd_file, 
                                                        min_read_length, colsum_out, 
                                                        gff_df, feature_of_interest,
                                                        expand_width, 
                                                        remove_overhang = TRUE,
                                                        snapdisp, filter_for_frame)
  
  # Run OverlayedTable to create an average of reads at positions at and around the feature of interest 
  
  print('Overlaying tibbles for feature of interest and calculating the average relative reads at each position')
  
  overlayed_tibbles <- OverlayedTibble(gene_poscodon_codon_i200, gene_names, dataset,
                                       hd_file, min_read_length, colsum_out, gff_df,
                                       feature_of_interest, expand_width,
                                       remove_overhang = TRUE, snapdisp, filter_for_frame)
  
  # Create a graph using ggplot
  
  print('Creating graph')
  
  overlayed_plot <- GeneratePlot(gene_poscodon_codon_i200, gene_names, dataset, hd_file,
                                         min_read_length, colsum_out, gff_df, feature_of_interest,
                                         expand_width, remove_overhang = TRUE, snapdisp,
                                         filter_for_frame, size = 12)
  
  # save plot as PDF
  
  print('Save plot as PDF')
  
  save_plot_pdf <- function(overlayed_plot, output_dir){
    overlayed_plot %>%
      ggsave(
        filename = file.path(output_dir, paste0("Meta_feature_plot", feature_of_interest,".pdf")),
        width = 6, height = 5
      )
  }
  
  save_plot_pdf(overlayed_plot, output_dir)
  
  print('Done')
  
  
} else{
  
  # FindAllFeatures is a function that contains ExpandFeatureRegionAllGenes, ExpandedRegionNormalization and OverlayedTable, which are defined above.
  
  FindAllFeatures <- function(gene_poscodon_codon_i200 = gene_poscodon_codon_i200, 
                              gene_names = gene_names, dataset, hd_file, 
                              min_read_length, colsum_out, 
                              gff_df, .x , 
                              expand_width, remove_overhang, filter_for_frame, snapdisp){
    
    # set .x to feature being studies, allowing purrr::map to iterate over different features of interest without having problems due to .x being a changing vector
    
    feature_being_studied <- .x
    
    
    # Run ExpandFeatureRegionAllGenes to get a list of occurrances of the codon of interest 
    
    print(paste0('Finding occurances of ', feature_being_studied))
    output_feature_info <- suppressMessages(ExpandFeatureRegion(gene_poscodon_codon_i200, 
                                                                gene_names, dataset, hd_file, 
                                                                min_read_length, colsum_out, 
                                                                gff_df, feature_of_interest,
                                                                expand_width, 
                                                                remove_overhang = TRUE,
                                                                snapdisp, filter_for_frame))
    
    
    # Check for the presence of the feature of interest. Output_feature_info being empty will cause problems with normalization
    
    if(length(output_feature_info) == 0){
      print('No occurrances of the feature of interest')
      
      if(expand_width > 1){
        
        print('Try script with an expand_width of 1L to check for occurances near to start or stop codon')
      }
      
      print('Done')
      stop()
    }
    
    
    # Run ExpandedRegionNormalization to calculate the relative number of reads mapping to each position arounf the feature of interest
    
    normalized_expand_list <- ExpandedRegionNormalisation(gene_poscodon_codon_i200, 
                                                          gene_names, dataset, hd_file, 
                                                          min_read_length, colsum_out, 
                                                          gff_df, feature_of_interest,
                                                          expand_width, 
                                                          remove_overhang = TRUE,
                                                          snapdisp, filter_for_frame)
    
    
    # Run OverlayedTable to create an average of reads at positions at and around the feature of interest 
    
    overlayed_tibbles <- OverlayedTibble(gene_poscodon_codon_i200, gene_names, dataset,
                                         hd_file, min_read_length, colsum_out, gff_df,
                                         feature_of_interest, expand_width,
                                         remove_overhang = TRUE, snapdisp, filter_for_frame)
    
    
    
    # Create a new tibble listing the feature being studeied, and the RelCount at position 0, ie RelCount at the feature of interest 
    
    feature_rel_use<- tibble(Feature = feature_being_studied, RelCount = filter(overlayed_tibbles, overlayed_tibbles$Rel_Pos == 0)$RelCount)
    
  }
  
  # Use purrr::map to extract the RelCounts at position 0 of all desired features of interest     
  
  feature_rel_use <- purrr::map_df(.x = feature_of_interest, .f = FindAllFeatures, yeast_codon_pos_i200 = yeast_codon_pos_i200, 
                                   gene_names = gene_names, dataset = dataset, hd_file = hd_file, 
                                   min_read_length = min_read_length, colsum_out = colsum_out, 
                                   gff_df = gff_df, expand_width = expand_width, remove_overhang = remove_overhang, filter_for_frame = filter_for_frame, snapdisp = snapdisp)
  
  # Rearrange feature_rel_use to be in descending order, so features with the highest relative use are listed at the top
  
  feature_rel_use <- arrange(feature_rel_use, desc(RelCount))
  
  # Save feature_rel_use as a tsv file 
  
  write.table(feature_rel_use, file = "Feature_Relativ_use.tsv", sep = "\t", row.names = F, quote = F)
  
  # Users can then look at feature_rel_use and see which features they want to investigate further, and can use as a single feature_of_interest input to produce a graph 
  
} 
















