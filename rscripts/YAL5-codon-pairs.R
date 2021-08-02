
# YAL5_h5 is at location $HOME/riboviz/riboviz/Mok-simYAL5/output/A/A.h5

# Given an h5 file, GFF file and .tsv file, this script creates a metafeature plot 
# for the codon pair of interest. 

print('Starting process')

# source packages and functions from rscripts 
suppressMessages(source(here::here("rscripts", "read_count_functions.R")))
suppressMessages(source(here::here("rscripts", "stats_figs_block_functions.R")))


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
parser$add_argument('--feature', help='Feature of interest, e.g. codon pair')
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
feature_of_interest <- args$feature
output_dir <- args$output
expand_width <- args$expand_width
startpos <- args$startpos
startlen <- args$startlen
filter_for_frame <- args$frame
min_read_length <- args$minreadlen
colsum_out <- args$colsum_out


hd_file <- here::here("Mok-simYAL5", "output", "A", "A.h5")
dataset <- "Mok-simYAL5"
feature_of_interest <- 'TCC AAG'
expand_width = 5L
startpos <-1
startlen <- 10
filter_for_frame <- TRUE
min_read_length <- 10
yeast_codon_table <- here::here("data", "yeast_codon_table.tsv")
gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3")
colsum_out <- TRUE


# hd_file <- here::here("tinysim", "A.h5")
# dataset <- "tinysim"
# gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "tiny_2genes_20utrs.gff3")
# gene <- "MAT"
# gene_names <- unique(gff_df$Name)
# 
# mike <- GetGeneDatamatrix(gene = "MIKE", dataset = "tinysim", hd_file = "tinysim/A.h5")


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
  CodonPos1 = yeast_codon_pos_i200$PosCodon, 
  CodonPos2 = dplyr::lead(yeast_codon_pos_i200$PosCodon),
  CodonPair = paste(yeast_codon_pos_i200$Codon, dplyr::lead(yeast_codon_pos_i200$Codon))
)
  
    # > str(gene_poscodon_codon_i200)
    # tibble [2,826,757 x 4] (S3: tbl_df/tbl/data.frame)
    # $ Gene      : chr [1:2826757] "YAL068C" "YAL068C" "YAL068C" "YAL068C" ...
    # $ CodonPos_1: num [1:2826757] 1 2 3 4 5 6 7 8 9 10 ...
    # $ CodonPos_2: num [1:2826757] 2 3 4 5 6 7 8 9 10 11 ...
    # $ CodonPair : chr [1:2826757] "ATG GTC" "GTC AAA" "AAA TTA" "TTA ACT" ...


#####
# gene <- "YAL003W"
# asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))


print('Fetching read counts')


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



#' FilterForFrameFunction(): fetches A-site assigned counts, filtered for a reading frame
#' 
#' Counts are fetched with GetGeneDatamatrix and A-site assignment is carried out with CalcAsiteFixed (in nt positions)
#' the reading frame is then filtered for and the counts are returned (in codon positions)
#' 
#' @param gene from gene_names to get read lengths for
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' 
#' @return a list of numeric values (counts) for a single reading frame
#' 
#' @example 
#' 
#' FilterForFrameFunction(gene, dataset, hd_file, min_read_length)
#' 
#' @export
FilterForFrameFunction <- function(gene, dataset, hd_file, min_read_length){
  
  asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))
  
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
# TEST: FilterForFrameFunction(): returns a list of numeric values = TRUE
# TEST: FilterForFrameFunction(): length(filtered_counts) = length(CDS) for gene of interest = TRUE
# gives:
# > str(frames)
# num [1:207] 811 429 488 102 994 146 173 762 13 176 ...




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
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if keeping and grouping all reading frames for each codon
#' 
#' @return a tibble which contains the columns "Gene", "PosCodon" and "Count" for a list of genes
#' 
#' @example 
#' 
#' gff_df <- readGFFAsDf(gff)
#' 
#' GetAllCodonPosCounts(gene_names, dataset, hd_file, min_read_length, snapdisp = 0L)
#' 
#' @export 
GetAllCodonPosCounts <- function(gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame){
  
  gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
  
  GetAllCodonPosCounts1Gene <- function(gene, dataset, hd_file, min_read_length, asite_displacement_length, snapdisp, filter_for_frame){
    
    subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
    
    left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start))
    
    right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end))
    
    asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", 
                                                                                "yeast_standard_asite_disp_length.txt"))
    
    reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
    
    reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length, asite_displacement_length)
    
    
    if(filter_for_frame == FALSE){
      
      codon_counts_1_gene <- GetGeneCodonPosReads1dsnap(gene, dataset, hd_file, 
                                                        left, right, min_read_length, 
                                                        asite_displacement_length, snapdisp = 0L)
      
        # > str(codon_counts_1_gene)
        # num [1:207] 4249 825 1017 1176 1116 ...
      
    } else {
      
      codon_counts_1_gene <- FilterForFrameFunction(gene, dataset, hd_file, min_read_length)
      
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
                                           snapdisp,
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
# for filter_for_frame = TRUE
# > str(total_codon_pos_counts)
# Classes 'tbl_df', 'tbl' and 'data.frame':   2749 observations of 3 variables:
#   $ Gene    : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ PosCodon: int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count   : num  811 429 488 102 994 146 173 762 13 176 ...

total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, dataset, hd_file, min_read_length, snapdisp = 0L, filter_for_frame)



#######



print('Creating information tibble')
 



#' AddCodonNamesToCodonPosCounts(): takes codon names from .tsv file, joined to the codon counts table 
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
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, gene_names, dataset = Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE, gff_df)
#' 
#' @export
AddCodonNamesToCodonPosCounts <- function(gene_poscodon_codon_i200, gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame){
  
  total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, dataset, hd_file, min_read_length, snapdisp = 0L, filter_for_frame = filter_for_frame)
  
  transcript_tibbles <- total_codon_pos_counts %>% left_join(gene_poscodon_codon_i200, by = c("PosCodon" = "CodonPos1", "Gene" = "Gene"), keep = FALSE, copy = TRUE)
  
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

transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame)




#####




print('Filtering for feature of interest')


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
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, gene, dataset = Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, snapdisp = 0L, filter_for_frame = TRUE, feature_of_interest = 'TCC AAG')
#' 
#' @export
FilterForFeatureOfInterestPositions <- function(gene_poscodon_codon_i200, gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame, feature_of_interest ){
  
  transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, 
                                                                      gene_names, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      snapdisp, 
                                                                      filter_for_frame)
  
  interesting_feature_table <- dplyr::filter(transcript_gene_pos_poscodon_frame, CodonPair == feature_of_interest)
  
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
  

interesting_feature_positions <- FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, gene_names, dataset, hd_file, min_read_length, snapdisp, filter_for_frame, feature_of_interest)
  


print('Expanding regions around features of interest')



### Expand frame around feature of interest ###



#' ExpandFeatureRegionAllGenes(): Generates a window around each feature of interest across all genes in gene_names
#' 
#' @return a list of tibbles, a tibble for each feature of interest with an expanded window around the feature
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
#' 
#' @example 
#' 
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' ExpandFeatureRegionAllGenes(gene_poscodon_codon_i200, gene_names, dataset, hd_file, min_read_length, colsum_out,
#'                             gff_df, feature_of_interest, expand_width = 5L,
#'                             remove_overhang = TRUE, snapdisp, filter_for_frame)
#' 
#' @export                                                           
ExpandFeatureRegionAllGenes <- function(gene_poscodon_codon_i200, 
                                        gene_names, dataset, hd_file, 
                                        min_read_length, colsum_out, 
                                        gff_df, feature_of_interest,
                                        expand_width = 5L, 
                                        remove_overhang = TRUE,
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
                                         min_read_length, colsum_out, 
                                         gff_df, feature_of_interest, transcript_gene_pos_poscodon_frame,
                                         filter_for_frame){ 
    print(gene)
    
    TranscriptForOneGene <- function(gene_poscodon_codon_i200, 
                                     gene, gene_names, dataset, hd_file, 
                                     min_read_length, colsum_out, 
                                     gff_df, feature_of_interest){
      
      interesting_feature_tibble <- dplyr::filter(transcript_gene_pos_poscodon_frame, CodonPair == feature_of_interest)
      
      transcript_for_one_gene <- dplyr::filter(interesting_feature_tibble, Gene == gene)
      
      return(transcript_for_one_gene)
      
    }
    
    transcript_for_one_gene <- TranscriptForOneGene(gene_poscodon_codon_i200, 
                                                    gene, gene_names, dataset, hd_file, 
                                                    min_read_length, colsum_out, 
                                                    gff_df, feature_of_interest)
    
    
    #if (remove_overhang) {
    # return an empty tibble if the desired region hangs over the edge of the coding region
    
    ExpandRegions <- function(transcript_for_one_gene, transcript_gene_pos_poscodon_frame, gene, dataset, hd_file, expand_width, remove_overhang = TRUE){
      
      interesting_features <- transcript_for_one_gene
      
      transcript_gene_pos_poscodon_gene_interest <- dplyr::filter(transcript_gene_pos_poscodon_frame, Gene == gene)
      
      gene_length <- GetGeneLength(gene, dataset, hd_file)
      
      if (interesting_features <= expand_width  |interesting_features + expand_width > gene_length/3) {
        return()
      } else {
        expand_feature_region <- tibble(
          dplyr::slice(transcript_gene_pos_poscodon_gene_interest, (interesting_features - expand_width):(interesting_features + expand_width), each = FALSE),
          Rel_Pos =  seq(- expand_width, expand_width)
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
    expand_feature_region <- purrr::map(.x = transcript_for_one_gene$CodonPos1, .f = ExpandRegions, 
                                      transcript_gene_pos_poscodon_frame,
                                      gene,
                                      dataset,
                                      hd_file,
                                      expand_width,
                                      remove_overhang = TRUE)
    
    
    return(expand_feature_region)
  }
  
  expand_feature_region <- purrr::map(.x = gene_names, .f = AllGeneInterestingFeatures,
                                    gene_names = gene_names,
                                    gene_poscodon_codon_i200 = gene_poscodon_codon_i200, 
                                    dataset,
                                    hd_file, 
                                    min_read_length,
                                    colsum_out, 
                                    gff_df,
                                    feature_of_interest,
                                    transcript_gene_pos_poscodon_frame)
  
  # produces a list for each gene, containing a list for each occurrence of the feature of interest
  # Unlist to produce one list, containing each occurrence of the feature of interest
  
  expand_feature_region <- unlist(expand_feature_region, recursive = F)
  
  # remove NULLS, which represent features of interest occuring within one expand width of the UTRs
  expand_feature_region <- expand_feature_region[!sapply(expand_feature_region, is.null)]
  
  return(expand_feature_region)
}
#' #TEST: ExpandFeatureRegionAllGenes(): output is a list of tidy format data frames (tibbles) = TRUE. type(expand_feature_region) = "list"
#' #TEST: ExpandFeatureRegionAllGenes(): number of tibbles in list matches the number of occurrences in "feature_of_interest" list = TRUE
#' #TEST: ExpandFeatureRegionAllGenes(): each tibble contains 6 columns = TRUE
#' #TEST: ExpandFeatureRegionAllGenes(): the column names are %in% c("Gene", "Pos_Codon1", "Pos_Codon2", "Count", "CodonPair", "Rel_Pos")
#' #TEST: ExpandFeatureRegionAllGenes(): number of observations in each output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#' #TEST: ExpandFeatureRegionAllGenes(): the position from "interesting_feature_positions" has "Rel_Pos" value 0 = TRUE
#' #TEST: ExpandFeatureRegionAllGenes(): the column "Rel_Pos" goes from -"expand_width to +"expand_width"
#gives:
# > str(expand_feature_region)
# List of 8
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ CodonPos2: num  3 4 5 6 7 8 9 10 11 12 ...
#   $ Count    : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ CodonPair: chr  "GCA TCC" "TCC ACC" "ACC GAT" "GAT TTC" ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  52 53 54 55 56 57 58 59 60 61 ...
#   $ CodonPos2: num  53 54 55 56 57 58 59 60 61 62 ...
#   $ Count    : num  42 53 648 293 121 92 519 79 765 196 ...
#   $ CodonPair: chr  "TTC AAC" "AAC CAC" "CAC ATC" "ATC GCT" ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
#   $ CodonPos1: num  378 379 380 381 382 383 384 385 386 387 ...
#   $ CodonPos2: num  379 380 381 382 383 384 385 386 387 388 ...
#   $ Count    : num  28 74 56 201 4 166 46 15 34 206 ...
#   $ CodonPair: chr  "ACT GGT" "GGT GAC" "GAC GAA" "GAA TCT" ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
#   $ CodonPos1: num  503 504 505 506 507 508 509 510 511 512 ...
#   $ CodonPos2: num  504 505 506 507 508 509 510 511 512 513 ...
#   $ Count    : num  266 468 224 550 42 19 133 31 115 31 ...
#   $ CodonPair: chr  "GAC AAG" "AAG GGT" "GGT AGA" "AGA TTG" ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL005C" "YAL005C" "YAL005C" "YAL005C" ...
#   $ CodonPos1: num  530 531 532 533 534 535 536 537 538 539 ...
#   $ CodonPos2: num  531 532 533 534 535 536 537 538 539 540 ...
#   $ Count    : num  21 32 47 3 7 85 123 40 124 39 ...
#   $ CodonPair: chr  "TCT CAA" "CAA AGA" "AGA ATT" "ATT GCT" ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL012W" "YAL012W" "YAL012W" "YAL012W" ...
#   $ CodonPos1: num  316 317 318 319 320 321 322 323 324 325 ...
#   $ CodonPos2: num  317 318 319 320 321 322 323 324 325 326 ...
#   $ Count    : num  179 47 260 63 24 29 176 20 33 99 ...
#   $ CodonPair: chr  "GGT GCT" "GCT GAA" "GAA GCT" "GCT GCT" ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL035W" "YAL035W" "YAL035W" "YAL035W" ...
#   $ CodonPos1: num  85 86 87 88 89 90 91 92 93 94 ...
#   $ CodonPos2: num  86 87 88 89 90 91 92 93 94 95 ...
#   $ Count    : num  20 57 81 19 39 17 133 144 39 49 ...
#   $ CodonPair: chr  "AAG CCT" "CCT ATA" "ATA CTA" "CTA AAG" ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene     : chr  "YAL038W" "YAL038W" "YAL038W" "YAL038W" ...
#   $ CodonPos1: num  407 408 409 410 411 412 413 414 415 416 ...
#   $ CodonPos2: num  408 409 410 411 412 413 414 415 416 417 ...
#   $ Count    : num  185 446 307 37 0 286 271 171 124 161 ...
#   $ CodonPair: chr  "ACC CCA" "CCA AGA" "AGA TTG" "TTG GTT" ...
#   $ Rel_Pos  : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...

expand_feature_region <- ExpandFeatureRegionAllGenes(gene_poscodon_codon_i200, 
                                              gene_names, dataset, hd_file, 
                                              min_read_length, colsum_out, 
                                              gff_df, feature_of_interest,
                                              expand_width = 5L, 
                                              remove_overhang = TRUE,
                                              snapdisp, filter_for_frame) 


### Normalization ###


print('Normalising data')

#' ExpandedRegionNormalization(): carries out normalization within each expanded frame so that they are comparable
#' 
#' Normalizes the ExpandRegion list generating a RelCount column with the normalization values 
#' As this function is not looped it will only generate one normalized tibble for one occurence 
#' 
#' @param .x which is the list of tidy format data frames (tibbles) generated by the function ExpandFeatureRegion
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' 
#' @example 
#' 
#' ExpandedRegionNormalization(expand_feature_region, expand_width = 5L)
#' 
#' @export 
ExpandedRegionNormalization <- function(.x, expand_width){
  
  # dplyr::mutate(.x, RelCount = PerCodonCounts / sum(PerCodonCounts) * (2 * expand_width + 1))
  dplyr::mutate(.x, RelCount = Count / sum(Count) * (2 * expand_width + 1))
  
}
#TEST: ExpandedRegionNormalization(): creates a tidy format data frame (tibble) = TRUE
#TEST: ExpandedRegionNormalization(): the tibble contains 6 columns = TRUE
#TEST: ExpandedRegionNormalization(): the column names are %in% c("Gene", "Pos_Codon1", "Pos_Codon2", "Rel_Count", "Rel_Pos", "RelCount")
#TEST: ExpandedRegionNormalization(): number of observations in the output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandedRegionNormalization(): the column "Rel_Pos" goes from -"expand_width to +"expand_width" 
#gives:
# > str(normalized_expanded_feature_region)
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos_Codon1: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Pos_Codon2: num  3 4 5 6 7 8 9 10 11 12 ...
#   $ Rel_Count : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ Rel_Pos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
#   $ RelCount  : num  1.347 1.532 0.32 3.12 0.458 ...

# normalized_expanded_feature_region <- ExpandedRegionNormalization(expand_feature_region[[1]], expand_width = 5L)



# Normalization carried out for all the tibbles within ExpandList 
normalized_expand_list <- purrr::map(
  .x = expand_feature_region,
  .f = ExpandedRegionNormalization,
  expand_width = 5L
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
# from NormalizedExpandList. Need to join by Rel_Pos, in RelCount need the mean for 
# each Rel_Pos (sum row(x) / number of row(x))

#' OverlayedTable(): overlays tidy format data frames (tibbles) to create a single overlayed tibble.
#' 
#' FIXME: Takes normalized_expand_list as its input.
#' 
#' @param normalized_expand_list the output from the looped function ExpandedRegionNormalization()
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' 
#' @example 
#' 
#' normalized_expand_list <- purrr::map(.x = expand_feature_region,
#'                                      .f = ExpandedRegionNormalization,
#'                                      expand_width)
#' 
#' OverlayedTable(normalized_expand_list, expand_width)
#' 
#' @export
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
#TEST: OverlayedTable(): creates a tidy format data frame (tibble) = TRUE
#TEST: OverlayedTable(): the tibble contains 2 columns = TRUE
#TEST: OverlayedTable(): the column names are %in% c("Rel_Pos", "RelCount")
#TEST: OverlayedTable(): number of observations in the output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: OverlayedTable(): the column "Rel_Pos" goes from -"expand_width to +"expand_width" 
#TEST: OverlayedTable(): RelCount is a numeric 
#gives:
# > str(overlayed_tibbles)
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 2 variables
# $ Rel_Pos : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
# $ RelCount: num [1:11] 1.029 1.07 0.987 1.45 1.151 ...

overlayed_tibbles <- OverlayedTable(normalized_expand_list, expand_width) 


print('Creating plot')


overlayed_plot <- ggplot(overlayed_tibbles, mapping = aes(x = Rel_Pos, y = RelCount)) + 
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        title = element_text(size = 12, face = 'bold'))+
  labs(title = paste0('Meta-feature plot of codon pair ', feature_of_interest),
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




