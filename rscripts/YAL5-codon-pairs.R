
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
# parser$add_argument('--frame', help='frame to be studied', default = 0)
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
# filtering_frame <- args$frame
min_read_length <- args$minreadlen
colsum_out <- args$colsum_out


hd_file <- here::here("Mok-simYAL5", "output", "A", "A.h5")
dataset <- "Mok-simYAL5"
feature <- 'TCC AAG'
expand_width = 5L
startpos <-1
startlen <- 10
# filtering_frame <- 0
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
gene <- "YAL003W"
asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))

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



FilterForFrameFunction <- function(reads_asitepos, left, right){
  
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
  
frames <- FilterForFrameFunction(reads_asitepos, left, right)
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
GetAllCodonPosCounts <- function(gene_names, dataset, hd_file, min_read_length, snapdisp){
  
  gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
  
  GetAllCodonPosCounts1Gene <- function(gene, dataset, hd_file, min_read_length, asite_displacement_length, snapdisp){
    
    subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
    
    left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start))
    
    right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end))
    
    asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", 
                                                                                "yeast_standard_asite_disp_length.txt"))
    
    
    if(SnapCodon == TRUE){
      
      codon_counts_1_gene <- GetGeneCodonPosReads1dsnap(gene, dataset, hd_file, 
                                                        left, right, min_read_length, 
                                                        asite_displacement_length, snapdisp = 0L)
      
        # > str(codon_counts_1_gene)
        # num [1:207] 4249 825 1017 1176 1116 ...
      
    } else {
      
      codon_counts_1_gene <- FilterForFrameFunction(reads_asitepos, left, right)
      
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
  )
  
  return (total_codon_pos_counts)
}
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

total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, dataset, hd_file, min_read_length, snapdisp = 0L)


#######





print('Creating information tibble')
 


         
        #' # The function CreateTranscriptInfoTibbleAllGenes applies all of the functions up until and including filtering for frame of interest to all of the genes 
        #' # it returns one tibble, containing all of the genes and information
        #' # to run for a single gene or to run each function separately, filter for a gene of interest and run
        #' 
        #' # Filter down the gene_poscodon_codon_i200 file to gene of interest, e.g. YAL003W
        #' YAL003W_pos <- dplyr::filter(gene_poscodon_codon_i200, Gene=="YAL003W")
        #' gene <- 'YAL003W'
        #' 
        #' 
        #' CreateTranscriptInfoTibbleAllGenes <- function(gene, dataset, hd_file, gff_df, startpos, startlen, colsum_out){
        #' 
        #'           
        #' ### Functions for A-site assignment of reads extracted from .h5 file ###
        #' 
        #'           
        #' #' TidyAsiteCountsByPosition(): Align A-site assigned counts to the nucleotide positions for a single gene
        #' #' 
        #' #' The function ReadAsiteDisplacementLengthFromFile() from rscripts: is used within this function.
        #' #' The function GetGeneDatamatrix() from rscripts: is used within this function. 
        #' #' The function CalcAsiteFixed() from rscripts: read_count_functions is used within this function.
        #' #' 
        #' #' @param dataset name of dataset stored in .h5 file.
        #' #' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
        #' #' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
        #' #' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
        #' #' 
        #' #' @return a tibble with the columns "Pos" and "Count"
        #' #' 
        #' #' @examples
        #' #' 
        #' #' TidyAsiteCountsByPosition(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE)
        #' #' 
        #' #' @export
        #' TidyAsiteCountsByPosition <- function(gene, dataset, hd_file, min_read_length, colsum_out){
        #'   
        #'   asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))
        #'   #This step extracts the asite reads from the reads_pos_length file  
        #'   
        #'     # > str(asite_displacement_length)
        #'     # spec_tbl_df [3 x 2] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
        #'     # $ read_length       : num [1:3] 28 29 30
        #'     # $ asite_displacement: num [1:3] 15 15 15
        #'             
        #'   reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
        #'   
        #'     # > str(reads_pos_length) (for YAL003W)
        #'     # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
        #'               
        #'   asite_counts_by_position <- CalcAsiteFixed(reads_pos_length,
        #'                                              min_read_length = min_read_length,
        #'                                              asite_displacement_length = asite_displacement_length,
        #'                                              colsum_out = colsum_out
        #'                                              )
        #'     # > str(asite_counts_by_position)(for YAL003W)
        #'     # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
        #'             
        #'   tidy_asite_counts <- tibble(Pos = 1:length(asite_counts_by_position), Count = asite_counts_by_position)
        #'             
        #'   return(tidy_asite_counts)
        #'             
        #' }
        #' #TEST: TidyAsiteCountsByPosition(): returns a a tidy format data frame (tibble). 
        #' #TEST: TidyAsiteCountsByPosition(): the tibble has 2 columns.
        #' #TEST: TidyAsiteCountsByPosition(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df for the gene.
        #' #TEST: TidyAsiteCountsByPosition(): the column names are %in% c("Pos", "Count")
        #' # gives:          
        #' # > str(tidy_asite_count_output)
        #' # Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 2 variables:
        #' #   $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
        #' #   $ Count: num  0 0 0 0 0 0 0 0 0 0 ...
        #'           
        #' tidy_asite_counts <- TidyAsiteCountsByPosition(gene, dataset, hd_file, min_read_length, colsum_out)
        #'           
        #'       # The end result here is that the A-site assigned counts are aligned to the gene of interest
        #'       # (in nucleotides, including UTRs and UTRs)
        #'           
        #' 
        #' #####
        #' 
        #' 
        #'  
        #' ### Functions for assembling tibbles consisting of counts, positions and codons ###
        #'   
        #' 
        #' 
        #' #' TranscriptPosToCodonPos(): Assigns the nucleotide positions of UTRs and CDS to codon positions
        #' #' 
        #' #' @param gene from gene_names to pull information from the gff_df file
        #' #' @param gff_df from which to extract the UTRs and CDS widths. 
        #' #' 
        #' #' @return Tidy data frame (tibble) containing the columns: "Gene", "Pos" 
        #' #' (position of nucleotides), "Pos_Codon" (leading codon), "Pos_Codon2" (lagging codon) 
        #' #' and "Frame" (reading frame). 
        #' #' 
        #' #' @example 
        #' #' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
        #' #' 
        #' #' TranscriptPosToCodonPos(gene = "YAL003W", gff_df)
        #' #' @export
        #' TranscriptPosToCodonPos <- function(gene, gff_df){
        #'     
        #'     subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
        #'     # the gff file is filtered down to the gene of interest
        #'     
        #'     UTR5_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR5") %>% select(width)
        #'     # the width of the UTR5 is defined from the gff file, in this case 250 
        #'     
        #'     CDS_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(width)
        #'     # the width of the CDS is defined from the gff file, in this case 621
        #'     
        #'     UTR3_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(width) 
        #'     # the width of the UTR3 is defined from the gff file, in this case 250 
        #'     
        #'     transcript_length <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(end)
        #'     # defines the length of the transcript including both UTRs and the CDS (in nucleotides), in this case 1121
        #'     
        #'     CDS_start <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start)
        #'     # define the start position of the CDS from the gff file, in this case 251
        #'     
        #'     CDS_end <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end)
        #'     # define the end position of the CDS from the gff file, in this case 871
        #'     
        #'     CDS_codon_positions <- rep(1:(CDS_width$width/3), each = 3)
        #'     # In order to add the CDS codon positions the length of the CDS in nucleotides is 
        #'     # divided by 3 (3 nucleotides = 1 codon) and repeated 3 times
        #'       
        #'       # > str(CDS_codon_positions)
        #'       # int [1:621] 1 1 1 2 2 2 3 3 3 4 ...
        #'       
        #'     NA_UTR5_width <- rep(NA, times = UTR5_width$width)
        #'     # NAs are repeated for the width of UTR5 
        #'     
        #'     NA_UTR3_width <- rep(NA, times = UTR3_width$width)
        #'     # NAs are repeated for the width of UTR3
        #'     
        #'     transcript_gene_pos_poscodon_frame <- tibble(
        #'       Gene = gene,
        #'       Pos = 1:transcript_length$end,
        #'       Pos_Codon1 = c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)),
        #'       Pos_Codon2 = dplyr::lead(c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)), n = 3),
        #'       Frame = seq(from = 2L, to = (transcript_length$end + 1L)) %% 3 # works as UTRS are 250, might not work with UTRS of others values
        #'       # add the count column, likely something we want here 
        #'     )
        #'     
        #'     # For codon pairs name of Pos_Codon became Pos_Codon1 and Pos_Codon2 was added 
        #'     # where the line for Pos_Codon1 was copied and dplyr::lead((), n = 3) was added 
        #'     
        #'     return(transcript_gene_pos_poscodon_frame)
        #' }
        #' #TEST: TranscriptPosToCodonPos(): creates a tibble = TRUE
        #' #TEST: TranscriptPosToCodonPos(): the tibble contains 5 columns = TRUE
        #' #TEST: TranscriptPosToCodonPos(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
        #' #TEST: TranscriptPosToCodonPos(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame") 
        #' #TEST: TranscriptPosToCodonPos(): The first and last 250 rows of the columns "Pos_codon1" and "Pos_Codon2" for each gene contain NA. 
        #' # gives:
        #' # > str(transcript_pos_to_codon_pos_output)
        #' # Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 5 variables:
        #' #   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
        #' #   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
        #' #   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
        #' #   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
        #' #   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
        #' 
        #' 
        #' transcript_gene_pos_poscodon_frame <- TranscriptPosToCodonPos(gene, gff_df)
        #' 
        #'   # The end result is a table with the columns Gene, Pos, Pos_Codon and Frame (reading frame).
        #'   # For the Pos_Codon column the UTR positions have NA, while the CDS has the codon positions 
        #'              
        #' 
        #' #' AddAsiteCountsToTranscriptPosToCodonPos(): merges A-site counts with transcript_pos_to_codon_pos_output
        #' #' 
        #' #' This is a helper function for AddAsiteCountsToTranscriptPosToCodonPosAllGenes 
        #' #' TidyAsiteCountsByPosition() and TranscriptPosToCodonPos() are used in this function. 
        #' #' 
        #' #' @param gene in gene_names to pull information from the gff_df file 
        #' #' @param dataset name of dataset stored in .h5 file.
        #' #' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
        #' #' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
        #' #' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
        #' #' @param gff_df from which to extract the UTRs and CDS widths.
        #' #' 
        #' #' @return a list of tibbles where the A-site assigned counts are aligned to the nucleotide positions of each gene of interest (UTRs and CDS)
        #' #' 
        #' #' @examples 
        #' #' 
        #' #' do TidyAsiteCountsByPosition() and TranscriptPosToCodonPos() need to be defined here?
        #' #' 
        #' #' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
        #' #' 
        #' #' AddAsiteCountsToTranscriptPosToCodonPos(gene = "YAL003W", dataset = Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE, gff_df)
        #' #' 
        #' #' @export
        #' AddAsiteCountsToTranscriptPosToCodonPos <- function(gene, dataset, hd_file, min_read_length, colsum_out, gff_df){
        #'     
        #'     tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, dataset, hd_file, min_read_length, colsum_out = colsum_out)
        #'     
        #'     transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene, gff_df)
        #'     
        #'     transcript_info_tibble <- dplyr::left_join(transcript_pos_to_codon_pos_output, tidy_asite_count_output, by = "Pos")
        #'     
        #'     return(transcript_info_tibble)
        #' }
        #' #TEST: AddAsiteCountsToTranscriptPosToCodonPos(): creates a tidy format data frame (tibble) = TRUE
        #' #TEST: AddAsiteCountsToTranscriptPosToCodonPos(): the tibble contains 6 columns = TRUE
        #' #TEST: AddAsiteCountsToTranscriptPosToCodonPos(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
        #' #TEST: AddAsiteCountsToTranscriptPosToCodonPos(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame", "Count") 
        #' # gives:
        #' # > str(transcript_info_tibble)
        #' # Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 6 variables:
        #' #   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
        #' #   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
        #' #   $ Pos_Codon1: int  NA NA NA NA NA NA NA NA NA NA ...
        #' #   $ Pos_Codon2: int  NA NA NA NA NA NA NA NA NA NA ...
        #' #   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
        #' #   $ Count     : num  0 0 0 0 0 0 0 0 0 0 ...
        #' 
        #'     # End result is a table which contains the asite counts, reading frame, codons and 
        #'     # codon positions for the gene of interest
        #'   
        #' transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
        #' 
        #' 
        #' 
        #' #' AddCodonNamesToTranscriptInfoTibble(): Provides codon names from the .tsv file to the 
        #' #' tidy format data frame (tibble) generated from AddAsiteCountsToTranscriptPosToCodonPos   
        #' #' 
        #' #' AddAsiteCountsToTranscriptPosToCodonPos() is used in this function
        #' #' 
        #' #' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
        #' #' @param gene from gene_names   
        #' #' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
        #' #' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
        #' #' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
        #' #' @param gff_df from which to extract the UTRs and CDS widths.
        #' #' 
        #' #' @example 
        #' #' 
        #' #' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
        #' #' 
        #' #' AddCodonNamesToTranscriptInfoTibble(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df)
        #' #' 
        #' #' @export                                                                     
        #' AddCodonNamesToTranscriptInfoTibble <- function(gene_poscodon_codon_i200, dataset = Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE, gff_df gff_df){
        #'     
        #'     codon_table <- dplyr::filter(gene_poscodon_codon_i200, Gene == gene) %>% 
        #'       dplyr::select("CodonPos_1", "CodonPos_2", "CodonPair")
        #'     # this provides the positions of the codons and the codon names 
        #'     
        #'     transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, dataset = dataset, hd_file = hd_file, min_read_length = min_read_length, colsum_out = colsum_out, gff_df)
        #'     # this provides the tibble with the columns Gene, Pos, Pos_Codon and Frame (reading frame)
        #'     
        #'     transcript_info_tibble <- left_join(transcript_info_tibble, codon_table, by = c("Pos_Codon1" = "CodonPos_1", "Pos_Codon2" = "CodonPos_2"), keep = FALSE)
        #' 
        #'     return(transcript_info_tibble)
        #' }
        #' #TEST: AddCodonNamesToTranscriptInfoTibble(): creates a tidy format data frame (tibble) = TRUE
        #' #TEST: AddCodonNamesToTranscriptInfoTibble(): the tibble contains 7 columns = TRUE
        #' #TEST: AddCodonNamesToTranscriptInfoTibble(): number of observations in the output tibble = width of UTR5+CDS+UTR3 from gff_df, for YAL003W = 1121.
        #' #TEST: AddCodonNamesToTranscriptInfoTibble(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame", "Count", "CodonPair") 
        #' #TEST: AddCodonNamesToTranscriptInfoTibble(): the first and last 250 rows of the columns "Pos_codon1", "Pos_Codon2" and "CodonPair" for each gene contain NA. 
        #' # gives:
        #' # > str(transcript_info_tibble)
        #' # Classes 'tbl_df', 'tbl' and 'data.frame':   1121 observations of 7 variables:
        #' #   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
        #' #   $ Pos       : int  1 2 3 4 5 6 7 8 9 10 ...
        #' #   $ Pos_Codon1: num  NA NA NA NA NA NA NA NA NA NA ...
        #' #   $ Pos_Codon2: num  NA NA NA NA NA NA NA NA NA NA ...
        #' #   $ Frame     : num  2 0 1 2 0 1 2 0 1 2 ...
        #' #   $ Count     : num  0 0 0 0 0 0 0 0 0 0 ...
        #' #   $ CodonPair : chr  NA NA NA NA ...
        #' 
        #' transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df)
        #'  
        #' 
        #'  
        #' }





#' AddCodonNamesToCodonPosCounts(): takes codon names from .tsv file, joined to the codon counts table 
#' 
#' Uses the function GetAllCodonPosCounts().
#' 
#' @return a  tidy format data frame (tibble) which contains the genes in gene_names, codon positions, counts and the codon pair. 
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
#' AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, gene_names, dataset = Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE, gff_df)
#' 
#' @export
AddCodonNamesToCodonPosCounts <- function(gene_poscodon_codon_i200, gene_names, dataset, hd_file, min_read_length, colsum_out, gff_df){
  
  total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, dataset, hd_file, min_read_length, snapdisp = 0L)
  
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
#TEST: AddCodonNamesToCodonPosCounts(): the column names are %in% c("Gene", "Pos_Codon1", "Pos_Codon2", "Count", "CodonPair") 
#TEST: AddCodonNamesToCodonPosCounts(): the unique gene names in column "Gene" match the genes in gene_names (unique(total_codon_pos_counts$Gene) = gene_names) = TRUE.
# gives:
# > str(transcript_gene_pos_poscodon_frame)
# Classes 'tbl_df', 'tbl' and 'data.frame':   2,749 observations of 5 variables:
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  1 2 3 4 5 6 7 8 9 10 ...
#   $ CodonPos2: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Count    : num  4249 825 1017 1176 1116 ...
#   $ CodonPair: chr  "ATG GCA" "GCA TCC" "TCC ACC" "ACC GAT" ...

transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df)

          
          
          #' # single_gene_info_tibble <- CreateTranscriptInfoTibbleAllGenes(gene, dataset, hd_file, gff_df, filtering_frame, startpos, startlen)
          #' 
          #' # CreateTranscriptInfoTibbleAllGenes (reads_pos_length, gene, dataset, hd_file, gff_df)
          #' transcript_info_tibble <- suppressMessages(purrr::map_dfr(.x = gene_names,
          #'                                                           .f = CreateTranscriptInfoTibbleAllGenes,
          #'                                                           dataset,
          #'                                                           hd_file,
          #'                                                           gff_df,
          #'                                                           colsum_out,
          #'                                                           startpos,
          #'                                                           startlen))
          #' 
          #' # this takes all of the functions defined above and applies them to all of the genes in the sample.
          #' 
          #' # For each gene, the first and last 250 nt may have NA in the codon column and Pos_Codon
          #' # This is due to the 250 nt utr buffer, present in many riboviz example datasets
          #' 
          #' # > str(transcript_info_tibble)
          #' # tibble [10,747 x 7] (S3: tbl_df/tbl/data.frame)
          #' # $ Gene      : chr [1:10747] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
          #' # $ Pos       : int [1:10747] 1 2 3 4 5 6 7 8 9 10 ...
          #' # $ Pos_Codon1: num [1:10747] NA NA NA NA NA NA NA NA NA NA ...
          #' # $ Pos_Codon2: num [1:10747] NA NA NA NA NA NA NA NA NA NA ...
          #' # $ Frame     : num [1:10747] 2 0 1 2 0 1 2 0 1 2 ...
          #' # $ Count     : num [1:10747] 0 0 0 0 0 0 0 0 0 0 ...
          #' # $ CodonPair : chr [1:10747] NA NA NA NA ...
          #' 
          #' 
          #' 
          #' #' FilterForFrame(): Filter for reading frame of interest
          #' #' 
          #' #' The default reading frame is 0
          #' #' FIXME: Don't want to have the output from a function as the input for another function
          #' #' 
          #' #' @param transcript_info_tibble from AddCodonNamesToTranscriptInfoTibble 
          #' #' @param filtering_frame which sets the reading frame that you filter for 
          #' #' 
          #' #' @example 
          #' #' 
          #' #' FilterForFrame(transcript_info_tibble, filtering_frame = 0)
          #' #' 
          #' #' @export
          #' FilterForFrame <- function(transcript_info_tibble, filtering_frame){
          #'   
          #'   transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df)
          #'     
          #'     transcript_info_tibble <- dplyr::filter(transcript_info_tibble, 
          #'                                             Frame == filtering_frame)
          #'     
          #'     return(transcript_info_tibble)
          #'     
          #'   }
          #' #TEST: FilterForFrame(): creates a tidy format data frame (tibble) = TRUE
          #' #TEST: FilterForFrame(): the tibble contains 7 columns = TRUE
          #' #TEST: FilterForFrame(): the number of objects in gene_names is equal to the number of unique objects in transcript_info_tibble$Gene
          #' #TEST: FilterForFrame(): The "Frame" column only contains the parameter you filtered for, e.g. 0.
          #' #TEST: AddCodonNamesToTranscriptInfoTibble(): the column names are %in% c("Gene", "Pos", "Pos_Codon1", "Pos_Codon2", "Frame", "Count", "CodonPair") 
          #' #TEST: AddCodonNamesToTranscriptInfoTibble(): the first and last ~83 rows of the columns "Pos_codon1", "Pos_Codon2" and "CodonPair" for each gene contain NA. 
          #' # gives: 
          #' # > str(transcript_info_tibble) 
          #' # Classes 'tbl_df', 'tbl' and 'data.frame':   3,584 observations of 7 variables
          #' #   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
          #' #   $ Pos       : int  2 5 8 11 14 17 20 23 26 29 ...
          #' #   $ Pos_Codon1: num  NA NA NA NA NA NA NA NA NA NA ...
          #' #   $ Pos_Codon2: num  NA NA NA NA NA NA NA NA NA NA ...
          #' #   $ Frame     : num  0 0 0 0 0 0 0 0 0 0 ...
          #' #   $ Count     : num  0 0 0 0 0 0 0 0 0 0 ...
          #' #   $ CodonPair : chr  NA NA NA NA ...
          #' 
          #' transcript_info_tibble <- FilterForFrame(transcript_info_tibble, filtering_frame)





###




print('Filtering for feature of interest')


#' FilterForFeatureOfInterestPositions(): Filters for the feature of interest 
#' 
#' Uses the AddCodonNamesToCodonPosCounts to get counts 
#' 
#' @return a tidy format data frame (tibble) which contains the gene, codon positions and counts for the feature of interest 
#' 
#' @param .tsv file from which to fetch the codon names associated with the CDS co-ordinates for each gene
#' @param gene from gene_names   
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param feature_of_interest character, each incidence of the feature will be extracted from transcript_info_tibble
#' 
#' @example 
#' 
#' gff_df <- readGFFAsDf(here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, gene, dataset = Mok-simYAL5, hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), min_read_length = 10, colsum_out = TRUE, gff_df, feature_of_interest)
#' 
#' @export
FilterForFeatureOfInterestPositions <- function(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df, feature_of_interest ){
  
  transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, 
                                                                      gene, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      colsum_out, 
                                                                      gff_df)
  
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
  
    # # Which refers to the position of the first codon in the codon pair as the second codon
    # # in the pair is the first codon position + 1

interesting_feature_positions <- FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, gene, dataset, hd_file, min_read_length, colsum_out, gff_df, feature_of_interest)
  


### Expand frame around feature of interest ###


#' ExpandFeatureRegion(): Expands the window around the feature position of interest
#' 
#' Helperfunction for ExpandFeatureRegionForList()
#' This function is not looped and will only generate one tidy format data frame (tibble)
#' 
#' @return an expanded window around the feature of interest which allows you to compare the rate of 
#' translation for the feature of interest to its surrounding translational environment 
#' 
#' @param gene in gene_names to pull information from the gff_df file 
#' @param transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes()
#' @param interesting_feature_positions list of numeric values describing position of the feature of interest
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' @param remove_overhang default = TRUE, removes features of interest who's positions do not allow for complete windows to be generated
#' 
#' @example 
#' 
#' ExpandFeatureRegion(.x = interesting_feature_positions, 
#'                       transcript_info_tibble, 
#'                       gene = "YAL003W", 
#'                       dataset = "Mok-simYAL5", 
#'                       hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"), 
#'                       expand_width = 5L, 
#'                       remove_overhang = TRUE)
#'                       
#' @export                       
ExpandFeatureRegion <- function(gene_poscodon_codon_i200, 
                                gene, dataset, hd_file, 
                                min_read_length, colsum_out, 
                                gff_df, feature_of_interest,
                                  expand_width = 5L, 
                                  remove_overhang = TRUE) {

  transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, 
                                                                      gene, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      colsum_out, 
                                                                      gff_df)

  
  TranscriptForOneGene <- function(gene_poscodon_codon_i200, 
                                   gene, dataset, hd_file, 
                                   min_read_length, colsum_out, 
                                   gff_df, feature_of_interest){
    
    interesting_feature_tibble <- FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, 
                                                                         gene, dataset, hd_file, 
                                                                         min_read_length, colsum_out, 
                                                                         gff_df, feature_of_interest)
    
    transcript_for_one_gene <- dplyr::filter(interesting_feature_tibble, Gene == gene)
    
    interesting_feature_positions <- transcript_for_one_gene$CodonPos1
    
    return(interesting_feature_positions)
    
  }
  
  interesting_features <- TranscriptForOneGene(gene_poscodon_codon_i200, 
                                                    gene, dataset, hd_file, 
                                                    min_read_length, colsum_out, 
                                                    gff_df, feature_of_interest)
  
  gene_length <- GetGeneLength(gene, dataset, hd_file)

  #if (remove_overhang) {
  # return an empty tibble if the desired region hangs over the edge of the coding region

  if (interesting_features <= expand_width  |interesting_features + expand_width > gene_length/3) {
    return()
  } else {
    output_feature_info <- tibble(
      dplyr::slice(transcript_gene_pos_poscodon_frame, (interesting_features - expand_width):(interesting_features + expand_width), each = FALSE),
      Rel_Pos =  seq(- expand_width, expand_width)
    )
    
    if(dim(output_feature_info)[1] == (2*expand_width + 1)){
      
    return(output_feature_info)
  }else{
    return()
    }
  }
  # The if statement ensures that feature positions that are less/more than the 
  # expand_width value are discarded 
}
#TEST: ExpandFeatureRegion(): creates a tidy format data frame (tibble) = TRUE
#TEST: ExpandFeatureRegion(): the tibble contains 5 columns = TRUE
#TEST: ExpandFeatureRegion(): the column names are %in% c("Gene", "CodonPos1", "CodonPos2", "Count", "Rel_Pos") 
#TEST: ExpandFeatureRegion(): number of observations in the output tibble = "expand_width" * 2 + 1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandFeatureRegion(): the position from "interesting_feature_positions" has "Rel_Pos" value 0 = TRUE
#TEST: ExpandFeatureRegion(): the column "Rel_Pos" goes from -"expand_width to +"expand_width"
#gives: 
# > str(output_feature_info)
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 6 variables
# $ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
# $ CodonPos1: num [1:11] 2 3 4 5 6 7 8 9 10 11 ...
# $ CodonPos2: num [1:11] 3 4 5 6 7 8 9 10 11 12 ...
# $ Count    : num [1:11] 825 1017 1176 1116 284 ...
# $ CodonPair: chr [1:11] "GCA TCC" "TCC ACC" "ACC GAT" "GAT TTC" ...
# $ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...

output_feature_info <- ExpandFeatureRegion(gene_poscodon_codon_i200, 
                                           gene, dataset, hd_file, 
                                           min_read_length, colsum_out, 
                                           gff_df, feature_of_interest, 
                                           expand_width, remove_overhang)








listed <- purrr::map(.x = gene_names,
                     .f = ExpandFeatureRegion,
                     gene_poscodon_codon_i200, 
                     dataset, hd_file, 
                     min_read_length, colsum_out, 
                     gff_df, feature_of_interest, 
                     expand_width, remove_overhang)




#' ExpandFeatureRegionForList(): Applies the ExpandFeatureRegion function to all of the features of interest. 
#' 
#' This function generates expanded tibbles for each position contained within the feature_of_interest list.
#' Each occurence in feature_of_interest list generates one tibble 
#' 
#' FIXME: Takes feature_of_interest and transcript_info_tibble as inputs 
#' FIXME: THIS STILL CREATES TIBBLES OF CODON PAIRS OF INTEREST IN THE WRONG GENES
#' 
#' @param transcript_info_tibble from CreateTranscriptInfoTibbleAllGenes()
#' @param gene in gene_names to pull information from the gff_df file 
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created from BAM files for dataset samples.
#' @param startpos position of the start codon, default = 1
#' @param startlen, smallest length of reads', default = 10
#' @param feature_of_interest character, each incidence of the codon pair will be extracted from transcript_info_tibble
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param expand_width integer which provides the number of positions on each side of the feature of interest to include in the window
#' 
#' @example 
#' 
#' ExpandFeatureRegionForList(transcript_info_tibble, 
#'                              gene = unique(transcript_info_tibble$Gene), 
#'                              dataset = "Mok-simYAL5", 
#'                              hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5") 
#'                              startpos = 1, 
#'                              startlen = 10, 
#'                              feature_of_interest = "TCC AAG", 
#'                              gff_df, 
#'                              expand_width = 5L)
#' 
#' @export
ExpandFeatureRegionForList <- function(gene_poscodon_codon_i200, 
                                       gene, 
                                       dataset, 
                                       hd_file, 
                                       min_read_length, 
                                       colsum_out, 
                                       gff_df, 
                                         startpos, 
                                         startlen, 
                                         feature_of_interest,
                                         expand_width){
  

  
  
  transcript_gene_pos_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_poscodon_codon_i200, 
                                                                      gene, 
                                                                      dataset, 
                                                                      hd_file, 
                                                                      min_read_length, 
                                                                      colsum_out, 
                                                                      gff_df)
  
  interesting_feature_positions <- FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, 
                                                                       gene, dataset, hd_file, 
                                                                       min_read_length, colsum_out, 
                                                                       gff_df, feature_of_interest)
  
  interesting_features <- interesting_feature_positions$CodonPos1
  
  expand_feature_region <- list()
  
  
  # A loop is used here so each unique gene can be processed separately. there were issues with using map and lapply 
  # on the filter function of FilterForFeatureOfInterestPositions, and filtering using only one gene would produce NULLs
  
  for(gene_name in unique(gene)){
  
    interesting_feature_positions <- FilterForFeatureOfInterestPositions(gene_poscodon_codon_i200, 
                                                                         gene, dataset, hd_file, 
                                                                         min_read_length, colsum_out, 
                                                                         gff_df, feature_of_interest) 
  
    transcript_info_tibble_gene <- dplyr::filter(transcript_gene_pos_poscodon_frame, Gene == gene_name)
    
    tmp_expand_feature_region <- purrr::map(
      .x = interesting_feature_positions,
      .f = ExpandFeatureRegion,
      transcript_gene_pos_poscodon_frame,
      gene = gene_name,
      dataset,
      hd_file,
      gff_df,
      expand_width = expand_width
    )
   
    # remove any tibbles from the list that are null, so are within one expand_width of the UTRs
    tmp_expand_feature_region <- tmp_expand_feature_region[!sapply(tmp_expand_feature_region,is.null)]
    
    
    expand_feature_region <- append(expand_feature_region, tmp_expand_feature_region)
    
  }
  return(expand_feature_region)
}
#TEST: ExpandFeatureRegionForList(): output is a list of tidy format data frames (tibbles) = TRUE. type(expand_feature_region) = "list"
#TEST: ExpandFeatureRegionForList(): number of tibbles in list matches the number of occurrences in "feature_of_interest" list = TRUE
#TEST: ExpandFeatureRegionForList(): each tibble contains 5 columns = TRUE
#TEST: ExpandFeatureRegionForList(): the column names are %in% c("Gene", "Pos_Codon1", "Pos_Codon2", "Rel_Count", "Rel_Pos")
#TEST: ExpandFeatureRegionForList(): number of observations in each output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandFeatureRegionForList(): the position from "interesting_feature_positions" has "Rel_Pos" value 0 = TRUE
#TEST: ExpandFeatureRegionForList(): the column "Rel_Pos" goes from -"expand_width to +"expand_width"
#gives:
# > str(expand_feature_region)
# List of 8
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 5 variables
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos_Codon1: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Pos_Codon2: num  3 4 5 6 7 8 9 10 11 12 ...
#   $ Rel_Count : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ Rel_Pos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes 'tbl_df', 'tbl' and 'data.frame':   11 observations of 5 variables
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos_Codon1: num  52 53 54 55 56 57 58 59 60 61 ...
#   $ Pos_Codon2: num  53 54 55 56 57 58 59 60 61 62 ...
#   $ Rel_Count : num  42 53 648 293 121 92 519 79 765 196 ...
#   $ Rel_Pos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...

expand_feature_region <- ExpandFeatureRegionForList(gene_poscodon_codon_i200, 
                                                    dataset, 
                                                    hd_file, 
                                                    min_read_length, 
                                                    colsum_out, 
                                                    gff_df, 
                                                         gene = unique(transcript_gene_pos_poscodon_frame$Gene),
                                                         startpos, 
                                                         startlen, 
                                                         feature_of_interest, 
                                                         expand_width)

# three fuctions called in one.returns a list containing a tibble for each occurrence of the feature of interest
# TEST:: expand_feature_region should be a list
# type(expand_feature_region)
# [1] "list"

# the dimensions of each item in the list shoud be [(2*expand_width+1)X4] as there are 4 rows; Gene, Pos_Codon, Rel_Count, Rel_Pos
# and one column for each position being included ie (-5 to 5) relative to the stop codon
# TEST:: Where an expand_width of =5 is given,the following would be observed:
# > dim(expand_feature_region[[1]])
# [1] 11  4



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
  dplyr::mutate(.x, RelCount = Rel_Count / sum(Rel_Count) * (2 * expand_width + 1))
  
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
  .f = ExpandedFeatureRegionNormalization,
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




