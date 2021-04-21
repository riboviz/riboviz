
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


# IF SINGLE CODONS: USE yeast_codon_pos_i200 THROUGHOUT
# IF CODON PAIRS: USE gene_poscodon_codon_i200 THROUGHOUT (derived from yeast_codon_pos_i200)



YAL5_h5 <- here::here("Mok-simYAL5", "output", "A", "A.h5")

YAL5_gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3")

gff_df <- readGFFAsDf(YAL5_gff)



# Import the .tsv file: (used the updated file which contains the first 
# 200 codons of each gene)

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



# Configured the yeast_codon_pos_i200 to show the positions lead and lag positions 
# for each PosCodon and configured the Codon column to CodonPair so that the outputs 
# looks like this: 
# 
# > gene_poscodon_codon_i200
#   # A tibble: 2,826,757 x 3
#   Gene    PosCodon CodonPair
# <chr>   <chr>    <chr>    
#   1 YAL068C 1 2      ATG GTC  
#   2 YAL068C 2 3      GTC AAA  
# 3 YAL068C 3 4      AAA TTA  
# 4 YAL068C 4 5      TTA ACT  
# 5 YAL068C 5 6      ACT TCA  
# 6 YAL068C 6 7      TCA ATC  
# 7 YAL068C 7 8      ATC GCC  
# 8 YAL068C 8 9      GCC GCT  
# 9 YAL068C 9 10     GCT GGT  
# 10 YAL068C 10 11    GGT GTC  
# # ... with 2,826,747 more rows

gene_poscodon_codon_i200 <- tibble::tibble(
  Gene = yeast_codon_pos_i200$Gene,
  CodonPos_1 = yeast_codon_pos_i200$PosCodon, 
  CodonPos_2 = dplyr::lead(yeast_codon_pos_i200$PosCodon),
  CodonPair = paste(yeast_codon_pos_i200$Codon, dplyr::lead(yeast_codon_pos_i200$Codon))
)

# Replaced yeast_codon_pos_i200 with gene_poscodon_codon_i200 throughout the script



# Filter down the gene_poscodon_codon_i200 file to the gene that you are 
# working with, in this case YAL003W

# For codon pairs yeast_codon_pos_i200 input was changed to gene_poscodon_codon_i200

YAL003W_pos <- dplyr::filter(gene_poscodon_codon_i200, Gene=="YAL003W")
print(YAL003W_pos)



# read the structure of the h5 file 
h5ls(YAL5_h5)

# group        name       otype  dclass       dim
# 0                           /     YAL003W   H5I_GROUP                  
# 1                    /YAL003W Mok-simYAL5   H5I_GROUP                  
# 2        /YAL003W/Mok-simYAL5       reads   H5I_GROUP                  
# 3  /YAL003W/Mok-simYAL5/reads        data H5I_DATASET INTEGER 41 x 1121
# 4                           /     YAL005C   H5I_GROUP                  
# 5                    /YAL005C Mok-simYAL5   H5I_GROUP                  
# 6        /YAL005C/Mok-simYAL5       reads   H5I_GROUP                  
# 7  /YAL005C/Mok-simYAL5/reads        data H5I_DATASET INTEGER 41 x 2429
# 8                           /     YAL012W   H5I_GROUP                  
# 9                    /YAL012W Mok-simYAL5   H5I_GROUP                  
# 10       /YAL012W/Mok-simYAL5       reads   H5I_GROUP                  
# 11 /YAL012W/Mok-simYAL5/reads        data H5I_DATASET INTEGER 41 x 1685
# 12                          /     YAL035W   H5I_GROUP                  
# 13                   /YAL035W Mok-simYAL5   H5I_GROUP                  
# 14       /YAL035W/Mok-simYAL5       reads   H5I_GROUP                  
# 15 /YAL035W/Mok-simYAL5/reads        data H5I_DATASET INTEGER 41 x 3509
# 16                          /     YAL038W   H5I_GROUP                  
# 17                   /YAL038W Mok-simYAL5   H5I_GROUP                  
# 18       /YAL038W/Mok-simYAL5       reads   H5I_GROUP                  
# 19 /YAL038W/Mok-simYAL5/reads        data H5I_DATASET INTEGER 41 x 2003



### NOT DIRECTLY USED YET ###

# GetAllGeneDatamatrix: Function to get the datamatrix for each of the genes 
# contained within the h5 file

# Loop list of gene names through GetGeneDatamatrix to extract all the 
# information

# Required parameters: gene_names, dataset, hd_file)

gene_names <- rhdf5::h5ls(YAL5_h5, recursive = 1)$name

GetAllGeneDatamatrix <- function(dataset, hd_file, gene_names){
  for (gene in gene_names){
    get_datamatrix_list <- GetGeneDatamatrix(gene, dataset, hd_file)
    return(get_datamatrix_list)
  }
}

# Result: GetDatamatrixList which contains the datamatrix for all of the genes 
# contained within the list 


# Extract the datamatrix for the five genes 
Mok_YAL5_data <- GetAllGeneDatamatrix(dataset="Mok-simYAL5", hd_file = YAL5_h5, 
                               gene_names)

### ###




Reads_pos_length <- GetGeneDatamatrix(gene="YAL003W", 
                                      dataset="Mok-simYAL5", 
                                      hd_file=YAL5_h5)

tidy_gene_datamatrix <- TidyDatamatrix(GetGeneDatamatrix(gene="YAL003W", 
                                                         dataset="Mok-simYAL5", 
                                                         hd_file=YAL5_h5), 
                                       startpos = 1, 
                                       startlen = 10) 

### CalcAsiteFixed

asite_counts_by_position <- CalcAsiteFixed(Reads_pos_length, 
                        min_read_length = 10, 
                        asite_displacement_length = data.frame(read_length = c(28, 29, 30), 
                                                               asite_displacement = c(15, 15, 15)), 
                        colsum_out = TRUE)

# > str(Asite)
# num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...


TidyAsiteCountsByPosition <- function(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE){
  
  asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", "yeast_standard_asite_disp_length.txt"))
  
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  
  asite_counts_by_position <- CalcAsiteFixed(reads_pos_length, 
                                             min_read_length = min_read_length, 
                                             asite_displacement_length = asite_displacement_length,
                                             colsum_out = colsum_out)
  
  tidy_asite_counts <- tibble(Pos = 1:length(asite_counts_by_position), Count = asite_counts_by_position)
  
  return(tidy_asite_counts)
  
}


tidy_asite_count_output <- TidyAsiteCountsByPosition(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = YAL5_h5, min_read_length = 10, colsum_out = TRUE)


# > str(tidy_asite_count_output)
# tibble [1,121 x 2] (S3: tbl_df/tbl/data.frame)
# $ Pos  : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
# $ Count: num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...



TranscriptPosToCodonPos <- function(gene, gff_df){
  
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
  
  UTR5_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR5") %>% select(width)
    
  CDS_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(width)
    
  UTR3_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(width) 
  
  transcript_length <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(end)
  
  CDS_start <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start)
  
  CDS_end <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end)
  
  CDS_codon_positions <- rep(1:(CDS_width$width/3), each = 3)
  
  NA_UTR5_width <- rep(NA, times = UTR5_width$width)
  
  NA_UTR3_width <- rep(NA, times = UTR3_width$width)
  
  transcript_gene_pos_poscodon_frame <- tibble(
    Gene = gene,
    Pos = 1:transcript_length$end,
    Pos_Codon = c(rep(NA, times = UTR5_width$width), CDS_codon_positions, rep(NA, times = UTR3_width$width)),
    Frame = seq(from = 2L, to = (transcript_length$end + 1L)) %% 3 # works as UTRS are 250, might not work with UTRS of others values
    # add the count column, likely something we want here 
  )
  
  # frame filter zero 
  
  return(transcript_gene_pos_poscodon_frame)
}


transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene = "YAL003W", gff_df)
# 
# # output:
# > TranscriptPosToCodonPosOutput
# # A tibble: 1,121 x 4
# Gene      Pos Pos_Codon Frame
# <chr>   <int>     <int> <dbl>
#   1 YAL003W     1        NA     2
# 2 YAL003W     2        NA     0
# 3 YAL003W     3        NA     1
# 4 YAL003W     4        NA     2
# 5 YAL003W     5        NA     0
# 6 YAL003W     6        NA     1
# 7 YAL003W     7        NA     2
# 8 YAL003W     8        NA     0
# 9 YAL003W     9        NA     1
# 10 YAL003W    10        NA     2
# # ... with 1,111 more rows

# create a codon_per_codons_count_table


AddAsiteCountsToTranscriptPosToCodonPos <- function(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
  
  tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE)
  
  transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene, gff_df)
  
  transcript_info_tibble <- dplyr::left_join(transcript_pos_to_codon_pos_output, tidy_asite_count_output, by = "Pos")
  
  return(transcript_info_tibble)
  
}
 


transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = YAL5_h5, min_read_length = 10, colsum_out = TRUE, gff_df)


AddCodonNamesToTranscriptInfoTibble <- function(yeast_codon_pos_i200, gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
  
  codon_table <- dplyr::filter(yeast_codon_pos_i200, Gene == "YAL003W") %>% 
    dplyr::select("PosCodon", "Codon")
  
  transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df)
  
  transcript_info_tibble <- left_join(transcript_info_tibble, codon_table, by = c("Pos_Codon" = "PosCodon"), keep = FALSE)
  
  return(transcript_info_tibble)
    
}

transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(yeast_codon_pos_i200, gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = YAL5_h5, min_read_length = 10, colsum_out = TRUE, gff_df)

 

FilterForFrame <- function(transcript_info_tibble, filtering_frame = 0){
  
  transcript_info_tibble <- dplyr::filter(transcript_info_tibble, Frame == filtering_frame)
  
  return(transcript_info_tibble)
  
}

transcript_info_tibble <- FilterForFrame(transcript_info_tibble, filtering_frame = 0)

# > transcript_info_tibble
# # A tibble: 374 x 6
# Gene      Pos Pos_Codon Frame Count Codon
# <chr>   <int>     <dbl> <dbl> <dbl> <chr>
#   1 YAL003W     2        NA     0     0 NA   
# 2 YAL003W     5        NA     0     0 NA   
# 3 YAL003W     8        NA     0     0 NA   
# 4 YAL003W    11        NA     0     0 NA   
# 5 YAL003W    14        NA     0     0 NA   
# 6 YAL003W    17        NA     0     0 NA   
# 7 YAL003W    20        NA     0     0 NA   
# 8 YAL003W    23        NA     0     0 NA   
# 9 YAL003W    26        NA     0     0 NA   
# 10 YAL003W    29        NA     0     0 NA   
# # ... with 364 more rows



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



# # Function to filter for codons of interest, generates a list of the positions of the codon of interest
# FilterForCodonOfInterestPositions <- function(gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10, codon_of_interest){
#   
#   joined_gene_per_codon_counts <- NucleotideToCodonPosition(gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10)
#   
#   interesting_codon_table <- dplyr::filter(joined_gene_per_codon_counts, CodonPair == codon_of_interest) 
#   # Changed Codon == codon_of_interest to CodonPair = codon_of_interest
#     
#           # # A tibble: 2 x 4
#           #      Gene    CodonPos CodonPair PerCodonCounts
#           #      <chr>   <chr>      <chr>              <int>
#           #   1 YAL003W   7 8      TCC AAG             2313
#           #   2 YAL003W  57 58     TCC AAG             2193
#   
#   # interesting_codon_positions <- interesting_codon_table$CodonPos
#   
#   # return(interesting_codon_positions)
#   
#     # might not use these last two steps for codon pairs as there are two columns that are needed
# }
# 
# interesting_codon_positions <- FilterForCodonOfInterestPositions(gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10, codon_of_interest = "TCC AAG")

# Output from FilterForCodonOfInterestPositions for single codons
# [1]   8  22  39  58  97  99 109 110 113 116 120 128 144 166 167 178

# Output from FilterForCodonOfInterestPositions for codon pairs 
# [1] "7 8"   "57 58"


# Function FilterForCodonOfInterestPositions has been edited for codon pairs 
FilterForCodonPairOfInterestPositions <- function(gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10, codon_of_interest){

    joined_gene_per_codon_counts <- NucleotideToCodonPosition(gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10)

    interesting_codon_table <- dplyr::filter(joined_gene_per_codon_counts, CodonPair == codon_of_interest)
    # Changed Codon == codon_of_interest to CodonPair = codon_of_interest
    
    interesting_first_codon_positions <- interesting_codon_table$CodonPos_1
    
    # return(interesting_codon_positions)
    return(interesting_first_codon_positions)
}
  
interesting_first_codon_positions <- FilterForCodonPairOfInterestPositions(gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10, codon_of_interest = "TCC AAG")

# Returns: 
# > interesting_first_codon_position
# [1]  7 57
# 
# Which refers to the position of the first codon in the codon pair as the second codon
# in the pair is the first codon position + 1



# # Function needed to go from one position to a region.
# # This function will be using codon position values, not nucleotide position values
# ExpandCodonRegion <- function(.x = interesting_codon_positions, gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10, gff_df, expand_width = 5L, remove_overhang = TRUE) {
#   
#   joined_gene_per_codon_counts <- NucleotideToCodonPosition(gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10)
#   
#   gene_length <- GetGeneLength(gene, dataset, hd_file)
#   
#   #if (remove_overhang) {
#   # return an empty tibble if the desired region hangs over the edge of the coding region
#   
#   if (.x < expand_width | .x + expand_width > gene_length ) {
#     return( tibble() )
#   } else {
#     output_codon_info <- tibble(
#       dplyr::slice(joined_gene_per_codon_counts, (.x - expand_width):(.x + expand_width), each = FALSE),
#       Rel_Pos =  seq(- expand_width, expand_width)
#     )
#     return(output_codon_info)
#   }
#   #  } # if(remove_overhang)
# }
# # potential test: length of output_codon_info == (2x expand_width + 1)


# ExpandCodonRegionOutput <- ExpandCodonRegion(.x = interesting_codon_positions, gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10, gff_df, expand_width = 5L)
# 
# # Gives the output:
# 
# # > ExpandCodonRegionOutput
# # # A tibble: 11 x 5
# # Gene    CodonPos Codon PerCodonCounts Rel_Pos
# # <chr>      <dbl> <chr>          <int>   <int>
# #   1 YAL003W        3 TCC             1289      -5
# # 2 YAL003W        4 ACC              678      -4
# # 3 YAL003W        5 GAT              605      -3
# # 4 YAL003W        6 TTC              364      -2
# # 5 YAL003W        7 TCC             1154      -1
# # 6 YAL003W        8 AAG             1159       0
# # 7 YAL003W        9 ATT             1058       1
# # 8 YAL003W       10 GAA              500       2
# # 9 YAL003W       11 ACT              974       3
# # 10 YAL003W       12 TTG              199       4
# # 11 YAL003W       13 AAA             1053       5


ExpandCodonPairRegion <- function(.x = interesting_first_codon_positions, gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10, gff_df, expand_width = 5L, remove_overhang = TRUE) {
  
  joined_gene_per_codon_counts <- NucleotideToCodonPosition(gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10)
  
  gene_length <- GetGeneLength(gene, dataset, hd_file)
  
  #if (remove_overhang) {
  # return an empty tibble if the desired region hangs over the edge of the coding region
  
  if (.x < expand_width | .x + expand_width > gene_length ) {
    return( tibble() )
  } else {
    output_codon_info <- tibble(
      dplyr::slice(joined_gene_per_codon_counts, (.x - expand_width):(.x + expand_width), each = FALSE),
      Rel_Pos =  seq(- expand_width, expand_width)
    )
    return(output_codon_info)
  }
  #  } # if(remove_overhang)
}


ExpandCodonPairRegionOutput <- ExpandCodonPairRegion(.x = interesting_first_codon_positions, gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10, gff_df, expand_width = 5L, remove_overhang = TRUE)

# > ExpandCodonPairRegionOutput
# # A tibble: 11 x 6
# Gene    CodonPos_1 CodonPos_2 CodonPair PerCodonCounts Rel_Pos
# <chr>        <dbl>      <dbl> <chr>              <int>   <int>
#   1 YAL003W          2          3 GCA TCC             2977      -5
# 2 YAL003W          3          4 TCC ACC             1967      -4
# 3 YAL003W          4          5 ACC GAT             1283      -3
# 4 YAL003W          5          6 GAT TTC              969      -2
# 5 YAL003W          6          7 TTC TCC             1518      -1
# 6 YAL003W          7          8 TCC AAG             2313       0
# 7 YAL003W          8          9 AAG ATT             2217       1
# 8 YAL003W          9         10 ATT GAA             1558       2
# 9 YAL003W         10         11 GAA ACT             1474       3
# 10 YAL003W         11         12 ACT TTG             1173       4
# 11 YAL003W         12         13 TTG AAA             1252       5




# Apply the ExpandCodonRegion function to the codons of interest to generate expanded tibbles for each position 

ExpandCodonPairRegionForList <- function(gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10, codon_of_interest, gff_df, expand_width = 5L){
  
  interesting_first_codon_positions <- FilterForCodonPairOfInterestPositions(gene_poscodon_codon_i200, 
                                        gene, dataset, hd_file, startpos = 1, startlen = 10, codon_of_interest)
  
  ExpandList <- purrr::map(
    # .x = interesting_codon_positions, # vector of codon positions to iterate over for single codons
    .x = interesting_first_codon_positions,
    # .f = ExpandCodonRegion,   # function to use at each codon position of interest for single codons
    .f = ExpandCodonPairRegion,
    gene_poscodon_codon_i200,
    gene,
    dataset,
    hd_file,
    startpos = 1,
    startlen = 10,
    gff_df,
    expand_width = 5L
  )
  
}

ExpandList <- ExpandCodonPairRegionForList(gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10, codon_of_interest = "TCC AAG", gff_df, expand_width = 5L)

# > ExpandList
# [[1]]
# # A tibble: 11 x 6
# Gene    CodonPos_1 CodonPos_2 CodonPair PerCodonCounts Rel_Pos
# <chr>        <dbl>      <dbl> <chr>              <int>   <int>
#   1 YAL003W          2          3 GCA TCC             2977      -5
# 2 YAL003W          3          4 TCC ACC             1967      -4
# 3 YAL003W          4          5 ACC GAT             1283      -3
# 4 YAL003W          5          6 GAT TTC              969      -2
# 5 YAL003W          6          7 TTC TCC             1518      -1
# 6 YAL003W          7          8 TCC AAG             2313       0
# 7 YAL003W          8          9 AAG ATT             2217       1
# 8 YAL003W          9         10 ATT GAA             1558       2
# 9 YAL003W         10         11 GAA ACT             1474       3
# 10 YAL003W         11         12 ACT TTG             1173       4
# 11 YAL003W         12         13 TTG AAA             1252       5
# 
# [[2]]
# # A tibble: 11 x 6
# Gene    CodonPos_1 CodonPos_2 CodonPair PerCodonCounts Rel_Pos
# <chr>        <dbl>      <dbl> <chr>              <int>   <int>
#   1 YAL003W         52         53 TTC AAC             2545      -5
# 2 YAL003W         53         54 AAC CAC             2480      -4
# 3 YAL003W         54         55 CAC ATC             3354      -3
# 4 YAL003W         55         56 ATC GCT             2510      -2
# 5 YAL003W         56         57 GCT TCC             1678      -1
# 6 YAL003W         57         58 TCC AAG             2193       0
# 7 YAL003W         58         59 AAG GCC             1315       1
# 8 YAL003W         59         60 GCC GAT             2335       2
# 9 YAL003W         60         61 GAT GAA             3061       3
# 10 YAL003W         61         62 GAA TTC             1495       4
# 11 YAL003W         62         63 TTC GAC              889       5



# Normalization carried out within each expanded frame so that they are comparable 
# Normalizes the ExpandCodonRegion list generating a RelCount column with the normalization values
ExpandedCodonRegionNormalization <- function(.x, expand_width = 5L){
  
  dplyr::mutate(.x, RelCount = PerCodonCounts / sum(PerCodonCounts) * (2 * expand_width + 1))
  
}


table_a <- ExpandedCodonRegionNormalization(ExpandCodonPairRegionOutput, expand_width = 5L)

# # A tibble: 11 x 6
# Gene    CodonPos Codon PerCodonCounts Rel_Pos RelCount
# <chr>      <dbl> <chr>          <int>   <int>    <dbl>
#   1 YAL003W        3 TCC             1289      -5    1.57 
# 2 YAL003W        4 ACC              678      -4    0.826
# 3 YAL003W        5 GAT              605      -3    0.737
# 4 YAL003W        6 TTC              364      -2    0.443
# 5 YAL003W        7 TCC             1154      -1    1.41 
# 6 YAL003W        8 AAG             1159       0    1.41 
# 7 YAL003W        9 ATT             1058       1    1.29 
# 8 YAL003W       10 GAA              500       2    0.609
# 9 YAL003W       11 ACT              974       3    1.19 
# 10 YAL003W       12 TTG              199       4    0.242
# 11 YAL003W       13 AAA             1053       5    1.28 


# > table_a <- ExpandedCodonRegionNormalization(ExpandCodonPairRegionOutput, expand_width = 5L)
# > table_a
# # A tibble: 11 x 7
# Gene    CodonPos_1 CodonPos_2 CodonPair PerCodonCounts Rel_Pos RelCount
# <chr>        <dbl>      <dbl> <chr>              <int>   <int>    <dbl>
#   1 YAL003W          2          3 GCA TCC             2977      -5    1.75 
# 2 YAL003W          3          4 TCC ACC             1967      -4    1.16 
# 3 YAL003W          4          5 ACC GAT             1283      -3    0.755
# 4 YAL003W          5          6 GAT TTC              969      -2    0.570
# 5 YAL003W          6          7 TTC TCC             1518      -1    0.893
# 6 YAL003W          7          8 TCC AAG             2313       0    1.36 
# 7 YAL003W          8          9 AAG ATT             2217       1    1.30 
# 8 YAL003W          9         10 ATT GAA             1558       2    0.916
# 9 YAL003W         10         11 GAA ACT             1474       3    0.867
# 10 YAL003W         11         12 ACT TTG             1173       4    0.690
# 11 YAL003W         12         13 TTG AAA             1252       5    0.736



# Normalization carried out for all the tibbles within ExpandList 
NormalizedExpandList <- purrr::map(
  .x = ExpandList,
  .f = ExpandedCodonRegionNormalization,
  expand_width = 5L
)



# Function to generate a plot for the expanded list of the codon
GenerateGraphs <- function(.x){
  
  plot <- ggplot(.x, mapping = aes(x = Rel_Pos, y = RelCount)) + geom_line()
  
}

plot <- GenerateGraphs(.x = table_a)




# Function to generate a plot for each of the tibbles within NormalizedExpandList
GenerateAllGraphs <- purrr::map(
  .x = NormalizedExpandList,
  .f = GenerateGraphs
)

# Generates a list consisting of all of the graphs, by running e.g. "GenerateAllGraphs[[2]]" 
# from the console you can open its graph 



# Function to overlay graphs into a single graph.Need to generate a single tibble 
# from NormalizedExpandList. Need to join by Rel_Pos, in RelCount need the mean for 
# each Rel_Pos (sum row(x) / number of row(x))






Overlayed <- function(NormalizedExpandList, expand_width = 5L){
  number_of_objects <- length(NormalizedExpandList)
  
  result = lapply(NormalizedExpandList, "[", c("Rel_Pos", "RelCount"))
  
  joined_result = result %>% reduce(full_join, by = c("Rel_Pos"), sum("RelCount"))
  
  Overlayed_tibbles <- tibble::tibble(
    Rel_Pos = seq(- expand_width, expand_width),
    RelCount = sum(joined_result["RelCount.y", "RelCount.x"])/number_of_objects 
  )
}


Over <- Overlayed(NormalizedExpandList, expand_width = 5L)  



  
  # input : NormalizedExpandList - list of tibbles for each occurence of the codon pair 
  # > NormalizedExpandList
  
  # # A tibble: 11 x 7
  # Gene    CodonPos_1 CodonPos_2 CodonPair PerCodonCounts Rel_Pos RelCount
  # <chr>        <dbl>      <dbl> <chr>              <int>   <int>    <dbl>
  #   1 YAL003W          2          3 GCA TCC             2977      -5    1.75 
  # 2 YAL003W          3          4 TCC ACC             1967      -4    1.16 
  # 3 YAL003W          4          5 ACC GAT             1283      -3    0.755
  # 4 YAL003W          5          6 GAT TTC              969      -2    0.570
  # 5 YAL003W          6          7 TTC TCC             1518      -1    0.893
  # 6 YAL003W          7          8 TCC AAG             2313       0    1.36 
  # 7 YAL003W          8          9 AAG ATT             2217       1    1.30 
  # 8 YAL003W          9         10 ATT GAA             1558       2    0.916
  # 9 YAL003W         10         11 GAA ACT             1474       3    0.867
  # 10 YAL003W         11         12 ACT TTG             1173       4    0.690
  # 11 YAL003W         12         13 TTG AAA             1252       5    0.736
  # 
  # [[2]]
  # # A tibble: 11 x 7
  # Gene    CodonPos_1 CodonPos_2 CodonPair PerCodonCounts Rel_Pos RelCount
  # <chr>        <dbl>      <dbl> <chr>              <int>   <int>    <dbl>
  #   1 YAL003W         52         53 TTC AAC             2545      -5    1.17 
  # 2 YAL003W         53         54 AAC CAC             2480      -4    1.14 
  # 3 YAL003W         54         55 CAC ATC             3354      -3    1.55 
  # 4 YAL003W         55         56 ATC GCT             2510      -2    1.16 
  # 5 YAL003W         56         57 GCT TCC             1678      -1    0.774
  # 6 YAL003W         57         58 TCC AAG             2193       0    1.01 
  # 7 YAL003W         58         59 AAG GCC             1315       1    0.606
  # 8 YAL003W         59         60 GCC GAT             2335       2    1.08 
  # 9 YAL003W         60         61 GAT GAA             3061       3    1.41 
  # 10 YAL003W         61         62 GAA TTC             1495       4    0.689
  # 11 YAL003W         62         63 TTC GAC              889       5    0.410
  
  



# how to apply to all genes?
# make some plots 
    # write some code to make a plot for a summarised version 



