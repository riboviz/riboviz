
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




# GetGeneDatamatrix for each of the genes contained within the Mok-simYAL5
#
# YAL003W_datamatrix <- GetGeneDatamatrix(gene="YAL003W", dataset="Mok-simYAL5",
#                                  hd_file = YAL5_h5) #data_mat
#
# # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...
# # where 1:41 represents read lengths 10:50 and 1:1121 represents the 1121 
# # nucleotides including UTRs and CDS
# 
# 
# # The GetGeneDatamatrix contains the CDS as well as the UTRs, results in 1121 
# # observations instead of the 207 which consists of the CDS
# 
# # YAL005C_data <- GetGeneDatamatrix(gene="YAL005C", dataset="Mok-simYAL5", 
# #                                   hd_file = YAL5_h5)
# # YAL012W_data <- GetGeneDatamatrix(gene="YAL012W", dataset="Mok-simYAL5", 
# #                                   hd_file = YAL5_h5)
# # YAL035W_data <- GetGeneDatamatrix(gene="YAL035W", dataset="Mok-simYAL5", 
# #                                   hd_file = YAL5_h5)
# # YAL038W_data <- GetGeneDatamatrix(gene="YAL038W", dataset="Mok-simYAL5", 
# #                                   hd_file = YAL5_h5)
# 
# 
# 
# # Make a table which contains the datamatrix from the h5 file, in this case only 
# # carried out for a single gene YAL003W
# 
# 
# 
#                                         
#   #   reads_pos_length <- GetGeneDatamatrix(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = YAL5_h5)
#   #           print(reads_pos_length)
#   #   reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE) 
#   #           print(reads_asitepos)
#   #   SumByFrame(reads_asitepos, left =251, right =871)
#   #   SnapToCodonList <- SnapToCodon(reads_asitepos, left=251, right=871, snapdisp = 0L)                                         
#   #                                        
#   # print(SnapToCodonList)                                       
# 
# # buffer <- 250
# # min_read_length <- 10
# # Codon_Pos_Reads <- GetCodonPositionReads(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = YAL5_h5, left = (buffer - 15), right = (buffer + 11), min_read_length = min_read_length)
# 
# 
# YAL003W_tidy <- TidyDatamatrix(data_mat = YAL003W_datamatrix, startpos = 1, startlen=10)
# 
# # A tibble: 45,961 x 3
# # ReadLen   Pos Counts
# # <int> <int>  <int>
# #   1      10     1      0
# # 2      11     1      0
# # 3      12     1      0
# # 4      13     1      0
# # 5      14     1      0
# # 6      15     1      0
# # 7      16     1      0
# # 8      17     1      0
# # 9      18     1      0
# # 10     19     1      0
# # # ... with 45,951 more rows
# 
# # Tibble contains the read lengths 10:50 for position 1-1121 
# # The number of observations (45,961) is equal to 41 read lengths multiplied by 1121 nucleotides
# 
# # Summaried the table so that the total count is shown for each of the positions 
# # summary_counts_for_positions <- summarise(group_by(YAL003W_tidy, Pos), Total_counts = sum(Counts))
# # A tibble: 1,121 x 2
# # Pos Total_counts
# # * <int>        <int>
# #   1     1            0
# #   2     2            0
# #   3     3            0
# #   4     4            0
# #   5     5            0
# #   6     6            0
# #   7     7            0
# #   8     8            0
# #   9     9            0
# #  10    10            0
# # ... with 1,111 more rows
# 
# # Tidydatamatrix is grouped by position, removes read lengths
# 
# 
# # start_codon_pos <- GetGeneStartCodonPos(gene="YAL003W", dataset="Mok-simYAL5", # 251 252 253
# #                      hd_file = YAL5_h5) 
# # 
# # stop_codon_pos <- GetGeneStopCodonPos(gene="YAL003W", dataset="Mok-simYAL5",  # 869, 870, 871
# #                     hd_file=YAL5_h5)
# 
# 
# # CDS start 251, CDS end 871 - CDS total: 621 nt -> 207 codons
# # 250 UTRs at each end for a total of 1121 nt



# KeepGeneCDS function: remove the UTRs at the 5' and 3' end, keep the CDS

KeepGeneCDS <- function(gene, dataset, hd_file, startpos = 1, startlen = 10){
  
  tidy_gene_datamatrix <- TidyDatamatrix(GetGeneDatamatrix(gene, dataset, hd_file), startpos, startlen) 
  
  summary_counts_for_positions <- dplyr::summarise(dplyr::group_by(tidy_gene_datamatrix, Pos), Total_counts = sum(Counts))
  
  start_codon_pos <- GetGeneStartCodonPos(gene, dataset, hd_file)
  stop_codon_pos <- GetGeneStopCodonPos(gene, dataset, hd_file)
  
  summary_counts_for_positions_CDS <- summary_counts_for_positions %>%
    dplyr::filter(dplyr::between(Pos, min(start_codon_pos), max(stop_codon_pos))) 

  return(summary_counts_for_positions_CDS)
  
}


gene_CDS <- KeepGeneCDS(gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10)

# KeepGeneCDS generates as an output:

# A tibble: 621 x 2
# Pos Total_counts
# <int>        <int>
#   1   251          225
#   2   252           90
#   3   253           77
#   4   254          246
#   5   255          291
#   6   256         1151
#   7   257         1272
#   8   258            8
#   9   259            9
#  10   260           25
# ... with 611 more rows




# NucleotideToCodonPosition function: change format from nucleotides to codons,
# combine the generated codon table and the yeast_codon_pos_i200 file  

NucleotideToCodonPosition <- function(gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10){
  
  codon_table <- dplyr::filter(gene_poscodon_codon_i200, Gene==gene)
  
  gene_CDS <- KeepGeneCDS(gene, dataset, hd_file, startpos = 1, startlen = 10)
  
  # per_codon_counts <- zoo::rollapply(data=gene_CDS$Total_counts, width = 3, sum, by = 3)
  
  per_codonpair_counts <- zoo::rollapply(data=gene_CDS$Total_counts, width = 6, sum, by = 3)
  # configured the code for per_codon_counts so that it sums 1+2, 2+3, 3+4...
  
  gene_per_codon_counts <- tibble::tibble(
    # CodonPos = paste(seq_len(length(per_codon_counts)), 
    CodonPos_1 = seq_len(length(per_codonpair_counts)),
    CodonPos_2 = seq_len(length(per_codonpair_counts)) +1,
        # Line 304 was replaced, so CodonPos became CodonPos_1 and CodonPos_2
    PerCodonCounts = per_codonpair_counts
  )
  
  joined_gene_per_codon_counts <- dplyr::full_join(
    x = gene_per_codon_counts, 
    y = codon_table,  
    by = c("CodonPos_1" = "CodonPos_1", "CodonPos_2" = "CodonPos_2")
  ) %>% 
    dplyr::select(Gene, CodonPos_1, CodonPos_2, CodonPair, PerCodonCounts) 
  
  
  return(joined_gene_per_codon_counts)
  
}


joined_gene_per_codon_counts_YAL003W <- NucleotideToCodonPosition(gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10)


# Output generated from NucleotideToCodonPosition function for single codons:
# 
# # A tibble: 207 x 4
#     Gene    CodonPos Codon PerCodonCounts
#     <chr>      <dbl> <chr>          <int>
#   1 YAL003W        1 ATG              392
#   2 YAL003W        2 GCA             1688
#   3 YAL003W        3 TCC             1289
#   4 YAL003W        4 ACC              678
#   5 YAL003W        5 GAT              605
#   6 YAL003W        6 TTC              364
#   7 YAL003W        7 TCC             1154
#   8 YAL003W        8 AAG             1159
#   9 YAL003W        9 ATT             1058
#  10 YAL003W       10 GAA              500
# # ... with 197 more rows
# 
# 
# Output generated from NucleotideToCodonPosition function for codon pairs:
# 
#   # A tibble: 207 x 4
#   Gene    CodonPos CodonPair PerCodonCounts
# <chr>   <chr>    <chr>              <int>
#   1 YAL003W 1 2      ATG GCA             2080
# 2 YAL003W 2 3      GCA TCC             2977
# 3 YAL003W 3 4      TCC ACC             1967
# 4 YAL003W 4 5      ACC GAT             1283
# 5 YAL003W 5 6      GAT TTC              969
# 6 YAL003W 6 7      TTC TCC             1518
# 7 YAL003W 7 8      TCC AAG             2313
# 8 YAL003W 8 9      AAG ATT             2217
# 9 YAL003W 9 10     ATT GAA             1558
# 10 YAL003W 10 11    GAA ACT             1474
# # ... with 197 more rows


# To solve problem caused by the output from FilterForCodonOfInterestPositions 
# being in string form I splot CodonPos into two columns consisting of CodonPos_1
# and CodonPos_2:
#   
#   # A tibble: 207 x 5
#   Gene    CodonPos_1 CodonPos_2 CodonPair PerCodonCounts
# <chr>        <dbl>      <dbl> <chr>              <int>
#   1 YAL003W          1          2 ATG GCA             2080
# 2 YAL003W          2          3 GCA TCC             2977
# 3 YAL003W          3          4 TCC ACC             1967
# 4 YAL003W          4          5 ACC GAT             1283
# 5 YAL003W          5          6 GAT TTC              969
# 6 YAL003W          6          7 TTC TCC             1518
# 7 YAL003W          7          8 TCC AAG             2313
# 8 YAL003W          8          9 AAG ATT             2217
# 9 YAL003W          9         10 ATT GAA             1558
# 10 YAL003W         10         11 GAA ACT             1474
# # ... with 197 more rows



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
}
  
interesting_codonpair_positions <- FilterForCodonPairOfInterestPositions(gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10, codon_of_interest = "TCC AAG")


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


ExpandCodonPairRegion <- function(.x = interesting_codonpair_positions, gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10, gff_df, expand_width = 5L, remove_overhang = TRUE) {
  
  joined_gene_per_codon_counts <- NucleotideToCodonPosition(gene_poscodon_codon_i200, gene, dataset, hd_file, startpos = 1, startlen = 10)
  
  gene_length <- GetGeneLength(gene, dataset, hd_file)
  
  #if (remove_overhang) {
  # return an empty tibble if the desired region hangs over the edge of the coding region
  
  if (.x$CodonPos_1 < expand_width | .x$CodonPos_2 + expand_width > gene_length ) {
    return( tibble() )
  } else {
    output_codon_info <- tibble(
      dplyr::slice(joined_gene_per_codon_counts, (.x$CodonPos_1 - expand_width):(.x$CodonPos_2 + expand_width -1), each = FALSE),
      Rel_Pos =  seq(- expand_width, expand_width)
    )
    return(output_codon_info)
  }
  #  } # if(remove_overhang)
}


ExpandCodonPairRegionOutput <- ExpandCodonPairRegion(.x = interesting_codonpair_positions, gene_poscodon_codon_i200, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10, gff_df, expand_width = 5L, remove_overhang = TRUE)






# Apply the ExpandCodonRegion function to the codons of interest to generate expanded tibbles for each position 

ExpandList <- purrr::map(
  # .x = interesting_codon_positions, # vector of codon positions to iterate over for single codons
  .x = interesting_codonpair_positions,
  # .f = ExpandCodonRegion,   # function to use at each codon position of interest for single codons
  .f = ExpandCodonPairRegion(),
  gene_poscodon_codon_i200,
  gene="YAL003W",
  dataset="Mok-simYAL5",
  hd_file=YAL5_h5,
  startpos = 1,
  startlen = 10,
  gff_df,
  expand_width = 5L
)


# error:
#   Error in NucleotideToCodonPosition(gene_poscodon_codon_i200, gene, dataset,  : 
#                                        argument "gene_poscodon_codon_i200" is missing, with no default 



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

# Generates a list consisting of all of the graphs, by running e.g. "GenerateAllGraphs[[3]]" 
# from the console you can open its graph 




# how to apply to all genes?
# make some plots 
    # write some code to make a plot for a summarised version 

# when does overlay happen 





################################################################################

# Need to do asite assignment. This will likley need to be carried out either as
# part of the tidydatamatrix function or right after? 


Reads_pos_length <- GetGeneDatamatrix(gene="YAL003W", 
                                      dataset="Mok-simYAL5", 
                                      hd_file=YAL5_h5)

tidy_gene_datamatrix <- TidyDatamatrix(GetGeneDatamatrix(gene="YAL003W", 
                                                         dataset="Mok-simYAL5", 
                                                         hd_file=YAL5_h5), 
                                       startpos = 1, 
                                       startlen = 10) 

# format generated from the tidy_gene_datamatrix - either need asite assignment to
# occur within that function or to give an output in the same format?

# > tidy_gene_datamatrix
# # A tibble: 45,961 x 3
# ReadLen   Pos Counts
# <int> <int>  <int>
#   1      10     1      0
# 2      11     1      0
# 3      12     1      0
# 4      13     1      0
# 5      14     1      0
# 6      15     1      0
# 7      16     1      0
# 8      17     1      0
# 9      18     1      0
# 10      19     1      0
# # ... with 45,951 more rows

matrix <- GetGeneDatamatrix(gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5)




### CalcAsiteFixedOneLength

fixed_one_length <- CalcAsiteFixedOneLength(Reads_pos_length, 
                                            min_read_length = 10, 
                                            read_length = 10:50, 
                                            asite_displacement = 15)

tidy_gene_datamatrix_fixed <- TidyDatamatrix(fixed_one_length, 
                                             startpos = 1, 
                                             startlen = 10)

# > tidy_gene_datamatrix_fixed
# # A tibble: 45,961 x 3
# ReadLen   Pos Counts
# <int> <int>  <int>
#   1      10     1      0
# 2      11     1      0
# 3      12     1      0
# 4      13     1      0
# 5      14     1      0
# 6      15     1      0
# 7      16     1      0
# 8      17     1      0
# 9      18     1      0
# 10      19     1      0
# # ... with 45,951 more rows

# Seems to be compatible with tidydatamatrix - might be the best option?



### CalcAsiteFixed

Asite <- CalcAsiteFixed(Reads_pos_length, 
                        min_read_length = 10, 
                        asite_displacement_length = data.frame(read_length = c(28, 29, 30), 
                                                               asite_displacement = c(15, 15, 15)), 
                        colsum_out = TRUE)

tidy_gene_datamatrix_asite <- TidyDatamatrix(Asite, startpos = 1, startlen = 10) 

# Error in startpos:(startpos + ncol(data_mat) - 1) : argument of length 0



### SnapToCodon

# reads_pos_length <- GetGeneDatamatrix(gene = "YAL003W", 
#                                       dataset = "vignette", 
#                                       hd_file = "vignette/output/WTnone/WTnone.h5")
#'  # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...

reads_asitepos <- CalcAsiteFixed(Reads_pos_length, 
                                 min_read_length = 10, 
                                 asite_displacement_length = data.frame(read_length = c(28, 29, 30), 
                                                                        asite_displacement = c(15, 15, 15)), 
                                 colsum_out = TRUE) 

#'  # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
SnapToCodon(reads_asitepos, left=251, right=871, snapdisp=0L)
#'  # num [1:207] 0 9 0 0 0 3 0 0 0 0 ...

# > SnapToCodon(reads_asitepos, left=251, right=871, snapdisp=0L)
# [1] 4249  825 1017 1176 1116  284 1553  776  607  491  260  876  837  648  421  739  149
# [18]  941  514  603 1442 1101 1917  233  623 1289  622  660   94  130  567  279 1004  328
# [35]  283  625  218  832  592  287  549  456  366  594 1051 1659  843  463  885 4259 1690
# [52]  437  857 1050  807  373 1147  816 1085 1312  318  982  757  277 1575  924  330  364
# [69]  440  409  427  282  843 1262  820 1024  781 1083  727  358  877  416  354 1824  667
# [86] 1092  992 1531  923  277 1072  464 1010  349 1035 1165  664  933  530  773 1098  521
# [103]   99 1055 1005  567  460  826 1460 1030  399 1792 1454 1069  679 1294  747  355 1188
# [120]  599  172  225  301  694 1020  950  977 1904 4744 1806  779 1083  805  554  245  887
# [137]  743  674  906  404  221  213  978 1340  285  666  725  765 1421 1171  843 2475 3188
# [154]  406  443  498  996  198 1104  900 1421  670  687 1891 1103 1514  311  953  199  314
# [171]  482  306  106  992 1618 1433 1090  547  469  503  797 1160  676  350 1069  924  452
# [188] 1161  693 2231 1062  938 1335  371  654  294  453 1096  274  357  489  508  955 1287
# [205] 1496   75    0


# Potential issue: the values generated from SnapToCodon does not match the values 
# generated from NucleotideToCodonPosition - why? 
# These values should be the same as far as I can tell - is one or the other taking 
# an extra input/output?


NucleotideToCodonPosition <- function(gene_poscodon_codon_i200, 
                                      gene="YAL003W", 
                                      dataset="Mok-simYAL5", 
                                      hd_file=YAL5_h5, 
                                      startpos = 1, 
                                      startlen = 10){
  
  codon_table <- dplyr::filter(gene_poscodon_codon_i200, Gene=="YAL003W")
  
  gene_CDS <- KeepGeneCDS(gene="YAL003W", 
                          dataset="Mok-simYAL5", 
                          hd_file=YAL5_h5, 
                          startpos = 1, 
                          startlen = 10)
  
  per_codon_counts <- zoo::rollapply(data=gene_CDS$Total_counts, 
                                     width = 3, sum, by = 3)
  
  # > per_codon_counts
  # [1]  392 1688 1289  678  605  364 1154 1159 1058  500  974  199 1053  634  690 1685 1420
  # [18] 2454  468  744 1659  918  822  136  146  681  384 1231  541  334  833  272  974  835
  # [35]  351  773  602  427  771 1416 2031  999  660 1037 5206 2951  473 1032 1514  967  451
  # [52] 1369 1176 1304 2050  460 1218  975  340 1995 1066  429  460  575  495  511  357  971
  # [69] 1581 1086 1374 1103 1302  957  475  954  587  428 1877  984 1360 1089 1716 1078  383
  # [86] 1217  649 1182  486 1181 1369  955 1138  713  912 1274  762  118 1122 1562  742  581
  # [103] 1006 1915 1339  442 2237 1871 1393  760 1690 1141  471 1348  855  284  242  388  850
  # [120] 1383 1245 1143 2309 6613 2485  984 1412 1027  691  299 1130  819  853 1359  494  296
  # [137]  266 1105 2010  353  777 1280 1041 1804 1596 1228 2730 5230  560  524  665 1394  226
  # [154] 1428 1154 1717  871  773 2890 1327 2147  431 1159  289  386  642  442  112 1164 2238
  # [171] 1783 1311  795  533  619  934 1570  898  424 1339 1128  662 1354  833 2570 1312 1229
  # [188] 1881  500  745  473  520 1317  454  412  568  668 1298 1528 2133  132    0    0    0
  # [205]    0    0    0
  
  per_codonpair_counts <- zoo::rollapply(data=gene_CDS$Total_counts, width = 6, sum, by = 3)
  # configured the code for per_codon_counts so that it sums 1+2, 2+3, 3+4...
  
  # > per_codon_counts
  # [1]  392 1688 1289  678  605  364 1154 1159 1058  500  974  199 1053  634  690 1685 1420
  # [18] 2454  468  744 1659  918  822  136  146  681  384 1231  541  334  833  272  974  835
  # [35]  351  773  602  427  771 1416 2031  999  660 1037 5206 2951  473 1032 1514  967  451
  # [52] 1369 1176 1304 2050  460 1218  975  340 1995 1066  429  460  575  495  511  357  971
  # [69] 1581 1086 1374 1103 1302  957  475  954  587  428 1877  984 1360 1089 1716 1078  383
  # [86] 1217  649 1182  486 1181 1369  955 1138  713  912 1274  762  118 1122 1562  742  581
  # [103] 1006 1915 1339  442 2237 1871 1393  760 1690 1141  471 1348  855  284  242  388  850
  # [120] 1383 1245 1143 2309 6613 2485  984 1412 1027  691  299 1130  819  853 1359  494  296
  # [137]  266 1105 2010  353  777 1280 1041 1804 1596 1228 2730 5230  560  524  665 1394  226
  # [154] 1428 1154 1717  871  773 2890 1327 2147  431 1159  289  386  642  442  112 1164 2238
  # [171] 1783 1311  795  533  619  934 1570  898  424 1339 1128  662 1354  833 2570 1312 1229
  # [188] 1881  500  745  473  520 1317  454  412  568  668 1298 1528 2133  132    0    0    0
  # [205]    0    0    0
  
  
  
  
  ### GetGeneCodonPosReads1dsnap
  
  GetGeneCodonPosReads1dsnap(gene = "YAL003W", 
                             dataset = "Mok-simYAL5", 
                             hd_file = YAL5_h5, 
                             left=251, right=871, 
                             min_read_length=10, 
                             asite_displacement_length = data.frame(read_length = c(28, 29, 30), 
                                                                    asite_displacement = c(15, 15, 15)), 
                             snapdisp=0L)
  
  # > GetGeneCodonPosReads1dsnap(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = YAL5_h5, 
  #   left=251, right=871, min_read_length=10, asite_displacement_length = data.frame(read_length 
  #   = c(28, 29, 30), asite_displacement = c(15, 15, 15)), snapdisp=0L)
  # [1] 4249  825 1017 1176 1116  284 1553  776  607  491  260  876  837  648  421  739  149
  # [18]  941  514  603 1442 1101 1917  233  623 1289  622  660   94  130  567  279 1004  328
  # [35]  283  625  218  832  592  287  549  456  366  594 1051 1659  843  463  885 4259 1690
  # [52]  437  857 1050  807  373 1147  816 1085 1312  318  982  757  277 1575  924  330  364
  # [69]  440  409  427  282  843 1262  820 1024  781 1083  727  358  877  416  354 1824  667
  # [86] 1092  992 1531  923  277 1072  464 1010  349 1035 1165  664  933  530  773 1098  521
  # [103]   99 1055 1005  567  460  826 1460 1030  399 1792 1454 1069  679 1294  747  355 1188
  # [120]  599  172  225  301  694 1020  950  977 1904 4744 1806  779 1083  805  554  245  887
  # [137]  743  674  906  404  221  213  978 1340  285  666  725  765 1421 1171  843 2475 3188
  # [154]  406  443  498  996  198 1104  900 1421  670  687 1891 1103 1514  311  953  199  314
  # [171]  482  306  106  992 1618 1433 1090  547  469  503  797 1160  676  350 1069  924  452
  # [188] 1161  693 2231 1062  938 1335  371  654  294  453 1096  274  357  489  508  955 1287
  # [205] 1496   75    0
  
  
  
  ### GetGeneReadFrame
  
  GetGeneReadFrame(gene = "YAL003W", 
                   dataset = "Mok-simYAL5", 
                   hd_file = YAL5_h5, 
                   left=251, right=871, 
                   min_read_length=10, 
                   asite_displacement_length = data.frame(read_length = c(28, 29, 30), 
                                                          asite_displacement = c(15, 15, 15)))
  
  # # A tibble: 1 x 7
  # gene    Ct_fr0 Ct_fr1 Ct_fr2 pval_fr0vs1 pval_fr0vs2 pval_fr0vsboth
  # <chr>    <dbl>  <dbl>  <dbl>       <dbl>       <dbl>          <dbl>
  #   1 YAL003W  61847  50043  65120       0.160       0.943          0.659
  
  
  # Not likely to be a useful function as it does not give you the reads 