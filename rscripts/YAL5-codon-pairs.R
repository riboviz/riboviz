
# YAL5_h5 is at location $HOME/riboviz/riboviz/Mok-simYAL5/output/A/A.h5

# source packages and functions from rscripts 
source(here::here("rscripts", "read_count_functions.R"))
source(here::here("rscripts", "stats_figs_block_functions.R"))


# other packages I loaded that I used 
library(tidyverse)
# install.packages("zoo") - new package required for function rollapply
library(zoo)

YAL5_h5 <- here::here("Mok-simYAL5", "output", "A", "A.h5")

YAL5_gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3")

gff_df <- readGFFAsDf(YAL5_gff)

# Import the .tsv file: (used the updated file which contains the first 
# 200 codons of each gene)

yeast_codon_pos_i200.tsv <- readr::read_tsv(file = here::here("data", "yeast_codon_table.tsv"))

# Confirm that the .tsv file has been imported and check structure:
head(yeast_codon_pos_i200.tsv) 



# Filter down the yeast_codon_pos_i200.tsv file to the gene that you are 
# working with, in this case YAL003W

YAL003W_pos <- dplyr::filter(yeast_codon_pos_i200.tsv, Gene=="YAL003W")
print(YAL003W_pos)



# read the structure of the h5 file 
h5ls(YAL5_h5)



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




# # GetGeneDatamatrix for each of the genes contained within the Mok-simYAL5
# 
# YAL003W_datamatrix <- GetGeneDatamatrix(gene="YAL003W", dataset="Mok-simYAL5", 
#                                   hd_file = YAL5_h5) #data_mat
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
# # 2     2            0
# # 3     3            0
# # 4     4            0
# # 5     5            0
# # 6     6            0
# # 7     7            0
# # 8     8            0
# # 9     9            0
# # 10    10            0
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

# A tibble: 621 x 2
# Pos Total_counts
# <int>        <int>
#   1   251          225
# 2   252           90
# 3   253           77
# 4   254          246
# 5   255          291
# 6   256         1151
# 7   257         1272
# 8   258            8
# 9   259            9
# 10   260           25
# ... with 611 more rows



# Make final table which contains Gene, Pos, Codon, Counts  
CodonTable=c()
CodonTable$Gene<-YAL003W_pos$Gene
CodonTable$Pos<-YAL003W_pos$PosCodon
CodonTable$Codon<-YAL003W_pos$Codon
Codon_table <- as_tibble(CodonTable)
print(Codon_table)



# Codon_table needs to be a function that creates a codon table for each gene, rename so that it is more consistent

NucleotideToCodonPosition <- function(Codon_table, gene, dataset, hd_file, startpos = 1, startlen = 10){
  
  gene_CDS <- KeepGeneCDS(gene, dataset, hd_file, startpos = 1, startlen = 10)
  
  per_codon_counts <- zoo::rollapply(data=gene_CDS$Total_counts, width = 3, sum, by = 3)
  
  gene_per_codon_counts <- tibble::tibble(
    CodonPos = seq_len(length(per_codon_counts)), 
    PerCodonCounts = per_codon_counts
  )
  
  joined_gene_per_codon_counts <- dplyr::full_join(
    x = gene_per_codon_counts, 
    y = Codon_table,  
    by = c("CodonPos" = "Pos")
  ) %>% 
    dplyr::select(Gene, CodonPos, Codon, PerCodonCounts) 
  
  
  return(joined_gene_per_codon_counts)
  
}


joined_gene_per_codon_counts_YAL003W <- NucleotideToCodonPosition(Codon_table, gene="YAL003W", dataset="Mok-simYAL5", hd_file=YAL5_h5, startpos = 1, startlen = 10)




# Filter down to the desired codon 
AAG_table <- dplyr::filter(joined_gene_per_codon_counts_YAL003W, Codon=="AAG")
# AAG_df<-as.data.frame(AAG_table) worth investigating removing this 


# Normalization - when? 
dplyr::mutate(RelCount = Count / sum(Count) * ( 2 * width + 1))
# goes from Pos_Codon - width to Pos_Codon + width



# Function needed to go from one position to a region.
ExpandCodonRegion <- function(gene, Pos_Codon, gff_df, width = 5L, remove_overhang = TRUE) {
  gene_length <- GetGeneLength(gene, dataset, hd_file)
  if (remove_overhang) {
    # return an empty tibble if the desired region hangs over the edge of the coding region
    if (Pos_codon < width | Pos_codon + width > gene_length ) {
      return( tibble() )
    }
    tibble(Gene      = rep(Gene, 2 * width ),
           Pos_Codon = seq(Pos_Codon - width, Pos_Codon + width),
           Rel_Pos   = seq(- width, width + 1L ) )
  }
}



# Gives you four tibble tables, one for each of the Pos in AAG_df, where width
# has been pre-set. However, it seems as though it is repeating the loop for
# the first position of AAG_df and not moving on to the next row. 

width <- 5L

for (Pos in AAG_df){
  Position <- AAG_df$Pos
  ExpandList <- dplyr::slice(MergedTable, (Position-width):(Position+width), each = FALSE)
  print(ExpandList)
}




# Assign reads to A-sites:
GetGeneReadFrame(gene, dataset, hd_file, left, right, min_read_length,
                 asite_displacement_length = data.frame(
                   read_length = c(28, 29, 30),
                   asite_displacement = c(15, 15, 15)
                 )) 
or 

CalculateGeneReadFrames(dataset, hd_file, gff_df, min_read_length, asite_displacement_length_from_file)




