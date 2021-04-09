
# source packages and functions from rscripts 
source(here::here("rscripts", "read_count_functions.R"))
source(here::here("rscripts", "stats_figs_block_functions.R"))
source(here::here("rscripts", "provenance.R"))


# other packages I loaded that I used 
library(tidyverse)
library(data.table) #used to load .tsv file 


# Import the .tsv file: (used the updated file which contains the first 
# 200 codons of each gene)
yeast_codon_pos_i200.tsv <- data.table::fread('https://github.com/riboviz/riboviz/raw/change-yeast_codon_pos_i200-194/data/yeast_codon_table.tsv')

# Confirm that the .tsv file has been imported and check structure:
head(yeast_codon_pos_i200.tsv) 



# Filter down the yeast_codon_pos_i200.tsv file to the gene that you are 
# working with, in this case YAL003W

YAL003W_pos <- dplyr::filter(yeast_codon_pos_i200.tsv, Gene=="YAL003W")
print(YAL003W_pos)




# read the structure of the h5 file 
h5ls("~/R/YAL5.h5")

# pathway to the h5 file 
YAL5 <- "~/R/YAL5.h5"



### NOT DIRECTLY USED YET ###

# GetDatamatrix: Function to get the datamatrix for each of the genes 
# contained within the h5 file

# Loop list of gene names through GetGeneDatamatrix to extract all the 
# information

# Required parameters: gene_names, dataset, hd_file)

gene_names <- rhdf5::h5ls(YAL5, recursive = 1)$name

GetDatamatrix <- function(dataset, hd_file, gene_names){
  for (gene in gene_names){
    GetDatamatrixList <- GetGeneDatamatrix(gene, dataset, hd_file)
    print(GetDatamatrixList)
  }
}

# Result: GetDatamatrixList which contains the datamatrix for all of the genes 
# contained within the list 


# Extract the datamatrix for the five genes 
Mok_YAL5_data <- GetDatamatrix(dataset="Mok-simYAL5", hd_file = "~/R/YAL5.h5", 
                               gene_names)

### ###




# GetGeneDatamatrix for each of the genes contained within the Mok-simYAL5

YAL003W_data <- GetGeneDatamatrix(gene="YAL003W", dataset="Mok-simYAL5", 
                                  hd_file = "~/R/YAL5.h5") #data_mat

# The GetGeneDatamatrix contains the CDS as well as the UTRs, results in 1121 
# observations instead of the 207 which consists of the CDS

# YAL005C_data <- GetGeneDatamatrix(gene="YAL005C", dataset="Mok-simYAL5", 
#                                   hd_file = "~/R/YAL5.h5")
# YAL012W_data <- GetGeneDatamatrix(gene="YAL012W", dataset="Mok-simYAL5", 
#                                   hd_file = "~/R/YAL5.h5")
# YAL035W_data <- GetGeneDatamatrix(gene="YAL035W", dataset="Mok-simYAL5", 
#                                   hd_file = "~/R/YAL5.h5")
# YAL038W_data <- GetGeneDatamatrix(gene="YAL038W", dataset="Mok-simYAL5", 
#                                   hd_file = "~/R/YAL5.h5")



# Make a table which contains the datamatrix from the h5 file, in this case only 
# carried out for a single gene YAL003W



GetGeneStartCodonPos(gene="YAL003W", dataset="Mok-simYAL5", # 251 252 253
                     hd_file = "~/R/YAL5.h5") 

GetGeneStopCodonPos(gene="YAL003W", dataset="Mok-simYAL5",  # 869, 870, 871
                    hd_file="~/R/YAL5.h5")

# CDS start 251, CDS end 871 - CDS total: 621 nt -> 207 codons
# 250 UTRs at each end for a total of 1121 nt


gff_df<-readGFFAsDf("~/R/Scer_YAL_5genes_w_250utrs.gff3")

                                        
  #   reads_pos_length <- GetGeneDatamatrix(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = "~/R/YAL5.h5")
  #           print(reads_pos_length)
  #   reads_asitepos <- CalcAsiteFixed(reads_pos_length, min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE) 
  #           print(reads_asitepos)
  #   SumByFrame(reads_asitepos, left =251, right =871)
  #   SnapToCodonList <- SnapToCodon(reads_asitepos, left=251, right=871, snapdisp = 0L)                                         
  #                                        
  # print(SnapToCodonList)                                       

# buffer <- 250
# min_read_length <- 10
# Codon_Pos_Reads <- GetCodonPositionReads(gene = "YAL003W", dataset = "Mok-simYAL5", hd_file = "~/R/YAL5.h5", left = (buffer - 15), right = (buffer + 11), min_read_length = min_read_length)


YAL003W_tidy <- TidyDatamatrix(data_mat = YAL003W_datamatrix, startpos = 1, startlen=10)

# Summaried the table so that the total count is shown for each of the positions 
Tidy_sum  <- summarise(group_by(YAL003W_tidy, Pos), sum(Counts))
Tidy_remove_UTR5 <- Tidy_sum [-c(1:250), ]
Tidy_remove_UTR3 <- Tidy_remove_UTR5[c(1:621), ]
#Tidy_nnt_to_codons <- 
  


# Tidy_sum still contains 1121 nucleotides, so still contains the UTRs as well
# Format is still then in the nucleotide format so need to divide by three to get codons 
# The YAL003W_pos table has the positions as 1,2,3... while Tidy is from 251





# Make final table which contains Gene, Pos, Codon, Counts  
CodonTable=c()
CodonTable$Gene<-YAL003W_pos$Gene
CodonTable$Pos<-YAL003W_pos$PosCodon
CodonTable$Codon<-YAL003W_pos$Codon
Codon_table <- as_tibble(CodonTable)
print(Codon_table)

Merged_table <- inner_join(Codon_table, Tidy_sum ) # Table which contains the columns Gene, Pos, Codon, Counts
print(Merged_table)

# the column of counts has the wrong name 



# Filter down to the desired codon 
AAG_table <- dplyr::filter(MergedTable, Codon=="AAG")
AAG_df<-as.data.frame(AAG_table)


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




