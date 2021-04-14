suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))
suppressMessages(library(rhdf5))
suppressMessages(library(parallel))
suppressMessages(library(optparse))
suppressMessages(library(RcppRoll))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(purrr))
suppressMessages(library(here))

## functions to be used in code - replace with source(here::here("rscripts", "read_count_functions.R"))

#Get a data matrix for the gene of interest 
GetGeneDatamatrix <- function(gene, dataset, hd_file){
  rhdf5::h5read(file = hd_file, name = paste0("/", gene, "/", dataset, "/reads/data")) %>%
    return()
}

# create tidy data matrix
TidyDatamatrix <- function(data_mat, startpos = 1, startlen = 1) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(data_mat) - 1)
  readlengths <- startlen:(startlen + nrow(data_mat) - 1)
  data_mat %>%
    set_colnames(positions) %>%
    as_tibble() %>%
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts))
}

## Actual code
# path to H5 file
file_url <- "wt.noAT.ribo.4_s/4c7158e339494e80dadcfd1a859f62/wt.noAT.ribo.4_s.h5"
Gene_of_interest <- 'SPCC188.09c.1'
Dataset <- 'D-Sp_2018'

# Create the gene data matrix 
gene_data_matrix <- GetGeneDatamatrix(gene= Gene_of_interest,
                                      dataset = Dataset,
                                      hd_file = file_url)

# extract the buffer so can identify stop codon later
buffer_left <- h5readAttributes(file_url, base::paste('/',Gene_of_interest,'/D-Sp_2018/reads', sep = ''))[['buffer_left']]
buffer_right <- h5readAttributes(file_url, base::paste('/',Gene_of_interest,'/D-Sp_2018/reads', sep = ''))[['buffer_right']]
nnt_buffer <- buffer_left

# Create a Tidy data matrix, using the the gene data matrix. set start as -nnt_buffer + 1
# so the actual start codon lies on 1
gene_of_interest_tidy_matrix <- TidyDatamatrix(data_mat = gene_data_matrix, startpos = -nnt_buffer +1 )


# make tidy data matrix, with all the counts of different read lengths at different positions
# aim: make a tidyDataMatrix with two columns: position and count
# input: Fil1TidyDataMatrix 
# process: for position i, take the sum of the counts at that position across read lengths
#          and return a new row for the new matrix with the position and the counts 
# for each position, take all the counts of reads of different sizes, and combine them into one 
# 'total read' column for each position in the gene
# output: TidyDataMatrix

# create an empty tibble with two columns; Pos and Counts.
gene_Total_reads_at_position <- tibble(Pos =integer(),
                                            Counts = integer())

# Loop through tidy data matrix, and for each position add up all of the counts for different 
# read lengths, storing them in the new tidy data matrix
for(i in gene_of_interest_tidy_matrix[1,]$Pos:max(gene_of_interest_tidy_matrix$Pos)){
  tmp_row <- gene_of_interest_tidy_matrix %>% filter(gene_of_interest_tidy_matrix$Pos ==i)
  new_row <- tibble(Pos = i, Counts = sum(tmp_row$Counts))
  gene_Total_reads_at_position <- gene_Total_reads_at_position %>% bind_rows(new_row)
}

# plot the data, so positions are along the x axis and number of counts is along the y axis
ggplot(gene_Total_reads_at_position,aes(x = Pos, y = Counts))+
  geom_col( width = 1, color = 'red')+
  labs(title = paste(file_url, ' - ', Gene_of_interest), x = 'Position relative to start codon', y = 'Number of reads')
