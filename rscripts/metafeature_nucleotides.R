# Uses the fasta file, gff and h5 file to produce a metafeature window averageing over positions of interest.
# The positions to be averaged across should be provided by the user as a TSV file with the headings Gene and Pos
# The TSV will provide info on the gene and nt position relative to the start codon for each of the features being averaged across

source(here::here("rscripts", "read_count_functions.R"))
source(here::here("rscripts", "stats_figs_block_functions.R"))

library(Biostrings)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(parallel)
library(rhdf5)
library(dplyr)


gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "tiny_2genes_20utrs.gff3")
fasta <- '../example-datasets/simulated/mok/annotation/tiny_2genes_20utrs.fa'
hd_file <- here::here("Mok-tinysim", "output", "A", "A.h5")
dataset <- 'Mok-tinysim'
min_read_length <- 10
expand_width <- 2
features_to_study <- read.delim('nt_testing.tsv')
output_dir <- '.'

gff_df <- readGFFAsDf(gff)
genome <- readDNAStringSet(fasta,format = "fasta")
names <- names(genome)

# CreateNtAnnotation produces a tibble containing the gene name, the CDS nucleotide position and the Nucleotide
# to be used to create meta feature plots based on nucleotide position rather than codon positon
# also allows use on species where a codon table tsv file (ie yeast_codon_table.tsv) is not available, as uses the fasta file

CreateNtAnnotation <- function(genome, names, gff_df){
  gene_seq_df <- data.frame(genome)
  
  # > head(gene_seq_df)
  # genome
  # MAT     AAAAAGAAAACAAAATAAAAATGGCCACATGATTTTTGTTTTCTTTTATTTT
  # MIKE ATAAAGTAAACTAAATTAAAATGATCAAGGAGTAATATTTGATTTCATTTAATTT
  
  GeneSequenceTibble <- tibble(Gene = names, sequence = gene_seq_df$genome)
  
  # head(GeneSequenceTibble)
  # # A tibble: 2 x 2
  # Gene  sequence                                               
  # <chr> <chr>                                                  
  #   1 MAT   AAAAAGAAAACAAAATAAAAATGGCCACATGATTTTTGTTTTCTTTTATTTT   
  # 2 MIKE  ATAAAGTAAACTAAATTAAAATGATCAAGGAGTAATATTTGATTTCATTTAATTT
  
  ConvertSequenceToNt <- function(.x, GeneSequenceTibble){
    gene <- .x
    
    subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene)
    
    left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start))
    # ie for gene 'MAT' left = 21
    
    right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end))
    # ie for gene 'MAT' right = 32
    
    print(gene)
    gene_seq <- dplyr::filter(GeneSequenceTibble, GeneSequenceTibble$Gene == gene)
    
    # > gene_seq
    # # A tibble: 1 x 2
    # Gene  sequence                                            
    # <chr> <chr>                                               
    # 1 MAT   AAAAAGAAAACAAAATAAAAATGGCCACATGATTTTTGTTTTCTTTTATTTT
    
    seq <- unlist(strsplit(gene_seq$sequence, ''))[left:right]
    # For gene 'MAT'
    # seq 
    # [1] "A" "T" "G" "G" "C" "C" "A" "C" "A" "T" "G" "A"
    
    nt_tibble <- tibble(Gene = gene, Pos = 1:length(seq), Nucleotide = seq)
    
  }
  
  nt_tibble <- purrr::map_df(.x = names, .f = ConvertSequenceToNt, GeneSequenceTibble)
  
  # >  head(nt_tibble)
  # # A tibble: 6 x 3
  # Gene    Pos Nucleotide
  # <chr> <int> <chr>     
  #   1 MAT       1 A         
  # 2 MAT       2 T         
  # 3 MAT       3 G         
  # 4 MAT       4 G         
  # 5 MAT       5 C         
  # 6 MAT       6 C 
  
}

# nt_tibble <- CreateNtAnnotation(genome, names, gff_df)

# nt_tibble
# # A tibble: 27 x 3
# Gene    Pos Nucleotide
# <chr> <int> <chr>     
#   1 MAT       1 A         
# 2 MAT       2 T         
# 3 MAT       3 G         
# 4 MAT       4 G         
# 5 MAT       5 C         
# 6 MAT       6 C         
# 7 MAT       7 A         
# 8 MAT       8 C         
# 9 MAT       9 A         
# 10 MAT      10 T         
# # ... with 17 more rows


# Get the location of reads on the transcript for an individual gene

GetReadPositions <- function(gene, dataset, hd_file, asite_displacement_length ,reads_asitepos, left, right){
  
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
  
}

# get the positions of reads for multiple genes and produce a tibble listing the gene, position and count. 

GetAllPosCounts <- function(gene_names, dataset, hd_file, min_read_length){
  
  gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
  
  GetAllPosCounts1Gene <- function(gene, dataset, hd_file, min_read_length, asite_displacement_length){
    
    subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
    
    left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start))
    
    right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end))
    
    asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", 
                                                                                "yeast_standard_asite_disp_length.txt"))
    
      
    nt_counts_1_gene <- GetReadPositions(gene,dataset, hd_file, asite_displacement_length, left, right)

    nt_pos_counts <- tibble(Gene = gene,
                               Pos = 1:length(nt_counts_1_gene),
                               Count = nt_counts_1_gene)
    
    as.data.frame(nt_pos_counts, row.names = NULL, optional = FALSE)
    
    return(nt_pos_counts)
    
  }
  
  total_nt_pos_counts <- purrr::map_dfr(.x = gene_names,
                                           .f = GetAllPosCounts1Gene,
                                           dataset,
                                           hd_file,
                                           min_read_length,
  )
  
  return (total_nt_pos_counts)
}

total_nt_pos_counts <- suppressMessages(GetAllPosCounts(gene_names, dataset, hd_file, min_read_length))

# total_nt_pos_counts
# # A tibble: 27 x 3
# Gene    Pos Count
# <chr> <int> <dbl>
#   1 MAT       1     0
# 2 MAT       2     0
# 3 MAT       3     0
# 4 MAT       4     0
# 5 MAT       5     1
# 6 MAT       6     1
# 7 MAT       7     2
# 8 MAT       8     0
# 9 MAT       9     0
# 10 MAT      10     0


# add the nucleotide identity to the total_nt_pos_counts tibble

AddNtToPosCounts <- function(gene_names, dataset, hd_file, min_read_length, gff_df){
  
  nt_tibble <- CreateNtAnnotation(genome, names, gff_df)
  
  total_nt_pos_counts <- GetAllPosCounts(gene_names, dataset, hd_file, min_read_length)
  
  transcript_tibbles <- left_join(total_nt_pos_counts, nt_tibble, by = c("Pos", "Gene"), keep = FALSE, copy = TRUE)
  
  transcript_gene_pos_nt_reads <- tibble(
    Gene = transcript_tibbles$Gene,
    Pos = transcript_tibbles$Pos,
    Count = transcript_tibbles$Count,
    Nucleotide = transcript_tibbles$Nucleotide
  )
  
  return(transcript_gene_pos_nt_reads)
}   


transcript_gene_pos_nt_reads <- suppressMessages(AddNtToPosCounts(gene_names, dataset, hd_file, min_read_length, gff_df))

# the output is a tibble with all nucleotide positions in the cds of a transcript, along with the number of reads mapping to each position

# transcript_gene_pos_nt_reads
# # A tibble: 27 x 4
# Gene    Pos Count Nucleotide
# <chr> <int> <dbl> <chr>     
#   1 MAT       1     0 A         
# 2 MAT       2     0 T         
# 3 MAT       3     0 G         
# 4 MAT       4     0 G         
# 5 MAT       5     1 C         
# 6 MAT       6     1 C         
# 7 MAT       7     2 A         
# 8 MAT       8     0 C         
# 9 MAT       9     0 A         
# 10 MAT      10     0 T    


# give a list of genes and positions to average over
# for each item in that list, take a window of a specified size around 




###  Slice out window around interesting features ###

# ExpandFeatureRegionAllGenes takes a slice out of transcript_gene_pos_poscodon_frame, centered on the position of the feature of interest
# It then sets the position of the feature of interest to 0, and changes the positions of adjacent codons to be relative to the feature of interest
# This is done for all occurrences of the feature of interest on all genes in the sample provided.



ExpandFeatureRegionAllGenes <- function(gene_names, dataset, hd_file, 
                                        min_read_length,   
                                        gff_df, features_to_study,
                                        expand_width) {
  
  # generate transcript_gene_pos_poscodon_frame, tests and descriptions above
  
  transcript_gene_pos_nt_reads <- suppressMessages(AddNtToPosCounts(gene_names, dataset, hd_file, min_read_length, gff_df))
  
  # take as inputs transcript_gene_pos_poscodon_frame and select for positions on separate genes
  
    # Use transcript one gen as an input to ExpandRegions
    # For each occurance of the codon of interest on the gene being studied, slice out a window of positions around the feature of interest from transcript_gene_pos_poscodon_frame 
    
    # return an empty tibble if the desired region hangs over the edge of the coding region
    
    ExpandRegions <- function(rows, features_to_study, transcript_gene_pos_nt_reads, gff_df, expand_width){
      
      # select the row of features
      
      features <- features_to_study[rows,]
      print(features)
      
      # filter to the relevant gene 
      
      transcript_gene_pos_nt_gene_interest <- dplyr::filter(transcript_gene_pos_nt_reads, Gene == features$Gene)
      
      gene_length <- filter(gff_df, gff_df$type == 'CDS' & gff_df$Name == features$Gene)$width
      
      # slice out a window based on the expand width, centered on the position listed in features, from transcript_gene_pos_nt_gene_interest
      
      if (features$Pos <= expand_width  |features$Pos + expand_width > gene_length) {
        return()
      } else {
        output_feature_info <- tibble(
          dplyr::slice(transcript_gene_pos_nt_gene_interest, (features$Pos - expand_width):(features$Pos + expand_width), each = FALSE),
          Rel_Pos =  seq(- expand_width, expand_width)
        )
        
        if(dim(output_feature_info)[1] == (2*expand_width + 1)){
          
          return(output_feature_info)
        }else{
          return()
        }
      }
    }
    # The if statement ensures that feature positions that are less/more than the expand_width value are discarded
    
    # get the number of rows being studied and create a vector of 1 to n rows.
    
    rows<- 1:nrow(features_to_study)
    
    # itterate ExpandRegions over each row in features_to_study
    
    output_feature_info <- purrr::map(.x = rows, .f = ExpandRegions, features_to_study,
                                      transcript_gene_pos_nt_reads,gff_df,
                                      expand_width)
  
  # produces a list containing an expanded window for each location listed in features_to_study
  
  # remove NULLS, which represent features of interest occuring within one expand width of the UTRs
  output_feature_info <- output_feature_info[!sapply(output_feature_info, is.null)]
  
  return(output_feature_info)
}

output_feature_info <- ExpandFeatureRegionAllGenes(gene_names, dataset, hd_file, 
                                                   min_read_length,   
                                                   gff_df, features_to_study,
                                                   expand_width)


# TO DO:: what happens if one of the features listed doesn't exist?


#TEST: ExpandFeatureRegionAllGenes(): creates an object of type "list" 
#TEST: ExpandFeatureRegionAllGenes(): length(list) == nrow(features_to_study)
#TEST: ExpandFeatureRegionAllGenes(): each tibble contains 5 columns = TRUE
#TEST: ExpandFeatureRegionAllGenes(): the column names are %in% c("Gene", "Pos", "Count", "Nucleotide", "Rel_Pos") 
#TEST: ExpandFeatureRegionAllGenes(): number of observations in the output tibble = "expand_width" * 2 + 1, so if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandFeatureRegionAllGenes(): the position from "interesting_feature_positions" has "Rel_Pos" value 0 = TRUE
#TEST: ExpandFeatureRegionAllGenes(): the column "Rel_Pos" goes from -"expand_width to +"expand_width"

## example, when features_to_study is

# Gene	Pos
# MIKE	5
# MAT	7
# MAT	5
# MAT	6

# produces:

# output_feature_info
# [[1]]
# # A tibble: 5 x 5
# Gene    Pos Count Nucleotide Rel_Pos
# <chr> <int> <dbl> <chr>        <int>
#   1 MIKE      3     0 G               -2
# 2 MIKE      4     0 A               -1
# 3 MIKE      5     1 T                0
# 4 MIKE      6     0 C                1
# 5 MIKE      7     0 A                2
# 
# [[2]]
# # A tibble: 5 x 5
# Gene    Pos Count Nucleotide Rel_Pos
# <chr> <int> <dbl> <chr>        <int>
#   1 MAT       5     1 C               -2
# 2 MAT       6     1 C               -1
# 3 MAT       7     2 A                0
# 4 MAT       8     0 C                1
# 5 MAT       9     0 A                2
# 
# [[3]]
# # A tibble: 5 x 5
# Gene    Pos Count Nucleotide Rel_Pos
# <chr> <int> <dbl> <chr>        <int>
#   1 MAT       3     0 G               -2
# 2 MAT       4     0 G               -1
# 3 MAT       5     1 C                0
# 4 MAT       6     1 C                1
# 5 MAT       7     2 A                2
# 
# [[4]]
# # A tibble: 5 x 5
# Gene    Pos Count Nucleotide Rel_Pos
# <chr> <int> <dbl> <chr>        <int>
#   1 MAT       4     0 G               -2
# 2 MAT       5     1 C               -1
# 3 MAT       6     1 C                0
# 4 MAT       7     2 A                1
# 5 MAT       8     0 C                2




### Normalization ###

# Normalization carried out within each expanded frame so that they are comparable 
# Normalizes the expand_feature_region list generating a RelCount column with the normalization values

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

normalized_expand_list <- purrr::map(
  .x = output_feature_info,
  .f = ExpandedRegionNormalization,
  expand_width
)


#TEST: ExpandedRegionNormalization(): creates a tidy format data frame (tibble) = TRUE
#TEST: ExpandedRegionNormalization(): the tibble contains 6 columns = TRUE
#TEST: ExpandedRegionNormalization(): the column names are %in% c("Gene", "Pos", "Count", "Nucleotide", "Rel_Pos", "RelCount")
#TEST: ExpandedRegionNormalization(): number of observations in the output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: ExpandedRegionNormalization(): the column "Rel_Pos" goes from -"expand_width to +"expand_width" 
#TEST: ExpandedRegionNormalization(): sum(normalized_expand_list[[1]]$RelCount)/nrow(normalized_expand_list[[1]]) == 1 
#TEST: ExpandedRegionNormalizetion(): None of the RelCount columns should contain NaN. These should all be set to 0

# TEST:: normalised_expand_list should be of type list
# type(normalized_expand_list)
# [1] "list"

# TEST:: Normalized_expand_list should be the same length as expand_feature_region
# > length(normalized_expand_list)==length(expand_feature_region)
# [1] TRUE

# TEST:: the dimensions of each item in the list shoud be [(2*expand_width+1) X 6] as there are now 5 rows; Gene, Pos, Count, Nucleotide, Rel_Pos, RelCount
# > dim(expand_feature_region[[1]])
# [1] 11  5

# TEST:: At each position, the sum of RelCount should be equal to (2*expand_width+1)
# ie if the expand width was 5:
# sum(normalized_expand_list[[1]]$RelCount)
# [1] 11

# Example: when features_to_study is
# 
# Gene	Pos
# MIKE	5
# MAT	7
# MAT	5
# MAT	6

# The output is

# > normalized_expand_list
# [[1]]
# # A tibble: 5 x 6
# Gene    Pos Count Nucleotide Rel_Pos RelCount
# <chr> <int> <dbl> <chr>        <int>    <dbl>
#   1 MIKE      3     0 G               -2        0
# 2 MIKE      4     0 A               -1        0
# 3 MIKE      5     1 T                0        5
# 4 MIKE      6     0 C                1        0
# 5 MIKE      7     0 A                2        0
# 
# [[2]]
# # A tibble: 5 x 6
# Gene    Pos Count Nucleotide Rel_Pos RelCount
# <chr> <int> <dbl> <chr>        <int>    <dbl>
#   1 MAT       5     1 C               -2     1.25
# 2 MAT       6     1 C               -1     1.25
# 3 MAT       7     2 A                0     2.5 
# 4 MAT       8     0 C                1     0   
# 5 MAT       9     0 A                2     0   
# 
# [[3]]
# # A tibble: 5 x 6
# Gene    Pos Count Nucleotide Rel_Pos RelCount
# <chr> <int> <dbl> <chr>        <int>    <dbl>
#   1 MAT       3     0 G               -2     0   
# 2 MAT       4     0 G               -1     0   
# 3 MAT       5     1 C                0     1.25
# 4 MAT       6     1 C                1     1.25
# 5 MAT       7     2 A                2     2.5 
# 
# [[4]]
# # A tibble: 5 x 6
# Gene    Pos Count Nucleotide Rel_Pos RelCount
# <chr> <int> <dbl> <chr>        <int>    <dbl>
#   1 MAT       4     0 G               -2     0   
# 2 MAT       5     1 C               -1     1.25
# 3 MAT       6     1 C                0     1.25
# 4 MAT       7     2 A                1     2.5 
# 5 MAT       8     0 C                2     0 




### Overlaying the normalized expanded tibbles ###

# Function to overlay graphs into a single graph. Need to generate a single tibble 
# from NormalizedExpandList. Join by Rel_Pos, in RelCount need the mean for 
# each Rel_Pos (sum row(x) / number of row(x))

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

overlayed_tibbles <- OverlayedTable(normalized_expand_list, expand_width) 

#TEST: OverlayedTable(): creates a tidy format data frame (tibble) = TRUE
#TEST: OverlayedTable(): the tibble contains 2 columns = TRUE
#TEST: OverlayedTable(): the column names are %in% c("Rel_Pos", "RelCount")
#TEST: OverlayedTable(): number of observations in the output tibble = "expand_width"*2+1, if "expand_width" = 5L the number of observations should be 11
#TEST: OverlayedTable(): the column "Rel_Pos" goes from -"expand_width to +"expand_width" 
#TEST: OverlayedTable(): RelCount is a numeric 

# Example: when features_to_study is
# 
# Gene	Pos
# MIKE	5
# MAT	7
# MAT	5
# MAT	6

# overlayed_tibbles
# # A tibble: 5 x 2
# Rel_Pos RelCount
# <int>    <dbl>
#   1      -2    0.312
# 2      -1    0.625
# 3       0    2.5  
# 4       1    0.938
# 5       2    0.625

overlayed_plot <- ggplot(overlayed_tibbles, mapping = aes(x = Rel_Pos, y = RelCount)) + 
  geom_line() +
  theme_bw() +
  theme(text=element_text(size=14),
        axis.title=element_text(size=14, face='bold'),
        title = element_text(size = 14, face='bold'))+
  labs(title = paste0('Relative read counts arof positions of interest'),
       x = 'Position relative to feature of interest',
       y = 'Relative read count', size = 2) + 
  scale_x_continuous(breaks = seq(-expand_width, expand_width, 2))

# save plot as PDF

print('Save plot as PDF')
save_plot_pdf <- function(overlayed_plot, output_dir){
  overlayed_plot %>%
    ggsave(
      filename = file.path(output_dir, paste0("Meta_feature_plot_positons_of_interest", '_', dataset,".pdf")),
      width = 6, height = 5
    )
}

save_plot_pdf(overlayed_plot, output_dir)
