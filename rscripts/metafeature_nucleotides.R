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
gff_df <- readGFFAsDf(gff)
fasta <- '../example-datasets/simulated/mok/annotation/tiny_2genes_20utrs.fa'
hd_file <- here::here("Mok-tinysim", "output", "A", "A.h5")
dataset <- 'Mok-tinysim'
min_read_length <- 10


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

nt_tibble <- CreateNtAnnotation(genome, names, gff_df)

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

AddNtToPosCounts <- function(nt_tibble, gene_names, dataset, hd_file, min_read_length,   gff_df){
  
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


transcript_gene_pos_nt_reads <- suppressMessages(AddNtToPosCounts(nt_tibble, gene_names, dataset, hd_file, min_read_length, gff_df))

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
