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
