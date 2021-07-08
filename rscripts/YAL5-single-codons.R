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

# Given an h5 file, GFF file and .tsv file, this script creates a metafeature plot 
# for the codon of interest. 

# FOR SINGLE CODONS: USE yeast_codon_pos_i200 THROUGHOUT


# Load the necessary input files 

YAL5_h5 <- here::here("Mok-simYAL5", "output", "A", "A.h5")
# The h5 file from the dataset of interest 

YAL5_gff <- here::here("..", "example-datasets", "simulated", "mok", "annotation", "Scer_YAL_5genes_w_250utrs.gff3")
gff_df <- readGFFAsDf(YAL5_gff)
# The GFF file for the simulated dataset, given the general name gff_df so that 
# the script does not have to change based on different gff_df files being used 

# Import the .tsv file: 
yeast_codon_pos_i200 <- readr::read_tsv(file = here::here("data", "yeast_codon_table.tsv"))

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


# Filter down the yeast_codon_pos_i200 file to the gene that you are 
# working with, in this case YAL003W
YAL003W_pos <- dplyr::filter(yeast_codon_pos_i200, Gene=="YAL003W")

# structure of the h5 file 
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


# Fetch the datamatrix for a single gene of interest , e.g YAL003W
reads_pos_length <- GetGeneDatamatrix(gene="YAL003W", 
                                      dataset="Mok-simYAL5", 
                                      hd_file=YAL5_h5)
  
  # > str(reads_pos_length)
  # int [1:41, 1:1121] 0 0 0 0 0 0 0 0 0 0 ...


# Fetch the tidydatamatrix for a single gene of interest, e.g. YAL003W 
tidy_gene_datamatrix <- TidyDatamatrix(GetGeneDatamatrix(gene="YAL003W", 
                                                         dataset="Mok-simYAL5", 
                                                         hd_file=YAL5_h5), 
                                       startpos = 1, 
                                       startlen = 10) 



### CalcAsiteFixed (Aim: A-site assignment) 

# CalcAsiteFixed for single gene of interest 
asite_counts_by_position <- CalcAsiteFixed(reads_pos_length, 
                                           min_read_length = 10, 
                                           asite_displacement_length = data.frame(read_length = c(28, 29, 30), 
                                                                                  asite_displacement = c(15, 15, 15)), 
                                           colsum_out = TRUE)

  # > str(asite_counts_by_position)
  # num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...

  # This step extracts the asite reads from the reads_pos_length file
  # This step is incorporated in the funciton below 


# Funtction with the aim of combining the asite count with the position of the nucleotides           
TidyAsiteCountsByPosition <- function(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE){
  
  asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(here::here("data", 
                                                                              "yeast_standard_asite_disp_length.txt"))
  
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  
  asite_counts_by_position <- CalcAsiteFixed(reads_pos_length, 
                                             min_read_length = min_read_length, 
                                             asite_displacement_length = asite_displacement_length,
                                             colsum_out = colsum_out)
  
  tidy_asite_counts <- tibble(Pos = 1:length(asite_counts_by_position), 
                              Count = asite_counts_by_position)
  
  return(tidy_asite_counts)
  
}

tidy_asite_count_output <- TidyAsiteCountsByPosition(gene = "YAL003W", 
                                                     dataset = "Mok-simYAL5", 
                                                     hd_file = YAL5_h5, 
                                                     min_read_length = 10, 
                                                     colsum_out = TRUE)

  # > str(tidy_asite_count_output)
  # tibble [1,121 x 2] (S3: tbl_df/tbl/data.frame)
  # $ Pos  : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
  # $ Count: num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...

  # The end result here is that the asite counts are aligned to the gene of interest
  # (in nucleotides, including UTRs and CDS)



# TranscriptPosToCodonPos generates a table with the columns Gene, Pos (position 
# of the nucleotide), Pos_Codon and frame (reading frame).

TranscriptPosToCodonPos <- function(gene, gff_df){
  
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene) 
  # the gff file is filtered down to the gene of interest
  
  UTR5_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR5") %>% select(width)
  # the width of the UTR5 is defined from the gff file, in this case 250 
  
  CDS_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(width)
  # the width of the CDS is defined from the gff file, in this case 621
  
  UTR3_width <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(width) 
  # the width of the UTR3 is defined from the gff file, in this case 250 
  
  transcript_length <- dplyr::filter(.data = subset_gff_df_by_gene, type == "UTR3") %>%  select(end)
  # defines the length of the transcript including both UTRs and the CDS (in nucleotides), in this case 1121
  
  CDS_start <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(start)
  # define the start position of the CDS from the gff file, in this case 251
  
  CDS_end <- dplyr::filter(.data = subset_gff_df_by_gene, type == "CDS") %>%  select(end)
  # define the end position of the CDS from the gff file, in this case 871
  
  CDS_codon_positions <- rep(1:(CDS_width$width/3), each = 3)
  # In order to add the CDS codon positions the length of the CDS in nucleotides is 
  # divided by 3 (3 nucleotides = 1 codon) and repeated 3 times
  
    # > str(CDS_codon_positions)
    # int [1:621] 1 1 1 2 2 2 3 3 3 4 ...
  
  NA_UTR5_width <- rep(NA, times = UTR5_width$width)
  # NAs are repeated for the width of UTR5 
  
  NA_UTR3_width <- rep(NA, times = UTR3_width$width)
  # NAs are repeated for the width of UTR3
  
  transcript_gene_pos_poscodon_frame <- tibble(
    Gene = gene,
    Pos = 1:transcript_length$end,
    Pos_Codon = c(rep(NA, times = UTR5_width$width), 
                  CDS_codon_positions, 
                  rep(NA, times = UTR3_width$width)),
    Frame = seq(from = 2L, to = (transcript_length$end + 1L)) %% 3 # works as UTRS are 250, might not work with UTRS of others values
    # add the count column, likely something we want here 
  )
  
  return(transcript_gene_pos_poscodon_frame)
}

transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene = "YAL003W", gff_df)

  # The end result is a table with the columns Gene, Pos, Pos_Codon and Frame (reading frame).
  # For the Pos_Codon column the UTR positions have NA, while the CDS has the codon positions 

  # > str(transcript_pos_to_codon_pos_output)
  # tibble [1,121 x 4] (S3: tbl_df/tbl/data.frame)
  # $ Gene     : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # $ Pos      : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
  # $ Pos_Codon: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
  # $ Frame    : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...



# Function joins the tables which are output by TidyAsiteCountsByPosition and 
# TranscriptPosToCodonPos functions to generate transcript_info_tibble

AddAsiteCountsToTranscriptPosToCodonPos <- function(gene, dataset, hd_file, min_read_length = 10, colsum_out = TRUE, gff_df){
  
  tidy_asite_count_output <- TidyAsiteCountsByPosition(gene, 
                                                       dataset, 
                                                       hd_file, 
                                                       min_read_length = 10, 
                                                       colsum_out = TRUE)
  
  transcript_pos_to_codon_pos_output <- TranscriptPosToCodonPos(gene, 
                                                                gff_df)
  
  transcript_info_tibble <- dplyr::left_join(transcript_pos_to_codon_pos_output, 
                                             tidy_asite_count_output, 
                                             by = "Pos")
  
  return(transcript_info_tibble)
  
}

transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene = "YAL003W", 
                                                                  dataset = "Mok-simYAL5", 
                                                                  hd_file = YAL5_h5, 
                                                                  min_read_length = 10, 
                                                                  colsum_out = TRUE, 
                                                                  gff_df)

  # > str(transcript_info_tibble)
  # tibble [1,121 x 5] (S3: tbl_df/tbl/data.frame)
  # $ Gene     : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # $ Pos      : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
  # $ Pos_Codon: int [1:1121] NA NA NA NA NA NA NA NA NA NA ...
  # $ Frame    : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...
  # $ Count    : num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...

  # End result is a table which contains the asite counts, reading frame, codons and 
  # codon positions for the gene of interest



# Function to add codon names to transcript_info_tibble
AddCodonNamesToTranscriptInfoTibble <- function(yeast_codon_pos_i200, 
                                                gene, dataset, 
                                                hd_file, 
                                                min_read_length = 10, 
                                                colsum_out = TRUE, 
                                                gff_df){
  
  codon_table <- dplyr::filter(yeast_codon_pos_i200, Gene == gene) %>%
    dplyr::select("PosCodon", "Codon")
  # this provides the positions of the codons and the codon names which can 
  # joined to the tibble generated by AddAsiteCountsToTranscriptPosToCodonPos 
  
  transcript_info_tibble <- AddAsiteCountsToTranscriptPosToCodonPos(gene, 
                                                                    dataset, 
                                                                    hd_file, 
                                                                    min_read_length = 10, 
                                                                    colsum_out = TRUE, 
                                                                    gff_df)

  transcript_info_tibble <- left_join(transcript_info_tibble, 
                                      codon_table, 
                                      by = c("Pos_Codon" = "PosCodon"), 
                                      keep = FALSE)
  
  return(transcript_info_tibble)
  
}

transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(yeast_codon_pos_i200, 
                                                              gene = "YAL003W", 
                                                              dataset = "Mok-simYAL5", 
                                                              hd_file = YAL5_h5, 
                                                              min_read_length = 10, 
                                                              colsum_out = TRUE, 
                                                              gff_df)

  # > str(transcript_info_tibble)
  # tibble [1,121 x 6] (S3: tbl_df/tbl/data.frame)
  # $ Gene     : chr [1:1121] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # $ Pos      : int [1:1121] 1 2 3 4 5 6 7 8 9 10 ...
  # $ Pos_Codon: num [1:1121] NA NA NA NA NA NA NA NA NA NA ...
  # $ Frame    : num [1:1121] 2 0 1 2 0 1 2 0 1 2 ...
  # $ Count    : num [1:1121] 0 0 0 0 0 0 0 0 0 0 ...
  # $ Codon    : chr [1:1121] NA NA NA NA ...



### Filter for reading frame of interest (default = 0)
FilterForFrame <- function(transcript_info_tibble, filtering_frame = 0){
  
  transcript_info_tibble <- dplyr::filter(transcript_info_tibble, 
                                          Frame == filtering_frame)
  
  return(transcript_info_tibble)
  
}

transcript_info_tibble <- FilterForFrame(transcript_info_tibble, 
                                         filtering_frame = 0)

  # removes reads which are in the reading frame 1 and 2 

  # > str(transcript_info_tibble)
  # tibble [374 x 6] (S3: tbl_df/tbl/data.frame)
  # $ Gene     : chr [1:374] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # $ Pos      : int [1:374] 2 5 8 11 14 17 20 23 26 29 ...
  # $ Pos_Codon: num [1:374] NA NA NA NA NA NA NA NA NA NA ...
  # $ Frame    : num [1:374] 0 0 0 0 0 0 0 0 0 0 ...
  # $ Count    : num [1:374] 0 0 0 0 0 0 0 0 0 0 ...
  # $ Codon    : chr [1:374] NA NA NA NA ...



# Function to filter for codons of interest, generates a list of the positions of the codon of interest
FilterForCodonOfInterestPositions <- function(transcript_info_tibble, codon_of_interest){

  # transcript_info_tibble <- AddCodonNamesToTranscriptInfoTibble(yeast_codon_pos_i200,
  #                                                               gene,
  #                                                               dataset,
  #                                                               hd_file,
  #                                                               min_read_length = 10,
  #                                                               colsum_out = TRUE,
  #                                                               gff_df)

  # generating the transcript_info_tibble inside of the function seems to result
  # in the list that is generated being triplicated. Skipping this step seems to
  # resolve that issue

  interesting_codon_table <- dplyr::filter(transcript_info_tibble, Codon == codon_of_interest)

  interesting_codon_positions <- interesting_codon_table$Pos_Codon

  return(interesting_codon_positions)
}

interesting_codon_positions <- FilterForCodonOfInterestPositions(transcript_info_tibble, codon_of_interest = "TCT")

  # > interesting_codon_positions
  # [1]  18  31  43  64  69 197



### Expand frame around codon of interest 

# Function needed to go from one position to an exapnded region.
# This function uses codon position values, not nucleotide position values
ExpandCodonRegion <- function(.x = interesting_codon_positions, transcript_info_tibble, gene, dataset, hd_file, expand_width = 5L, remove_overhang = TRUE) {

  gene_poscodon_count <- na.omit(tibble(
    Gene = transcript_info_tibble$Gene,
    Pos_Codon = transcript_info_tibble$Pos_Codon,
    Rel_Count = transcript_info_tibble$Count
  ))

  # added na.omit() as the slice function goes by the row number instead of Pos_codon
  # na.omit() removes the UTRs so only CDS remains

  gene_length <- GetGeneLength(gene, dataset, hd_file) # giving me 618??

  if (.x < expand_width | .x + expand_width > gene_length ) {
    return( tibble() )
  } else {
    output_codon_info <- tibble(
      dplyr::slice(gene_poscodon_count, (.x - expand_width):(.x + expand_width), each = FALSE),
      Rel_Pos =  seq(- expand_width, expand_width)
    )
    return(output_codon_info)
  }
  
  # The if statement ensures that codon positions that are less/more than the 
  # expand_width value are discarded 
}

output_codon_info <- ExpandCodonRegion(.x = interesting_codon_positions, 
                                       transcript_info_tibble , 
                                       gene="YAL003W", 
                                       dataset="Mok-simYAL5", 
                                       hd_file=YAL5_h5, 
                                       expand_width = 5L, 
                                       remove_overhang = TRUE)

  # # Gives the output:
  # 
  # > ExpandCodonRegionOutput
  # # A tibble: 11 x 4
  # Gene    Pos_Codon Rel_Count Rel_Pos
  # <chr>       <dbl>     <dbl>   <int>
  #   1 YAL003W        13       490      -5
  # 2 YAL003W        14       561      -4
  # 3 YAL003W        15        51      -3
  # 4 YAL003W        16       567      -2
  # 5 YAL003W        17        67      -1
  # 6 YAL003W        18        65       0
  # 7 YAL003W        19       158       1
  # 8 YAL003W        20        71       2
  # 9 YAL003W        21       323       3
  # 10 YAL003W        22       288       4
  # 11 YAL003W        23      1301       5



# Apply the ExpandCodonRegion function to the codons of interest to generate expanded tibbles for each position 

ExpandCodonRegionForList <- function(transcript_info_tibble, gene, dataset, hd_file, startpos = 1, startlen = 10, codon_of_interest, gff_df, expand_width = 5L){
  
  interesting_codon_positions <- FilterForCodonOfInterestPositions(transcript_info_tibble, codon_of_interest = "TCT")
  
  expand_codon_region <- purrr::map(
    .x = interesting_codon_positions, # vector of codon positions to iterate over for single codons
    .f = ExpandCodonRegion,   # function to use at each codon position of interest for single codons
    transcript_info_tibble,
    gene = "YAL003W",
    dataset = "Mok-simYAL5",
    hd_file = YAL5_h5,
    gff_df,
    expand_width = 5L
  )
  
  return(expand_codon_region)
  
}

expand_codon_region <- ExpandCodonRegionForList(transcript_info_tibble, 
                                                gene = "YAL003W", 
                                                dataset = "Mok-simYAL5", 
                                                hd_file = YAL5_h5, 
                                                startpos = 1, 
                                                startlen = 10, 
                                                codon_of_interest = "TCT", 
                                                gff_df, 
                                                expand_width = 5L)

  # > str(expand_codon_region)
  # List of 6
  # $ : tibble [11 x 4] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 13 14 15 16 17 18 19 20 21 22 ...
  # ..$ Rel_Count: num [1:11] 490 561 51 567 67 65 158 71 323 288 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ : tibble [11 x 4] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 26 27 28 29 30 31 32 33 34 35 ...
  # ..$ Rel_Count: num [1:11] 224 429 201 94 0 54 118 485 242 50 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ : tibble [11 x 4] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 38 39 40 41 42 43 44 45 46 47 ...
  # ..$ Rel_Count: num [1:11] 64 251 128 191 232 89 128 279 417 389 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ : tibble [11 x 4] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 59 60 61 62 63 64 65 66 67 68 ...
  # ..$ Rel_Count: num [1:11] 79 765 196 184 626 84 62 395 104 72 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ : tibble [11 x 4] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 64 65 66 67 68 69 70 71 72 73 ...
  # ..$ Rel_Count: num [1:11] 84 62 395 104 72 214 99 163 128 69 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # $ : tibble [11 x 4] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 192 193 194 195 196 197 198 199 200 201 ...
  # ..$ Rel_Count: num [1:11] 364 960 162 52 115 58 89 235 7 86 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...



### Normalization 

# Normalization carried out within each expanded frame so that they are comparable 
# Normalizes the expand_codon_region list generating a RelCount column with the normalization values
ExpandedCodonRegionNormalization <- function(.x, expand_width = 5L){
  
  # dplyr::mutate(.x, RelCount = PerCodonCounts / sum(PerCodonCounts) * (2 * expand_width + 1))
  dplyr::mutate(.x, RelCount = Rel_Count / sum(Rel_Count) * (2 * expand_width + 1))
  
}


normalized_expanded_codon_region <- ExpandedCodonRegionNormalization(expand_codon_region[[1]], expand_width = 5L)

  # # A tibble: 11 x 5
  # Gene    Pos_Codon Rel_Count Rel_Pos RelCount
  # <chr>       <dbl>     <dbl>   <int>    <dbl>
  #   1 YAL003W        13       490      -5    1.37 
  # 2 YAL003W        14       561      -4    1.57 
  # 3 YAL003W        15        51      -3    0.142
  # 4 YAL003W        16       567      -2    1.58 
  # 5 YAL003W        17        67      -1    0.187
  # 6 YAL003W        18        65       0    0.181
  # 7 YAL003W        19       158       1    0.441
  # 8 YAL003W        20        71       2    0.198
  # 9 YAL003W        21       323       3    0.901
  # 10 YAL003W        22       288       4    0.804
  # 11 YAL003W        23      1301       5    3.63 



# Normalization carried out for all the tibbles within ExpandList 
normalized_expand_list <- purrr::map(
  .x = expand_codon_region,
  .f = ExpandedCodonRegionNormalization,
  expand_width = 5L
)

  # > str(normalized_expand_list)
  # List of 6
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 13 14 15 16 17 18 19 20 21 22 ...
  # ..$ Rel_Count: num [1:11] 490 561 51 567 67 65 158 71 323 288 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 1.367 1.565 0.142 1.582 0.187 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 26 27 28 29 30 31 32 33 34 35 ...
  # ..$ Rel_Count: num [1:11] 224 429 201 94 0 54 118 485 242 50 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 1.067 2.043 0.957 0.448 0 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 38 39 40 41 42 43 44 45 46 47 ...
  # ..$ Rel_Count: num [1:11] 64 251 128 191 232 89 128 279 417 389 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 0.293 1.148 0.585 0.874 1.061 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 59 60 61 62 63 64 65 66 67 68 ...
  # ..$ Rel_Count: num [1:11] 79 765 196 184 626 84 62 395 104 72 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 0.312 3.026 0.775 0.728 2.476 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 64 65 66 67 68 69 70 71 72 73 ...
  # ..$ Rel_Count: num [1:11] 84 62 395 104 72 214 99 163 128 69 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 0.435 0.321 2.044 0.538 0.373 ...
  # $ : tibble [11 x 5] (S3: tbl_df/tbl/data.frame)
  # ..$ Gene     : chr [1:11] "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
  # ..$ Pos_Codon: num [1:11] 192 193 194 195 196 197 198 199 200 201 ...
  # ..$ Rel_Count: num [1:11] 364 960 162 52 115 58 89 235 7 86 ...
  # ..$ Rel_Pos  : int [1:11] -5 -4 -3 -2 -1 0 1 2 3 4 ...
  # ..$ RelCount : num [1:11] 1.686 4.446 0.75 0.241 0.533 ...



### Overlaying the normalized expanded tibbles 

# Function to overlay graphs into a single graph. Need to generate a single tibble 
# from NormalizedExpandList. Need to join by Rel_Pos, in RelCount need the mean for 
# each Rel_Pos (sum row(x) / number of row(x))

OverlayedTable <- function(normalized_expand_list, expand_width = 5L){
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

# something is wrong within the overlaying function, need to figure out why this isnt working anymore  

Over <- OverlayedTable(normalized_expand_list, expand_width = 5L) 

# output from OverlayedTable function for single codons:
# # > Over
# # # A tibble: 11 x 2
# # Rel_Pos RelCount
# # <int>    <dbl>
# #   1      -5    0.860
# # 2      -4    2.09
# # 3      -3    0.876
# # 4      -2    0.735
# # 5      -1    0.772
# # 6       0    0.426
# # 7       1    0.460
# # 8       2    1.21
# # 9       3    0.845
# # 10       4    0.644
# # 11       5    2.08

overlayed_plot <- ggplot(Over, mapping = aes(x = Rel_Pos, y = RelCount)) + geom_line()
