#' This script contains functions called for in YAL5-codon-pairs.R


### Function to filter for a reading frame ###


#' FilterForFrame(): extracts A-site assigned counts, filters for a reading frame
#' for one gene
#'
#' Counts are fetched with GetGeneDatamatrix() and A-site assignment is carried
#' out with CalcAsiteFixed() (in nt positions)
#'
#' The reading frame is then filtered for and the counts are returned
#' (in codon positions)
#'
#' @param gene to get read lengths for
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer, additional displacement in the snapping, default = 0L
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data format,
#' as created by readGFFAsDf() from which to extract the UTRs and CDS widths.
#'
#' @return a list of numeric values (read counts) for a reading frame (0/1/2)
#'
#' @example
#'
#' gff_df <- readGFFAsDf(here::here("..",
#'                                  "example-datasets",
#'                                  "simulated",
#'                                  "mok",
#'                                  "annotation",
#'                                  "Scer_YAL_5genes_w_250utrs.gff3"))
#'
#' FilterForFrame(gene = "YAL003W",
#'                dataset = "Mok-simYAL5",
#'                hd_file,
#'                min_read_length = 10,
#'                snapdisp = 0L,
#'                asite_disp_path = here::here(
#'                  "data", "yeast_standard_asite_disp_length.txt"),
#'                gff_df)
#'
#' @export
FilterForFrame <- function(gene, dataset, hd_file, min_read_length,
                           snapdisp, asite_disp_path, gff_df) {
  
  # Fetch values used for A-site assignment
  asite_disp_length <- suppressMessages(ReadAsiteDisplacementLengthFromFile(asite_disp_path))
  
  # Fetch the read counts for a single gene
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  
  # Calculate A-site assigned counts
  reads_asitepos <- CalcAsiteFixed(reads_pos_length,
                                   min_read_length,
                                   asite_disp_length)
  
  # Fetch GFF values for a gene, e.g. gene = "YAL003W"
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene)
  
  # Assign start position of the CDS, e.g. for YAL003W: 251
  left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                   type == "CDS")
                     %>%  select(start))
  
  # Assign end position of the CDS, e.g. for YAL003W: 871
  right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                    type == "CDS")
                      %>%  select(end))
  
  # Select the reads from the CDS, discards the UTR counts,
  # e.g. for YAL003W: 621 (num [1:621] 811 460 2978 429 251 ...)
  cds <- reads_asitepos[left:right]
  
  # Convert CDS nucleotide length to codon length, e.g. for YAL003W: 621 -> 207
  cds_length <- length(cds) / 3
  
  # Align CDS counts (in nucleotides) to their corresponding reading frame
  cds_frames <- tibble(Count = cds,
                       Frame = rep(c(0, 1, 2),
                                   times = cds_length))
  # > str(cds_frames)
  # tibble [621 x 2] (S3: tbl_df/tbl/data.frame)
  # $ Count: num [1:621] 811 460 2978 429 251 ...
  # $ Frame: num [1:621] 0 1 2 0 1 2 0 1 2 0 ...
  
  # Filter for the frame of interest
  filtered <- dplyr::filter(cds_frames, Frame == snapdisp)
  
  # Select the filtered counts (length(filtered_counts) == length(cds_length))
  filtered_counts <- filtered$Count
  
  return(filtered_counts)
}
# TEST: FilterForFrame(): returns a list of numeric values = TRUE
# TEST: FilterForFrame(): length(filtered_counts) = length(cds_length) for gene
#                         of interest = TRUE
# gives:
# > str(filtered_counts)
# num [1:207] 811 429 488 102 994 146 173 762 13 176 ...



### Functions to fetch and format read counts from the h5 file ###


#' GetAllCodonPosCounts1Gene(): helper function for GetAllCodonPosCounts(),
#' extracts the A-site assigned counts for a single gene 
#'
#' If filter_for_frame = FALSE the GetGeneCodonPosReads1dsnap() function is
#' applied to a gene and generates a tidy data frame (tibble) which contains the
#' counts from all reading frames.
#'
#' If filter_for_frame = TRUE the FilterForFrame() function is applied instead,
#' this generates a tibble which contains the read counts for a reading frame
#' (snapdisp = 0L or 1L or 2L).
#' 
#' @param gene for which to get read counts for
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping;
#' Default = 0L
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data format,
#' as created by readGFFAsDf() from which to extract the UTRs and CDS widths.
#' 
#' @return a tibble which contains the codon positions and counts for a gene  
#' 
#' @export
GetAllCodonPosCounts1Gene <- function(gene,
                                      dataset,
                                      hd_file,
                                      min_read_length,
                                      asite_disp_path,
                                      snapdisp,
                                      filter_for_frame,
                                      gff_df) {
  
  # Fetch values used for A-site assignment
  asite_disp_length <- suppressMessages(ReadAsiteDisplacementLengthFromFile(asite_disp_path))
  
  # Fetch gff values for a gene, e.g. gene = YAL003W
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene)
  
  # Assign start position of the CDS, e.g. for YAL003W: 251
  left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                   type == "CDS")
                     %>%  select(start))
  
  # Assign end position of the CDS, e.g. for YAL003W: 871
  right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                    type == "CDS")
                      %>%  select(end))
  
  # If all reads are to be kept
  if (filter_for_frame == FALSE) {
    
    # Groups counts for all reading frames to their respective codon
    codon_counts_1_gene <- GetGeneCodonPosReads1dsnap(gene,
                                                      dataset,
                                                      hd_file,
                                                      left,
                                                      right,
                                                      min_read_length,
                                                      asite_disp_length,
                                                      snapdisp = snapdisp)
  } else {
    
    # Fetch counts for the reading frame of interest
    codon_counts_1_gene <- FilterForFrame(gene,
                                          dataset,
                                          hd_file,
                                          min_read_length,
                                          snapdisp,
                                          asite_disp_path)
    
  }
  
  # Make a tibble which contains "Gene", "Pos" and "Count"
  codon_pos_counts <- tibble(Gene = gene,
                             Pos = 1:length(codon_counts_1_gene),
                             Count = codon_counts_1_gene)
  
  as.data.frame(codon_pos_counts,
                row.names = NULL,
                optional = FALSE)
  
  return(codon_pos_counts)
  
}
#TEST: GetAllCodonPosCounts1Gene(): retuns a tibble = TRUE
#TEST: GetAllCodonPosCounts1Gene(): the tibble has 3 columns, the names are %in% 
#      c("Gene", "Pos" and "Count").
#TEST: GetAllCodonPosCounts1Gene(): number of observations in the output tibble =
#      sum of CDS (codon co-ordinates) for the gene.
#gives:
# > str(counts_for_one_gene)
# Classes "tbl_df", "tbl" and "data.frame":   207 observations of 3 variables:
# $ Gene : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
# $ Pos  : int  1 2 3 4 5 6 7 8 9 10 ...
# $ Count: num  4249 825 1017 1176 1116 ...



#' GetAllCodonPosCounts(): extracts A-site assigned counts for a list of genes
#'
#' This function iterates over the function GetAllCodonPosCounts1Gene()
#'
#' If filter_for_frame = FALSE the GetGeneCodonPosReads1dsnap() function is
#' applied to a list of genes and generates a tibble which contains the counts
#' from all reading frames.
#'
#' If filter_for_frame = TRUE the FilterForFrame() function is applied instead,
#' this generates a tibble which contains the read counts for a reading frame
#' (snapdisp = 0L or 1L or 2L).
#'
#' @param gene from gene_names to get read lengths for
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping;
#' Default = 0L
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data format,
#' as created by readGFFAsDf() from which to extract the UTRs and CDS widths.
#'
#' @return tibble containing columns "Gene", "Pos" and "Count" for gene_names
#'
#' @example
#' 
#' gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
#' 
#' gff_df <- readGFFAsDf(here::here("..",
#'                                  "example-datasets",
#'                                  "simulated",
#'                                  "mok",
#'                                  "annotation",
#'                                  "Scer_YAL_5genes_w_250utrs.gff3"))
#'
#' GetAllCodonPosCounts(gene_names,
#'                      dataset = "Mok-simYAL5",
#'                      hd_file = here::here("Mok-simYAL5",
#'                                           "output",
#'                                           "A",
#'                                           "A.h5"),
#'                      min_read_length = 10,
#'                      snapdisp = 0L,
#'                      asite_disp_path = here::here(
#'                        "data", "yeast_standard_asite_disp_length.txt"),
#'                      filter_for_frame = TRUE,
#'                      gff_df)
#'
#' @export
GetAllCodonPosCounts <- function(gene_names, dataset, hd_file, min_read_length,
                                 snapdisp, asite_disp_path,
                                 filter_for_frame, gff_df) {
  
  # Apply GetAllCodonPosCounts1Gene() to genes contained within gene_names
  all_codon_pos_counts <- purrr::map_dfr(.x = gene_names,
                                         .f = GetAllCodonPosCounts1Gene,
                                         dataset,
                                         hd_file,
                                         min_read_length,
                                         snapdisp = snapdisp,
                                         asite_disp_path = asite_disp_path,
                                         filter_for_frame = filter_for_frame,
                                         gff_df
  )
  
  return(all_codon_pos_counts)
}
#TEST: GetAllCodonPosCounts(): returns a tibble = TRUE.
#TEST: GetAllCodonPosCounts(): the tibble has 3 columns, the names are %in% 
#      c("Gene", "Pos" and "Count").
#TEST: GetAllCodonPosCounts(): number of observations in the output tibble =
#      sum of CDS (codon co-ordinates) for all genes in gene_names.
#TEST: GetAllCodonPosCounts(): the unique gene names in column "Gene" match
#      the genes in gene_names (unique(all_codon_pos_counts$Gene) =
#      gene_names) = TRUE.
#gives:
# for filter_for_frame = FALSE
# > str(all_codon_pos_counts)
# Classes "tbl_df", "tbl" and "data.frame":   2749 observations of 3 variables:
#   $ Gene  : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos   : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count : num  4249 825 1017 1176 1116 ...
# for filter_for_frame = TRUE, snapdisp = 0L
# > str(all_codon_pos_counts)
# Classes "tbl_df", "tbl" and "data.frame":   2749 observations of 3 variables:
#   $ Gene  : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos   : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count : num  811 429 488 102 994 146 173 762 13 176 ...



#' AddCodonNamesToCodonPosCounts(): takes codon names from the annotation .tsv
#' file, joined to the codon counts table
#'
#' Uses the function GetAllCodonPosCounts().
#'
#' @param .tsv file from which to fetch the codon names associated with the CDS
#' co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created
#' from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data format,
#' as created by readGFFAsDf() from which to extract the UTRs and CDS widths.
#'
#' @return a  tidy format data frame (tibble) which contains the genes in
#' gene_names, codon positions, counts and the codon pair.
#'
#' @example
#' 
#' gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
#'
#' gff_df <- readGFFAsDf(here::here("..",
#'                                  "example-datasets",
#'                                  "simulated",
#'                                  "mok",
#'                                  "annotation",
#'                                  "Scer_YAL_5genes_w_250utrs.gff3"))
#' AddCodonNamesToCodonPosCounts(yeast_codon_pos_i200,
#'                               gene_names,
#'                               dataset = "Mok-simYAL5",
#'                               hd_file = here::here("Mok-simYAL5",
#'                                                    "output",
#'                                                    "A",
#'                                                    "A.h5"),
#'                               min_read_length = 10,
#'                               snapdisp = 0L,
#'                               asite_disp_path = here::here(
#'                                "data", "yeast_standard_asite_disp_length.txt"),
#'                               filter_for_frame = TRUE,
#'                               gff_df)
#'
#' @export
AddCodonNamesToCodonPosCounts <- function(yeast_codon_pos_i200, gene_names,
                                          dataset, hd_file, min_read_length,
                                          snapdisp, asite_disp_path,
                                          filter_for_frame, gff_df) {
  
  yeast_codon_pos_i200 <- suppressMessages(readr::read_tsv
                                           (file = codon_pos_table))
  
  # The yeast_codon_pos_i200 file is configured to show codon pair positions:
  codon_pos_pair_i200 <- tibble::tibble(
    Gene = yeast_codon_pos_i200$Gene,
    Pos = yeast_codon_pos_i200$PosCodon,
    CodonPair = paste(yeast_codon_pos_i200$Codon,
                      (dplyr::lead(yeast_codon_pos_i200$Codon)
                       %>% str_replace_all("ATG", "NA")))
  )
  
  # Fetch read counts and positions for each gene contained within gene_names
  all_codon_pos_counts <- GetAllCodonPosCounts(gene_names,
                                               dataset,
                                               hd_file,
                                               min_read_length,
                                               snapdisp,
                                               filter_for_frame = filter_for_frame,
                                               asite_disp_path,
                                               gff_df)
  
  # Join .tsv file with all_codon_pos_counts
  transcript_tibbles <- all_codon_pos_counts %>% left_join(codon_pos_pair_i200,
                                                           by = c("Pos" = "Pos",
                                                                  "Gene" = "Gene"),
                                                           keep = FALSE,
                                                           copy = TRUE)
  
  # Make final tibble with correct column names
  gene_pos_codon_counts <- tibble(Gene = transcript_tibbles$Gene,
                                  Pos = transcript_tibbles$Pos,
                                  CodonPair = transcript_tibbles$CodonPair,
                                  Count = transcript_tibbles$Count)
                                  
  
  return(gene_pos_codon_counts)
}
#TEST: AddCodonNamesToCodonPosCounts(): creates a tibble = TRUE
#TEST: AddCodonNamesToCodonPosCounts(): the tibble contains 4 columns = TRUE
#TEST: AddCodonNamesToCodonPosCounts(): number of observations in the output
#      tibble = sum of CDS (codon co-ordinates) for all genes in gene_names.
#TEST: AddCodonNamesToCodonPosCounts(): the column names are %in%
#      c("Gene", "Pos", "Count", "CodonPair")
#TEST: AddCodonNamesToCodonPosCounts(): the unique gene names in column "Gene"
#      match the genes in gene_names (unique(all_codon_pos_counts$Gene) =
#      gene_names) = TRUE.
# gives:
# > str(gene_pos_codon_counts)
# Classes "tbl_df", "tbl" and "data.frame":   2,749 observations of 4 variables:
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos      : num  1 2 3 4 5 6 7 8 9 10 ...
#   $ CodonPair: chr  "ATG GCA" "GCA TCC" "TCC ACC" "ACC GAT" ...
#   $ Count    : num  4249 825 1017 1176 1116 ...



#' FilterForFeaturePositions(): Filters for the feature of interest
#'
#' Uses the AddCodonNamesToCodonPosCounts() to get counts aligned to positions
#' and codons
#'
#' @param .tsv file from which to fetch the codon names associated with the CDS
#' co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon
#' @param feature_of_interest character, each incidence of the feature will be
#' extracted from transcript_info_tibble
#'
#' @return a tidy format data frame (tibble) which contains the columns "Gene",
#' "CodonPos1", "CodonPos2", "Count" and "CodonPair" for occurrences of the
#' feature of interest.
#'
#' @example
#'
#' FilterForFeaturePositions(yeast_codon_pos_i200,
#'                           gene,
#'                           dataset = "Mok-simYAL5",
#'                           hd_file = here::here("Mok-simYAL5",
#'                                                "output",
#'                                                "A",
#'                                                "A.h5"),
#'                           min_read_length = 10,
#'                           snapdisp = 0L,
#'                           asite_disp_path = here::here(
#'                              "data", "yeast_standard_asite_disp_length.txt"),
#'                           filter_for_frame = TRUE,
#'                           feature_of_interest = "TCC AAG")
#'
#' @export
FilterForFeaturePositions <- function(yeast_codon_pos_i200,
                                      gene_names,
                                      dataset,
                                      hd_file,
                                      min_read_length,
                                      snapdisp,
                                      asite_disp_path,
                                      filter_for_frame,
                                      feature_of_interest,
                                      gff_df) {
  
  # make a tibble containing the counts for each gene in gene names and aligned
  # to their positions and codons
  gene_pos_codon_counts <- AddCodonNamesToCodonPosCounts(yeast_codon_pos_i200,
                                                         gene_names,
                                                         dataset,
                                                         hd_file,
                                                         min_read_length,
                                                         snapdisp,
                                                         asite_disp_path,
                                                         filter_for_frame,
                                                         gff_df = gff_df)
  
  # filter for the occurrences of the feature of interest
  interesting_feature_table <- dplyr::filter(gene_pos_codon_counts,
                                             CodonPair == feature_of_interest)
  
  return(interesting_feature_table)
}
#TEST: FilterForFeaturePositions(): returns a tibble = TRUE
#TEST: FilterForFeaturePositions(): the tibble contains 4 columns = TRUE
#TEST: FilterForFeaturePositions(): the column names are %in%
#      c("Gene", "Pos", "Count", "CodonPair")
#TEST: FilterForFeaturePositions(): the unique feature of interest in
#      column "CodonPair" (unique(interesting_feature_positions$CodonPair)) = 1
#gives:
# > str(interesting_feature_positions)
# Classes "tbl_df", "tbl" and "data.frame":   8 observations of 5 variables:
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL005C" "YAL005C" ...
#   $ Pos      : num  7 57 383 508 535 321 90 412
#   $ Count    : num  1553 1147 358 317 498 ...
#   $ CodonPair: chr  "TCC AAG" "TCC AAG" "TCC AAG" "TCC AAG" ...



### Slice out and expand a frame around interesting features ###



#' TranscriptForOneGene(): helper function for  ExpandFeatureRegion()
#' 
#' Takes gene_pos_codon_counts, generated by AddCodonNamesToCodonPosCounts()
#' as an input in addition to a gene the feature of interest 
#' 
#' @param gene in gene names used to extract read counts from h5 file
#' @param gene_pos_codon_counts generated by AddCodonNamesToCodonPosCounts()
#' @param feature_of_interest character, each of incidence of the feature is extracted
#' 
#' @return a tibble for 
TranscriptForOneGene <- function(gene,
                                 gene_pos_codon_counts,
                                 feature_of_interest) {
  
  # Filter for the occurrences of the feature_of_interest
  interesting_feature_tibble <- dplyr::filter(gene_pos_codon_counts,
                                              CodonPair == feature_of_interest) 
  
  # Filters for the gene of interest from interesting_feature_tibble
  transcript_for_one_gene <- dplyr::filter(interesting_feature_tibble,
                                           Gene == gene)
  
  return(transcript_for_one_gene)
}
#TEST: TranscriptForOneGene(): returns a tibble = TRUE
#TEST: TranscriptForOneGene(): the tibble contains 4 columns = TRUE
#TEST: TranscriptForOneGene(): the column names are %in% "Gene", "Pos", "CodonPair"
#      and "Count"
#gives:
# > str(transcript_gene)
# Classes "tbl_df", "tbl" and "data.frame":   3 observations of 4 variables
#   $ Gene     : chr  "YAL035W" "YAL035W" "YAL035W"
#   $ Pos      : num  539 714 794
#   $ CodonPair: chr  "AAA GCC" "AAA GCC" "AAA GCC"
#   $ Count    : num  80 53 52



#' ExpandFeatureRegion(): Generates a window around each feature_of_interest
#' across all genes in gene_names
#' 
#' Feature_of_interest has position 0, all adjacent codons are assigned relative
#' positions.
#'
#' @param .tsv file from which to fetch the codon names associated with the CDS
#' co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon
#' @param feature_of_interest character, each incidence of the feature will be
#' extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions on each
#' side of the feature_of_interest to include in the window
#'
#' @return a list of tibbles, one tibble for each occurrence of the feature of
#' interest with an expanded window around the feature
#'
#' @example
#' 
#' gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
#'
#' gff_df <- readGFFAsDf(here::here("..",
#'                                  "example-datasets",
#'                                  "simulated",
#'                                  "mok",
#'                                  "annotation",
#'                                  "Scer_YAL_5genes_w_250utrs.gff3"))
#'
#' ExpandFeatureRegion(yeast_codon_pos_i200,
#'                     gene_names,
#'                     dataset = "Mok-simYAL5",
#'                     hd_file = here::here("Mok-simYAL5",
#'                                          "output",
#'                                          "A",
#'                                          "A.h5"),
#'                     gff_df,
#'                     min_read_length = 10,
#'                     snapdisp = 0L,
#'                     asite_disp_path = here::here(
#'                       "data", "yeast_standard_asite_disp_length.txt"),
#'                     filter_for_frame = TRUE,
#'                     feature_of_interest = "TCC AAG",
#'                     expand_width = 5L)
#'
#' @export
ExpandFeatureRegion <- function(yeast_codon_pos_i200, gene_names, dataset,
                                hd_file, gff_df, min_read_length, snapdisp,
                                asite_disp_path, filter_for_frame,
                                feature_of_interest, expand_width) {
  
  # Fetch the assigned read counts for genes in gene_names
  gene_pos_codon_counts <- AddCodonNamesToCodonPosCounts(yeast_codon_pos_i200,
                                                         gene_names,
                                                         dataset,
                                                         hd_file,
                                                         min_read_length,
                                                         snapdisp,
                                                         asite_disp_path,
                                                         filter_for_frame,
                                                         gff_df)
  
  # Takes gene_pos_codon_counts as inputs and select for positions
  # on separate genes
  AllGeneInterestingFeatures <- function(yeast_codon_pos_i200,
                                         gene, gene_names, dataset, hd_file,
                                         min_read_length, feature_of_interest,
                                         gene_pos_codon_counts,
                                         gff_df) {
    
    transcript_for_one_gene <- TranscriptForOneGene(gene,
                                                    gene_pos_codon_counts,
                                                    feature_of_interest)
    
    # Expand the region around an occurrence of the feature_of_interest
    ExpandRegions <- function(transcript_for_one_gene,
                              gene_pos_codon_counts,
                              gene_names,
                              hd_file,
                              gff_df,
                              expand_width) {
      
      # assigns new name to the transcript_for_one_gene tibble
      interesting_features <- transcript_for_one_gene
      
      # filters gene_pos_codon_counts down to the gene being processed
      gene_pos_codon_feature <- dplyr::filter(gene_pos_codon_counts,
                                              Gene == gene)
      
      # fetch gene length for gene being processed from the gff file
      gene_length <- dplyr::filter(gff_df,
                                   gff_df$type == "CDS" &
                                     gff_df$Name == gene) %>% select(width)
      
      # if the window cannot be formed around the feature_of_interest
      if (interesting_features <= expand_width  |
          interesting_features + expand_width > gene_length / 3) {
        
        return()
        
        # if the window can be formed around the feature_of_interest
      } else {
        
        # slice and return the window around the feature_of_interest
        expand_feature_region <- tibble(dplyr::slice(gene_pos_codon_feature,
                                                     (interesting_features - expand_width):(interesting_features + expand_width),
                                                     each = FALSE),
                                        RelPos =  seq(- expand_width, expand_width)
        )
        
        if (dim(expand_feature_region)[1] == (2 * expand_width + 1)) {
          
          return(expand_feature_region)
          
        } else {
          
          return()
        }
      }
    }
    # The if statement ensures that feature positions that are less/more than
    # the expand_width value are discarded
    
    expand_feature_region <- purrr::map(.x = transcript_for_one_gene$CodonPos1,
                                        .f = ExpandRegions,
                                        gene_pos_codon_counts,
                                        gene,
                                        gene_names = gene_names,
                                        gff_df = gff_df,
                                        expand_width = expand_width)
    
    return(expand_feature_region)
  }
  
  expand_feature_region <- purrr::map(.x = gene_names,
                                      .f = AllGeneInterestingFeatures,
                                      gene_names = gene_names,
                                      yeast_codon_pos_i200 = yeast_codon_pos_i200,
                                      dataset,
                                      hd_file,
                                      min_read_length,
                                      feature_of_interest,
                                      gene_pos_codon_counts,
                                      gff_df
  )
  
  expand_feature_region <- unlist(expand_feature_region, recursive = FALSE)
  
  # remove NULLS, which represent features of interest occurring within one
  # expand width of the UTRs
  expand_feature_region <- expand_feature_region[!sapply(expand_feature_region,
                                                         is.null)]
  
  return(expand_feature_region)
}
#TEST: ExpandFeatureRegion(): output is a list of tidy format data frames
#      (tibbles) = TRUE. type(expand_feature_region) = "list"
#TEST: ExpandFeatureRegion(): number of tibbles in list matches the number
#      of occurrences in "feature_of_interest" list = TRUE
#TEST: ExpandFeatureRegion(): each tibble contains 6 columns = TRUE
#TEST: ExpandFeatureRegion(): the column names are %in% c("Gene",
#      "Pos", "Count", "CodonPair", "RelPos")
#TEST: ExpandFeatureRegion(): number of observations in each output tibble =
#      "expand_width"*2+1, if "expand_width" = 5L the number of observations
#      should be 11
#TEST: ExpandFeatureRegion(): the position from "interesting_feature_positions"
#      has "RelPos" value 0 = TRUE
#TEST: ExpandFeatureRegion(): the column "RelPos" goes from -"expand_width
#      to +"expand_width"
#gives:
# > str(expand_feature_region)
# List of 8 (showing 2)
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 6 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos      : num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Count    : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ CodonPair: chr  "GCA TCC" "TCC ACC" "ACC GAT" "GAT TTC" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 6 variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos      : num  52 53 54 55 56 57 58 59 60 61 ...
#   $ Count    : num  42 53 648 293 121 92 519 79 765 196 ...
#   $ CodonPair: chr  "TTC AAC" "AAC CAC" "CAC ATC" "ATC GCT" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# ...


### Normalisation ###


#' ExpandedRegionNormalisation(): carries out normalization within each expanded
#' frame so that they are comparable
#'
#' Normalizes the ExpandFeatureRegion() list, generates a RelCount column with
#' the normalization values
#'
#' @param .tsv file from which to fetch the codon names associated with the CDS
#' co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data format,
#' as created by readGFFAsDf() from which to extract the UTRs and CDS widths.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon.
#' @param feature_of_interest character, each incidence of the feature will be
#' extracted from transcript_info_tibble
#' @param expand_width integer, provides number of positions each side of the 
#' feature_of_interest to include in the window
#'
#' @return a list of normalized tibbles, where each tibble is an occurrence of
#' the feature_of_interest
#'
#' @example
#'
#' gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
#' 
#' gff_df <- readGFFAsDf(here::here("..",
#'                                  "example-datasets",
#'                                  "simulated",
#'                                  "mok",
#'                                  "annotation",
#'                                  "Scer_YAL_5genes_w_250utrs.gff3"))
#'
#' ExpandedRegionNormalisation(yeast_codon_pos_i200,
#'                             gene_names,
#'                             dataset = "Mok-simYAL5",
#'                             hd_file = here::here("Mok-simYAL5",
#'                                                   "output",
#'                                                   "A",
#'                                                   "A.h5"),
#'                             gff_df,
#'                             min_read_length = 10,
#'                             snapdisp = 0L,
#'                             asite_disp_path = here::here(
#'                               "data", "yeast_standard_asite_disp_length.txt"),
#'                             filter_for_frame = TRUE,
#'                             feature_of_interest = "TCC AAG",
#'                             expand_width = 5L)
#'
#' @export
ExpandedRegionNormalisation <- function(yeast_codon_pos_i200,
                                        gene_names, dataset, hd_file,
                                        gff_df, min_read_length, snapdisp,
                                        asite_disp_path,
                                        filter_for_frame,
                                        feature_of_interest,
                                        expand_width) {
  
  # fetch the expanded tibbles for each occurrence of the feature_of_interest
  expand_feature_region <- ExpandFeatureRegion(yeast_codon_pos_i200,
                                               gene_names,
                                               dataset,
                                               hd_file,
                                               gff_df,
                                               min_read_length,
                                               snapdisp,
                                               asite_disp_path,
                                               filter_for_frame,
                                               feature_of_interest,
                                               expand_width)
  
  # function to normalise each tibble before they can be overlayed
  Normalization <- function(.x, expand_width) {
    dplyr::mutate(.x, RelCount = Count / sum(Count) * (2 * expand_width + 1))
  }
  
  normalized_expand_list <- purrr::map(.x = expand_feature_region,
                                       .f = Normalization,
                                       expand_width)
  
  # function to remove the occurrences of tibbles which contain NaN
  CheckForNaN <- function(.x) {
    
    normalized_expand_list <- .x
    
    Relcount_values <- unlist(normalized_expand_list$RelCount)
    
    SetNaNToZero <- function(Relcount_values) {
      
      if (is.nan(Relcount_values)) {
        print("NaN present")
        Relcount_values <- 0
      } else {
        Relcount_values <- Relcount_values
      }
    }
    Relcount_values <- unlist(purrr::map(.x = Relcount_values,
                                         .f = SetNaNToZero))
    
    dplyr::mutate(.x, RelCount = Relcount_values)
  }
  
  normalized_expand_list <- purrr::map(.x = normalized_expand_list,
                                       .f = CheckForNaN)
  
}
#TEST: ExpandedRegionNormalisation(): creates a list of tidy format data frame
#      (tibble) (type(normalized_expand_list)) = TRUE
#TEST: ExpandedRegionNormalisation(): the tibble contains 6 columns = TRUE
#TEST: ExpandedRegionNormalisation(): the column names are %in% c("Gene",
#      "Pos", "Count", "RelPos", "RelCount")
#TEST: ExpandedRegionNormalisation(): number of observations in the output
#      tibble = "expand_width"*2+1, if "expand_width" = 5L the number of
#      observations should be 11
#TEST: ExpandedRegionNormalisation(): the column "RelPos" goes from
#      -"expand_width to +"expand_width"
#gives:
# > str(normalized_expanded_feature_region)
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 6 variables
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ Pos       : num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Count     : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ RelPos    : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
#   $ RelCount  : num  1.347 1.532 0.32 3.12 0.458 ...


### Overlaying the normalized expanded tibbles ###


#' OverlayedTibble: overlays the expanded tibbles generated by
#' ExpandedRegionNormalisation()
#'
#' Uses the ExpandedRegionNormalisation() function to fetch counts before each
#' tibble is normalised internally
#'
#' @param .tsv file from which to fetch the codon names associated with the CDS
#' co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created
#' from BAM files for dataset samples.
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data format,
#' as created by readGFFAsDf() from which to extract the UTRs and CDS widths.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon
#' @param feature_of_interest character, each incidence of the feature will be
#' extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions on each
#' side of the feature_of_interest to include in the window
#'
#' @return  an overlayed tidy format data frame (tibble) consisting of all the
#' normalised tibbles.
#'
#' @example
#' 
#' gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
#' 
#' gff_df <- readGFFAsDf(here::here("..",
#'                                  "example-datasets",
#'                                  "simulated",
#'                                  "mok",
#'                                  "annotation",
#'                                  "Scer_YAL_5genes_w_250utrs.gff3"))
#'
#' OverlayedTibble(yeast_codon_pos_i200,
#'                  gene_names,
#'                  dataset = "Mok-simYAL5",
#'                  hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"),
#'                  gff_df,
#'                  min_read_length = 10,
#'                  snapdisp = 0L,
#'                  asite_disp_path = here::here(
#'                     "data", "yeast_standard_asite_disp_length.txt"),
#'                  filter_for_frame = TRUE,
#'                  feature_of_interest = "TCC AAG",
#'                  expand_width = 5L)
#'
#' @export
OverlayedTibble <- function(yeast_codon_pos_i200,
                            gene_names, dataset, hd_file,
                            gff_df, min_read_length,
                            snapdisp, asite_disp_path,
                            filter_for_frame, feature_of_interest,
                            expand_width
) {
  
  # fetch a list of normalised tibbles consisting of each occurrence of the
  # feature_of_interest
  normalized_expand_list <- ExpandedRegionNormalisation(yeast_codon_pos_i200,
                                                        gene_names, dataset,
                                                        hd_file, gff_df,
                                                        min_read_length,
                                                        snapdisp,
                                                        asite_disp_path,
                                                        filter_for_frame,
                                                        feature_of_interest,
                                                        expand_width)
  
  # the number of objects inside normalized_expand_list
  number_of_objects <- length(normalized_expand_list)
  
  # reduces normalized_expand_list to the columns RelPos and RelCount
  result <- lapply(normalized_expand_list, "[", c("RelPos", "RelCount"))
  
  # sums the columns RelPos and RelCount for all the tibbles into one tibble
  joined_result <- purrr::reduce(.x = result, .f = `+`)
  
  # average the tibble by dividing by the number of objects summed
  joined_rows <- mutate(joined_result / number_of_objects)
  
  # create the final overlayed tibble with the relative positions and the
  # normalised and overlayed relative counts around the feature_of_interest
  overlayed_tibbles <- tibble::tibble(RelPos = seq(- expand_width, expand_width),
                                      RelCount = joined_rows$RelCount)
  
  return(overlayed_tibbles)
}
#TEST: OverlayedTibble(): creates a tidy format data frame (tibble) = TRUE
#TEST: OverlayedTibble(): the tibble contains 2 columns = TRUE
#TEST: OverlayedTibble(): the column names are %in% c("RelPos", "RelCount")
#TEST: OverlayedTibble(): number of observations in the output tibble =
#      "expand_width"*2+1, if "expand_width" = 5L the number of observations
#      should be 11
#TEST: OverlayedTibble(): the column "RelPos" goes from -"expand_width to +
#      "expand_width"
#TEST: OverlayedTibble(): RelCount is a numeric
#gives:
# > str(overlayed_tibbles)
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 2 variables
#  $ RelPos   : int -5 -4 -3 -2 -1 0 1 2 3 4 ...
#  $ RelCount : num  1.029 1.07 0.987 1.45 1.151 ...



#' FindAllFeatures(): extracts the RelCount values for multiple features of interest 
#' 
#' The feature contains the functions ExpandFeatureRegionAllGenes(),
#' ExpandedRegionNormalization() and OverlayedTibble(), which are defined above.
#' 
#' @param .tsv file from which to fetch the codon names associated with the CDS
#' co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes, created
#' from BAM files for dataset samples. 
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml) 
#' @param gff_df from which to extract the UTRs and CDS widths. 
#' @param asite_disp_path integer, lengths used for A-site assignment 
#' @param expand_width integer which provides the number of positions on each
#' side of the feature_of_interest to include in the window 
#' @param filter_for_frame logical, TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon
#' @param feature_of_interest character, each incidence of the feature will be
#' extracted from transcript_info_tibble
#' @param snapdisp integer any additional displacement in the snapping 
#' 
#' @return a tidy format data frame (tibble) consisting of the columns "Feature" 
#' and "RelCount" containing all objects listed in feature_of_interest and their 
#' relative counts
#' 
#' @example 
#' 
#' gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
#' 
#' gff_df <- readGFFAsDf(here::here("..",
#'                                  "example-datasets",
#'                                  "simulated",
#'                                  "mok",
#'                                  "annotation",
#'                                  "Scer_YAL_5genes_w_250utrs.gff3"))
#' 
#' FindAllFeatures(yeast_codon_pos_i200,
#'                 gene_names,
#'                 dataset = "Mok-simYAL5",
#'                 hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"),
#'                 gff_df,
#'                 min_read_length = 10,
#'                 snapdisp = 0L,
#'                 asite_disp_path = here::here(
#'                     "data", "yeast_standard_asite_disp_length.txt"),
#'                 filter_for_frame = TRUE,
#'                 feature_of_interest = "TCC AAG",
#'                 expand_width = 5L)
#'                 
#'@export
FindAllFeatures <- function(yeast_codon_pos_i200,
                            gene_names, dataset, hd_file,
                            min_read_length, feature_of_interest,
                            gff_df, asite_disp_path,
                            expand_width, filter_for_frame, snapdisp) {
  
  
  # Run ExpandFeatureRegionAllGenes() to get a list of occurrences of the
  # feature_of_interest
  
  print(paste0("Finding occurrences of ", feature_of_interest))
  
  expanded_feature <- suppressMessages(ExpandFeatureRegion(yeast_codon_pos_i200,
                                                           gene_names,
                                                           dataset,
                                                           hd_file,
                                                           gff_df,
                                                           min_read_length,
                                                           snapdisp,
                                                           asite_disp_path,
                                                           filter_for_frame,
                                                           feature_of_interest,
                                                           expand_width))
  
  # Check for the presence of the feature_of_interest. expanded_feature
  # being empty will cause problems with normalization
  
  if (length(expanded_feature) == 0) {
    
    print(paste("No occurrences of", feature_of_interest))
    
    if (expand_width > 1) {
      
      print("Use expand_width = 1L to check for occurrences near to start or stop codon")
    }
    return()
  }
  
  # Run ExpandedRegionNormalisation() to calculate the relative number of
  # reads mapping to each position around the feature_of_interest
  
  normalized_expand_list <- ExpandedRegionNormalisation(yeast_codon_pos_i200,
                                                        gene_names, dataset,
                                                        hd_file, gff_df,
                                                        min_read_length,
                                                        snapdisp,
                                                        asite_disp_path,
                                                        filter_for_frame,
                                                        feature_of_interest,
                                                        expand_width)
  
  # Run OverlayedTibble() to create an average of reads at positions at and
  # around the feature_of_interest
  
  overlayed_tibbles <- OverlayedTibble(yeast_codon_pos_i200,
                                       gene_names, dataset, hd_file,
                                       gff_df, min_read_length,
                                       snapdisp, asite_disp_path,
                                       filter_for_frame,
                                       feature_of_interest,
                                       expand_width)
  
  # Create a new tibble listing the feature being studeied, and the RelCount
  # at position 0, ie RelCount at the feature_of_interest
  
  feature_rel_use <- tibble(Feature = feature_of_interest,
                            RelCount = filter(overlayed_tibbles,
                                              overlayed_tibbles$RelPos == 0)$RelCount)
}
#TEST: FindAllFeatures(): creates a tibble = TRUE
#TEST: FindAllFeatures(): the tibble contains 2 columns = TRUE
#TEST: FindAllFeatures(): the column names are %in% c("Feature", "RelCount")
#TEST: FindAllFeatures(): number of observations in the output tibble:
#      length(feature_of_interest) == length(feature_rel_use) 
#gives:
# > str(feature_rel_use)
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 2 variables
# $ Feature : chr [1:11] "CCA TGG" "AGA TGG" "GTA GTG" "TCA TAC" ...
# $ RelCount: num [1:11] 3.41 3.35 2.19 2.12 2.11 ...




### Generate plot around feature_of_interest ###


#' GeneratePlot(): Generates a metafeature plot around the feature_of_interest
#'
#' Fetches the overlayed tibble using the function overlayed_tibbles()
#'
#' @param .tsv file from which to fetch the codon names associated with the CDS
#' co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame, FALSE if
#' keeping and grouping all reading frames for each codon
#' @param feature_of_interest character, each incidence of the feature will be
#' extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions on each
#' side of the feature_of_interest to include in the window
#' @param size the size of the text on the metafeature plot
#'
#' @return A plot which shows the ribosomal occupancy around a feature of
#' interest
#'
#' @example
#'
#' gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name
#' 
#' gff_df <- readGFFAsDf(here::here("..",
#'                                  "example-datasets",
#'                                  "simulated",
#'                                  "mok",
#'                                  "annotation",
#'                                  "Scer_YAL_5genes_w_250utrs.gff3"))
#'
#' GeneratePlot(yeast_codon_pos_i200,
#'              gene_names,
#'              dataset = "Mok-simYAL5",
#'              hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"),
#'              gff_df,
#'              min_read_length = 10,
#'              snapdisp = 0L,
#'              asite_disp_path = here::here(
#'                "data", "yeast_standard_asite_disp_length.txt"),
#'              filter_for_frame = TRUE,
#'              feature_of_interest = "TCC AAG",
#'              expand_width = 5L,
#'              size = 12)
#'
#' @export
GeneratePlot <- function(yeast_codon_pos_i200, gene_names, dataset, hd_file,
                         gff_df, min_read_length, snapdisp, asite_disp_path,
                         filter_for_frame, feature_of_interest,
                         expand_width, size = 12) {
  
  overlayed_tibbles <- OverlayedTibble(yeast_codon_pos_i200,
                                       gene_names, dataset, hd_file,
                                       gff_df, min_read_length,
                                       snapdisp, asite_disp_path,
                                       filter_for_frame, feature_of_interest,
                                       expand_width)
  
  overlayed_plot <- ggplot(overlayed_tibbles,
                           mapping = aes(x = RelPos,
                                         y = RelCount)) +
    geom_line() +
    theme_bw() +
    theme(text = element_text(size = size),
          axis.title = element_text(size = size, face = "bold"),
          title = element_text(size = size, face = "bold")) +
    labs(title = paste0("Meta-feature plot of codon pair ",
                        feature_of_interest),
         x = "Position relative to feature of interest",
         y = "Normalised ribosomal occupancy", size = 2) +
    scale_x_continuous(breaks = seq(-expand_width, expand_width, 1))
}
#TEST: GeneratePlot(): produces a single plot =TRUE
#TEST: GeneratePlot(): title is "Meta-feature plot of codon pair
#      <feature_of_interest>
#TEST: GeneratePlot(): the x-axis is "Distance from codon pair (3nt) and y-axis
#      is "Normalised ribosomal occupancy"
