# This script uses a H5, GFF and TSV file to create either a
# metafeature plot for a feature of interest, or a table comparing the
# RelCount of multiple features of interest. 
# Currently designed for single codons.

## TEST::run on TinySim Dataset

suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))

# Load local dependencies.
if (interactive()) {
  # Use hard-coded script name and assume script is in "rscripts"
  # directory. This assumes that interactive R is being run within
  # the parent of rscripts/ but imposes no other constraints on
  # where rscripts/ or its parents are located.
  self <- "YAL5-single-codons.R"
  path_to_self <- here("rscripts", self)
  source(here::here("rscripts", "provenance.R"))
  source(here::here("rscripts", "read_count_functions.R"))
  source(here::here("rscripts", "stats_figs_block_functions.R"))
} else {
  # Deduce file name and path using reflection as before.
  self <- getopt::get_Rscript_filename()
  path_to_self <- self
  source(file.path(dirname(self), "provenance.R"))
  source(file.path(dirname(self), "read_count_functions.R"))
  source(file.path(dirname(self), "stats_figs_block_functions.R"))
}

### Functions to map number of reads to codon positions ###

#' Map reads to codons.
#'
#' Unlike `SnapToCodon` in `read_count_functions.R`, this function
#' only returns reads mapping to the first nucleotide of a codon for
#' the desired frame.
#'
#' @param gene Gene name to pull out read counts for.
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param min_read_length Minimum read length in H5 file (integer).
#' @param asite_displacement_length Length from the end of a read to
#' the A-site.
#' @param filtering_frame Frame from which to select reads from.
#' @return List of numeric values (read counts) for a reading frame
#' (0/1/2).
#'
#' @examples
#' FilterForFrameFunction(gene, gff_df, dataset, hd_file,
#'   min_read_length, asite_displacement_length, filtering_frame)
FilterForFrameFunction <- function(
  gene, gff_df, dataset, hd_file, min_read_length,
  asite_displacement_length, filtering_frame) {

  # Get the matrix of read counts.
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  # Assign reads to their A site nt, determined by read length.
  reads_asitepos <- CalcAsiteFixed(reads_pos_length,
                                   min_read_length,
                                   asite_displacement_length)
  # Get the gff rows for the gene being studied.
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df,
                                         seqnames == gene)
  # Get the position of the start codon.
  left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                   type == "CDS") %>%  select(start))
  # Get the position of the stop codon.
  right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                    type == "CDS") %>%  select(end))
  # Extract the reads that map to the CDS.
  cds_reads <- reads_asitepos[left:right]
  cds_length <- length(cds_reads) / 3
  # Create a tibble, assigning a frame to each nt, so the first nt in
  # each frame has the corresponding frame identity.
  cds_reads_frames <- tibble(Count = cds_reads,
                       Frame = rep(c(0, 1, 2), times = cds_length)
  )
  filtered_for_frame <- dplyr::filter(cds_reads_frames,
                                      Frame == filtering_frame)
  filtered_counts <- filtered_for_frame$Count
  return(filtered_counts)
}

#' Calculate read counts at each codon position for a gene.
#'
#' @param gene Gene in the sample.
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param asite_disp_path Path to A-site file.
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param min_read_length Minimum read length in H5 file (integer).
#' @param filter_for_frame Method to use for assigning reads (logical).
#' @param filtering_frame Frame from which to select reads from.
#' @param snapdisp Frame to filter to.
#' @return Read counts at each codon position.
GetAllCodonPosCounts1Gene <- function(
  gene, gff_df, asite_disp_path, dataset, hd_file, min_read_length,
  filter_for_frame, filtering_frame, snapdisp) {

  # Get the gff rows for the gene being studied.
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene)
  # Get the position of the start codon.
  left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                   type == "CDS") %>%  select(start))
  # Get the position of the stop codon.
  right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                    type == "CDS") %>%  select(end))
  # Get the A-site displacement for reads of different lengths.
  asite_displacement_length <- suppressMessages(
    ReadAsiteDisplacementLengthFromFile(asite_disp_path))
  # Calculate reads at codon positions in the desired way.
  if (! filter_for_frame) {
    # Assign reads to positions, adds the reads from all three
    # nucleotides in the codon.
    codon_counts_1_gene <- GetGeneCodonPosReads1dsnap(
      gene, dataset, hd_file, min_read_length, asite_displacement_length,
      left, right, snapdisp)
  } else {
    # Assign reads to positions, using only the reads from the
    # nucleotide in the codon of a designated frame.
    codon_counts_1_gene <- FilterForFrameFunction(
      gene, gff_df, dataset, hd_file, min_read_length,
      asite_displacement_length, filtering_frame)
  }
  # Add codon positions.
  codon_pos_counts <- tibble(Gene = gene,
                             PosCodon = 1:length(codon_counts_1_gene),
                             Count = codon_counts_1_gene)
  as.data.frame(codon_pos_counts, row.names = NULL, optional = FALSE)
  return(codon_pos_counts)
}

#' Calculate read counts at each codon position for a list of genes.
#'
#' Applies `GetGeneCodonPosReads1dsnap` to a list of genes and
#' generates a tidy data frame (tibble) which contains the counts for
#' all genes
#'
#' @param gene_names List of all the genes in the sample,
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param asite_disp_path Path to A-site file.
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param min_read_length Minimum read length in H5 file (integer).
#' @param filter_for_frame Method to use for assigning reads (logical).
#' @param filtering_frame Frame from which to select reads from.
#' @param snapdisp Frame to filter to.
#' @return A tibble listing the `Gene`, `PosCodon` and `Count`.
#' @export
#'
#' @examples
#'  total_codon_pos_counts <-
#'    GetAllCodonPosCounts(gene_names, gff_df, dataset, hd_file,
#'    min_read_length, filter_for_frame, filtering_frame,  snapdisp)
GetAllCodonPosCounts <- function(
  gene_names, gff_df, asite_disp_path, dataset, hd_file,
  min_read_length, filter_for_frame, filtering_frame, snapdisp) {

  # Apply GetAllCodonPosCounts1Gene to all genes in the sample.
  total_codon_pos_counts <- purrr::map_dfr(
    .x = gene_names,
    .f = GetAllCodonPosCounts1Gene,
    gff_df,
    asite_disp_path,
    dataset,
    hd_file,
    min_read_length,
    filtering_frame = filtering_frame,
    filter_for_frame = filter_for_frame,
    snapdisp = snapdisp
  )
  return(total_codon_pos_counts)
}

# TEST::GetAllCodonPosCounts(): returns a tibble.
# TEST::GetAllCodonPosCounts(): the tibble has 3 columns.
# TEST::GetAllCodonPosCounts(): the column names are %in%
# c("Gene", "PosCodon" and "Count").
# TEST::GetAllCodonPosCounts(): number of observations in the output
# tibble = sum of CDS (codon co-ordinates) for all genes in
# gene_names.
# TEST::GetAllCodonPosCounts(): the unique gene names in column "Gene"
# match the genes in gene_names (unique(total_codon_pos_counts$Gene) =
# gene_names) = TRUE.
# gives:
# > str(total_codon_pos_counts)
# Classes "tbl_df", "tbl" and "data.frame": 2749 observations of 3 variables:
#   $ Gene    : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ PosCodon: int  1 2 3 4 5 6 7 8 9 10 ...
#   $ Count   : num  4249 825 1017 1176 1116 ...
#
# TEST::GetAllCodonPosCounts example with tinysim, filter_for_frame =
# FALSE
#
# Gene      PosCodon  Count
# 1 MAT       1         0
# 2 MAT       2         2
# 3 MAT       3         2
# 4 MAT       4         0
# 5 MIKE      1         0
# 6 MIKE      2         1
# 7 MIKE      3         0
# 8 MIKE      4         0
# 9 MIKE      5         0
#
# TEST::GetAllCodonPosCounts example with tinysim, filter_for_frame =
# TRUE
#
# A tibble: 9 x 3
# Gene  PosCodon Count
# <chr>    <int> <dbl>
#   1 MAT          1     0
# 2 MAT          2     0
# 3 MAT          3     2
# 4 MAT          4     0
# 5 MIKE         1     0
# 6 MIKE         2     0
# 7 MIKE         3     0
# 8 MIKE         4     0
# 9 MIKE         5     0
#
# This is expected as only reads mapping to the first nucleotide of
# the codon are retained.

### Functions to add codon identities to positions ###

#' Add codon identities to position and counts.
#'
#' @param gene_names List of all the genes in the sample,
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param asite_disp_path Path to A-site file.
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param min_read_length Minimum read length in H5 file (integer).
#' @param filter_for_frame Method to use for assigning reads (logical).
#' @param filtering_frame Frame from which to select reads from.
#' @param snapdisp Frame to filter to.
#' @param yeast_codon_pos_i200 List of all codons and positions from
#' the sample being studied
#' @return A tibble with the columns `Gene`, "PosCodon", `Count`, and
#' `codon` for a list of genes.
#' @export
#'
#' @examples
#' AddCodonNamesToCodonPosCounts(
#'  gene_names, gff_df, asite_disp_path, dataset, hd_file,
#'  min_read_length, filter_for_frame,filtering_frame,
#'  snapdisp,yeast_codon_pos_i200)
AddCodonNamesToCodonPosCounts <- function(
  gene_names, gff_df, asite_disp_path, dataset, hd_file,
  min_read_length, filter_for_frame, filtering_frame, snapdisp,
  yeast_codon_pos_i200) {

  # Create tibble of read counts and positions.
  total_codon_pos_counts <- GetAllCodonPosCounts(
    gene_names, gff_df, asite_disp_path, dataset, hd_file,
    min_read_length, filter_for_frame, filtering_frame, snapdisp)
  # Add codon identities.
  transcript_tibbles <- left_join(
    total_codon_pos_counts, yeast_codon_pos_i200,
    by = c("PosCodon", "Gene"), keep = FALSE, copy = TRUE)
  # Turn into a tibble.
  transcript_gene_poscodon_frame <- tibble(
    Gene = transcript_tibbles$Gene,
    PosCodon = transcript_tibbles$PosCodon,
    Count = transcript_tibbles$Count,
    Codon = transcript_tibbles$Codon
  )
  return(transcript_gene_poscodon_frame)
}

# TEST::Expect to produce a tibble with each position in the CDS
# having the correct codon beside it.
# TEST::AddCodonNamesToCodonPosCounts(): returns a tibble.
# TEST::AddCodonNamesToCodonPosCounts(): the tibble has 4 columns.
# TEST::AddCodonNamesToCodonPosCounts(): the column names are %in% c("Gene",
# "PosCodon", "Count", "Codon").
# TEST::AddCodonNamesToCodonPosCounts(): number of observations in the output
# tibble = sum of CDS (codon co-ordinates) for all genes in gene_names.
# TEST::AddCodonNamesToCodonPosCounts(): the unique gene names in
# column "Gene" match the genes in gene_names
# (unique(total_codon_pos_counts$Gene) = gene_names) = TRUE.
#
# Example: using Mok-tinysim data, the following tibble is returned
# A tibble: 9 x 4
# Gene  PosCodon Count Codon
# <chr>    <dbl> <dbl> <chr>
#   1 MAT          1     0 ATG
# 2 MAT          2     2 GCC
# 3 MAT          3     2 ACA
# 4 MAT          4     0 TGA
# 5 MIKE         1     0 ATG
# 6 MIKE         2     1 ATC
# 7 MIKE         3     0 AAG
# 8 MIKE         4     0 GAG
# 9 MIKE         5     0 TAA

###  Slice out window around interesting features ###

#' Take an individual gene as an input, then filter for the codon of
#' interest on the gene being investigated.
#'
#' @param gene Gene.
#' @param feature_of_interest Feature of interest.
#' @param transcript_gene_poscodon_frame
#' @return transcript for one gene.
TranscriptForOneGene <- function(
  gene, feature_of_interest, transcript_gene_poscodon_frame) {

  # Filter transcript_gene_poscodon_frame to only include features of
  # interest.
  interesting_feature_tibble <- dplyr::filter(
    transcript_gene_poscodon_frame,
    Codon == feature_of_interest)
  # Select only the gene that the function is being applied to.
  transcript_for_one_gene <- dplyr::filter(
    interesting_feature_tibble, Gene == gene)
  return(transcript_for_one_gene)
}

#' Slice out a window of interest for each position of a feature of
#' interest in a gene (`transcript_gene_poscodon_frame`).
#' Returns an empty tibble if the desired region hangs over the edge
#' of the coding region
#' 
#' @param transcript_for_one_gene Transcript for one gene.
#' @param transcript_gene_poscodon_frame
#' @param gene Gene.
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param expand_width Number of codons to take a slice of either side
#' of occurrences of the feature_of_interest
#' @return feature information.
ExpandRegions <- function(
  transcript_for_one_gene, transcript_gene_poscodon_frame, gene,
  gff_df, dataset, hd_file, expand_width) {

  # Save applied transcript to another object remove problems with
  # purrr::map.
  interesting_features <- transcript_for_one_gene
  # Filter transcript_gene_pos_poscodon_gene_interest to just the gene
  # being used in this purrr::map cycle.
  transcript_gene_pos_poscodon_gene_interest <- dplyr::filter(
    transcript_gene_poscodon_frame, Gene == gene)
  # Get the length of the gene being studied
  gene_length <- filter(gff_df,
                        gff_df$type == "CDS" & gff_df$Name == gene)$width
  # If the position of the feature_of_interest is near the start or
  # stop codon return a NULL in resulting list of tibbles.
  if (interesting_features <= expand_width | interesting_features +
      expand_width > gene_length / 3) {
    return()
  } else {
    # Slice out a window of positions around the feature of interest.
    output_feature_info <- tibble(
      dplyr::slice(transcript_gene_pos_poscodon_gene_interest,
                   (interesting_features - expand_width):(interesting_features + expand_width),
                   each = FALSE),
      Rel_Pos =  seq(- expand_width, expand_width)
    )
    if (dim(output_feature_info)[1] == (2 * expand_width + 1)) {
      return(output_feature_info)
    } else {
      return()
    }
  }
}

#' Select positions of the feature of interest on separate genes.
#'
#' @param gene Gene.
#' @param gene_names List of all the genes in the sample,
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param min_read_length Minimum read length in H5 file (integer).
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param yeast_codon_pos_i200 List of all codons and positions from
#' the sample being studied
#' @param feature_of_interest Feature being studied.
#' @param transcript_gene_poscodon_frame
#' @return feature information.
AllGeneInterestingFeatures <- function(
  gene, gene_names, dataset, hd_file, min_read_length, gff_df,
  yeast_codon_pos_i200, feature_of_interest,
  transcript_gene_poscodon_frame) {

  # Select features of interest occurring on the gene being that the
  # function is being applied to in this purrr::map cycle.
  transcript_for_one_gene <- TranscriptForOneGene(
    gene, feature_of_interest, transcript_gene_poscodon_frame)
  # Apply ExpandRegion() to each occurrence of the feature_of_interest
  # on the gene being studied.
  output_feature_info <- purrr::map(.x = transcript_for_one_gene$PosCodon,
                                    .f = ExpandRegions,
                                    transcript_gene_poscodon_frame,
                                    gene,
				    gff_df,
                                    dataset,
                                    hd_file,
                                    expand_width)
  return(output_feature_info)
}

#' Slice out frames centred on the position of the feature of
#' interest.
#'
#' Given codon identities, positions and counts, slice out a frame
#' centred on the position of the feature of interest. Then, set the
#' position of the feature_of_interest to 0, and change the positions
#' of adjacent codons to be relative to the feature of interest
#'
#' This is done for all occurrences of the feature of interest on all
# 'genes in the sample provided.
#'
#' @param gene_names List of all the genes in the sample,
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param asite_disp_path Path to A-site file.
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param min_read_length Minimum read length in H5 file (integer).
#' @param filter_for_frame Method to use for assigning reads (logical).
#' @param filtering_frame Frame from which to select reads from.
#' @param snapdisp Frame to filter to.
#' @param yeast_codon_pos_i200 List of all codons and positions from
#' the sample being studied.
#' @param feature_of_interest Feature being studied.
#' @param expand_width Number of codons to take a slice of either side
#' of occurrences of the feature_of_interest
#' @return List of tibbles, each containing the columns `Gene`,
#' `PosCodon`, `Count", `Codon`, `Rel_Pos`, for a window of positions
#' determined by `expand_width` around an occurence of the
#' `feature_of_interest`.
#' @export
#'
#' @examples
#' ExpandFeatureRegionAllGenes(gene_names, gff_df, dataset, hd_file,
#' min_read_length, filter_for_frame, filtering_frame, snapdisp,
#' yeast_codon_pos_i200,feature_of_interest,expand_width)
ExpandFeatureRegionAllGenes <- function(
  gene_names, gff_df, asite_disp_path, dataset, hd_file,
  min_read_length,  filter_for_frame, filtering_frame, snapdisp,
  yeast_codon_pos_i200, feature_of_interest, expand_width) {

  # Generate transcript_gene_poscodon_frame.
  transcript_gene_poscodon_frame <- AddCodonNamesToCodonPosCounts(
    gene_names, gff_df, asite_disp_path, dataset, hd_file,
    min_read_length, filter_for_frame, filtering_frame, snapdisp,
    yeast_codon_pos_i200)
  # Apply AllGeneInterestingFeatures to each gene in the sample
  output_feature_info <- purrr::map(.x = gene_names,
                                    .f = AllGeneInterestingFeatures,
                                    gene_names,
                                    dataset,
                                    hd_file,
                                    min_read_length,
                                    gff_df,
                                    yeast_codon_pos_i200,
                                    feature_of_interest,
                                    transcript_gene_poscodon_frame)
  # Givem a list for each gene, containing a list for each occurrence
  # of feature_of_interest, unlist to produce one list, containing
  # each occurrence of feature_of_interest.
  output_feature_info <- unlist(output_feature_info, recursive = F)
  # Remove NULLs, which represent the feature_of_interest occurring
  # within one expand_width of the UTRs.
  output_feature_info <- output_feature_info[
    !sapply(output_feature_info, is.null)]
  return(output_feature_info)
}

# TEST::ExpandFeatureRegionAllGenes(): creates an object of type "list"
# TEST::ExpandFeatureRegionAllGenes(): Returns an empty list if there
# are no occurrences of the feature_of_interest
# TEST::ExpandFeatureRegionAllGenes(): length(list) ==
# nrow(interesting_feature_table)
# TEST::ExpandFeatureRegionAllGenes(): the tibble contains 5 columns =
# TRUE
# TEST::ExpandFeatureRegionAllGenes(): the column names are %in%
#  c("Gene", "PosCodon", "Count", "Codon, "Rel_Pos")
# TEST::ExpandFeatureRegionAllGenes(): number of observations in the output
# tibble = "expand_width" * 2 + 1,
# so if "expand_width" = 5L the number of observations should be 11
# TEST::ExpandFeatureRegionAllGenes(): the position from
# "interesting_feature_positions" has "Rel_Pos" value 0 = TRUE
# TEST::ExpandFeatureRegionAllGenes(): the column "Rel_Pos" goes
# from -"expand_width to +"expand_width"
# Example, tidysim with feature_of_interest == "GCC" and expand_width
# == "1L"
# [[1]]
# A tibble: 3 x 5
# Gene  PosCodon Count Codon Rel_Pos
# <chr>    <dbl> <dbl> <chr>   <int>
# 1 MAT          1     0 ATG        -1
# 2 MAT          2     2 GCC         0
# 3 MAT          3     2 ACA         1

### Functions for normalization ###

#' Change NaNs to zero.
#'
#' @param values Values.
#' @return 0 if NaN else `values`.
SetNaNToZero <- function(values) {
  if (is.nan(values)) {
    values <- 0
  } else {
    values <- values
  }
}

#' Change NaNs to 0s.
#'
#' Having NaNs in any table will cause problems during the overlay
#' step. Changing them to 0 removes this problem, and takes into
#' account that the NaNs are caused by no reads occurring.
CheckForNaN <- function(expanded_features, normalized_expand_tibble) {
  relcount_values <- unlist(normalized_expand_tibble$RelCount)
  relcount_values <- unlist(
    purrr::map(.x = relcount_values, .f = SetNaNToZero))
  dplyr::mutate(expanded_features, RelCount = relcount_values)
}

#' Carry out normalization within each expanded frame.
#' Normalization carried out within each expanded frame so that they
#' are comparable.
#' Normalizes the expand_feature_region list generating a RelCount
#' column with the normalization values
#' Normalizes the ExpandFeatureRegion list generating a RelCount
#' column with the normalization values.
#' As this function is not looped it will only generate one normalized
#' tibble for each occurrence of the feature_of_interest.
#'
#' @param expanded_features List of tidy format data frames (tibbles)
#' generated by the function ExpandFeatureRegion
#' @param expand_width Integer which provides the number of positions
#' on each side of the feature of interest to include in the window.
#' @return List of tibbles which contain the normalized counts within
#' the window so that the feature is comparable despite overall varying
#' levels of expression between genes
#' @export
#' @example
#' ExpandedRegionNormalization(expanded_features, expand_width)
ExpandedRegionNormalization <- function(expanded_features, expand_width) {

  # Normalize the read counts.
  normalized_expand_tibble <- dplyr::mutate(
    expanded_features, RelCount = Count / sum(Count) * (2 * expand_width + 1))
  CheckForNaN(expanded_features, normalized_expand_tibble)
}

# TEST::ExpandedRegionNormalization(): creates a tidy format data frame
# (tibble) = TRUE
# TEST::ExpandedRegionNormalization(): the tibble contains 5 columns =
# TRUE
# TEST::ExpandedRegionNormalization(): the column names are %in%
# c("Gene", "Pos_Codon", "Rel_Count", "Rel_Pos", "RelCount")
# TEST::ExpandedRegionNormalization(): number of observations in the output
# tibble = "expand_width"*2+1,
# if "expand_width" = 5L the number of observations should be 11
# TEST::ExpandedRegionNormalization(): the column "Rel_Pos" goes from
# -"expand_width to +"expand_width"
# TEST::ExpandedRegionNormalization(): sum(normalized_expand_list[[1]]$RelCount)
# /nrow(normalized_expand_list[[1]]) == 1
# TEST::ExpandedRegionNormalizetion(): None of the RelCount columns
# should contain NaN. These should all be set to 0
#
# Example: using tinysim
# A tibble: 3 x 6
# Gene  PosCodon Count Codon Rel_Pos RelCount
# <chr>    <dbl> <dbl> <chr>   <int>    <dbl>
# 1 MAT          1     0 ATG        -1      0
# 2 MAT          2     2 GCC         0      1.5
# 3 MAT          3     2 ACA         1      1.5
#
# After editing tinysim codon position 4 to inculde a GCC, the
# following output_feature_info object is produced
#
# > output_feature_info
# [[1]]
# A tibble: 3 x 5
# Gene  PosCodon Count Codon Rel_Pos
# <chr>    <dbl> <dbl> <chr>   <int>
# 1 MAT          1     0 ATG        -1
# 2 MAT          2     2 GCC         0
# 3 MAT          3     2 ACA         1
#
# [[2]]
# A tibble: 3 x 5
# Gene  PosCodon Count Codon Rel_Pos
# <chr>    <dbl> <dbl> <chr>   <int>
# 1 MIKE         3     0 AAG        -1
# 2 MIKE         4     0 GCC         0
# 3 MIKE         5     0 TAA         1
#
# This is processed by ExpandedRegionNormalization to produce the
# following output:
#
# > normalized_expand_list
# [[1]]
# A tibble: 3 x 6
# Gene  PosCodon Count Codon Rel_Pos RelCount
# <chr>    <dbl> <dbl> <chr>   <int>    <dbl>
# 1 MAT          1     0 ATG        -1      0
# 2 MAT          2     2 GCC         0      1.5
# 3 MAT          3     2 ACA         1      1.5
#
# [[2]]
# A tibble: 3 x 6
# Gene  PosCodon Count Codon Rel_Pos RelCount
# <chr>    <dbl> <dbl> <chr>   <int>    <dbl>
# 1 MIKE         3     0 AAG        -1        0
# 2 MIKE         4     0 GCC         0        0
# 3 MIKE         5     0 TAA         1        0
#
# Normalization carried out for all the tibbles within ExpandList.
# normalized_expand_list <- purrr::map(
#   .x = output_feature_info,
#   .f = ExpandedRegionNormalization,
#   expand_width
# )
#
# TEST:: normalised_expand_list should be of type list
# type(normalized_expand_list)
# [1] "list"
#
# TEST:: Normalized_expand_list should be the same length as
# expand_feature_region
# > length(normalized_expand_list)==length(expand_feature_region)
# [1] TRUE
#
# TEST:: the dimensions of each item in the list shoud be
# [(2*expand_width+1) X 5]
# as there are now 5 rows; Gene, Pos_Codon, Rel_Count, Rel_Pos, RelCount
# > dim(expand_feature_region[[1]])
# [1] 11  5
#
# TEST:: At each position, the sum of RelCount should be equal to
# (2*expand_width+1)
# ie if expand_width was 5:
# sum(normalized_expand_list[[1]]$RelCount)
# [1] 11

### Functions for overlaying the normalized expanded tibbles ###

#' Overlay tidy format data frames (tibbles) to create a single
#' overlayed tibble.
#'
#' FIXME: Takes normalized_expand_list as its input.
#'
#' @param normalized_expand_list Output from the looped function
#' ExpandedRegionNormalization().
#' @param expand_width Integer which provides the number of positions
#' on each side  of the feature of interest to include in the window.
#' @return Tibble which contains the mean counts for each position
#' from the normalized tibbles.
#' @export
#'
#' @example
#' normalized_expand_list <- purrr::map(.x = expand_feature_region,
#'                                      .f = ExpandedRegionNormalization,
#'                                      expand_width)
#' OverlayedTable(normalized_expand_list, expand_width)
OverlayedTable <- function(normalized_expand_list, expand_width) {

  number_of_objects <- length(normalized_expand_list)
  # Reduce normalized_expand_list to the columns Rel_Pos and
  # RelCount.
  result <- lapply(normalized_expand_list, "[", c("Rel_Pos", "RelCount"))
  joined_result <- result %>% purrr::reduce(
    full_join, by = c("Rel_Pos"), sum("RelCount"))
  # Rejoin values.
  joined_rows <- joined_result %>%
    mutate(SumRows = rowSums(
      select(joined_result, -"Rel_Pos")) / number_of_objects)
  # Recreate final tibble
  overlayed_tibbles <- tibble::tibble(
    Rel_Pos = seq(- expand_width, expand_width),
    RelCount = joined_rows$SumRows
  )
}

# TEST::OverlayedTable(): creates a tidy format data frame (tibble) =
# TRUE
# TEST::OverlayedTable(): the tibble contains 2 columns = TRUE
# TEST::OverlayedTable(): the column names are %in% c("Rel_Pos",
# "RelCount")
# TEST::OverlayedTable(): number of observations in the output tibble
# = "expand_width"*2+1,
# if "expand_width" = 5L the number of observations should be 11
# TEST::OverlayedTable(): the column "Rel_Pos" goes from
# -"expand_width to +"expand_width"
# TEST::OverlayedTable(): RelCount is a numeric
# gives:
# > str(overlayed_tibbles)
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 2
# variables
#   $ Rel_Pos : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
#   $ RelCount: num  0.893 1.125 0.992 0.998 0.779 ...

### Run Functions ###

#' Save plot as a PDF.
#'
#' @param overlayed_plot Plot.
#' @param feature_of_interest Feature of interest.
#' @param dataset Dataset.
#' @param output_dir Output directory.
SavePlotPdf <- function(overlayed_plot, feature_of_interest, dataset,
  output_dir) {
  overlayed_plot %>%
    ggsave(
      filename = file.path(output_dir,
                           paste0("Meta_feature_plot_",
                                  feature_of_interest, "_",
                                  dataset,
                                  ".pdf")),
      width = 6,
      height = 5
    )
}

#' Create a tibble of features and Relcounts at position 0 of an
#' overlayed tibble which can be used to compare the RelCount between
#' different features.
#'
#' @param gene_names List of all the genes in the sample,
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param asite_disp_path Path to A-site file.
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param min_read_length Minimum read length in H5 file (integer).
#' @param filter_for_frame Method to use for assigning reads (logical).
#' @param filtering_frame Frame from which to select reads from.
#' @param snapdisp Frame to filter to.
#' @param yeast_codon_pos_i200 List of all codons and positions from
#' the sample being studied
#' @param feature_being_studied List of features being studied.
#' @param expand_width Number of codons to take a slice of either side of.
#' occurrences of the feature_of_interest
#' @return Tibble of features and their RelCounts at position 0.
#' @export
#'
#' @examples
#' feature_rel_use <- purrr::map_df(
#'   .x = feature_of_interest,
#'   .f = FindAllFeatures,
#'   gene_names = gene_names,
#'   gff_df = gff_df,
#'   asite_disp_path = asite_disp_path,
#'   dataset = dataset,
#'   hd_file = hd_file,
#'   min_read_length = min_read_length,
#'   filter_for_frame = filter_for_frame,
#'   filtering_frame = filtering_frame,
#'   snapdisp = snapdisp,
#'   expand_width = expand_width,
#'   yeast_codon_pos_i200 = yeast_codon_pos_i200)
#'
FindAllFeatures <- function(
  gene_names, gff_df, asite_disp_path, dataset, hd_file,
  min_read_length, filter_for_frame, filtering_frame, snapdisp,
 yeast_codon_pos_i200, feature_being_studied, expand_width) {

  # Run ExpandFeatureRegionAllGenes to get a list of occurrences of
  # the feature of interest
  print(paste0("Finding occurences of ", feature_being_studied))
  output_feature_info <- suppressMessages(
    ExpandFeatureRegionAllGenes(
      gene_names, gff_df, asite_disp_path, dataset, hd_file,
      min_read_length, filter_for_frame, filtering_frame, snapdisp,
      yeast_codon_pos_i200, feature_of_interest = feature_being_studied,
      expand_width))
  # Check for the presence of the feature of interest.
  # Output_feature_info being empty will cause problems with
  # normalization.
  if (length(output_feature_info) == 0) {
    print(paste("No occurrences of", feature_being_studied))
    if (expand_width > 1) {
      print("Use an expand_width of 1L to check for occurrences near start or stop codons")
    }
    return()
  }
  # Run ExpandedRegionNormalization to calculate the relative number
  # of reads mapping to each position around the feature of interest.
  normalized_expand_list <- purrr::map(
    .x = output_feature_info,
    .f = ExpandedRegionNormalization,
    expand_width
  )
  # Run OverlayedTable to create an average of reads at positions at
  # and around the feature of interest.
  overlayed_tibbles <- OverlayedTable(normalized_expand_list, expand_width)
  # Create a new tibble listing the feature being studied, and the
  # RelCount at position 0, i.e. RelCount at the feature_of_interest.
  feature_rel_use <- tibble(
    Feature = feature_being_studied,
    RelCount = filter(overlayed_tibbles,
                      overlayed_tibbles$Rel_Pos == 0)$RelCount)
}

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              help = "Path input to h5 file"),
  make_option(c("-d", "--dataset"),
              type = "character",
              help = "Name of the dataset being studied"),
  make_option(c("-g", "--gff"),
              type = "character",
              help = "Path to the GFF3 file of the organism being studied"),
  make_option(c("-a", "--annotation"),
              type = "character",
              help = "Path to codon table for organism"),
  make_option("--feature",
              type = "character",
              help = "Feature of interest, e.g. codon"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "Path to output directory",
              default = "."),
  make_option(c("--expand_width"),
              type = "integer",
              help = "the desired range either side of the feature of
                interest",
              default = 5),
  make_option(c("--frame"),
              type = "integer",
              help = "reading frame to be studied",
              default = 0L),
  make_option(c("--minreadlen"),
              type = "integer",
              help = "minimum read length",
              default = 10),
  make_option(c("--filter_for_frame"),
              type = "logical",
              help = "counts for all reading frames per codon are
                summed and assigned to their corresponding codon.
                You either keep all by not filtering (FALSE) or
                filter for a specific reading frame (TRUE)",
              default = TRUE),
  make_option(c("--snapdisp"),
              type = "integer",
              help = "frame to filter to when using SnapToCodon",
              default = 0L),
  make_option(c("--asite_length"),
              type = "character",
              help = "Path to asite_disp_length. Default is specific
                for yeast when code is run from riboviz directory")
)

opt <- optparse::parse_args(OptionParser(option_list = option_list))

hd_file <- opt$input
dataset <- opt$dataset
gff <- opt$gff
yeast_codon_table <- opt$annotation
feature_of_interest <- opt$feature
output_dir <- opt$output
expand_width <- opt$expand_width
startlen <- opt$startlen
filtering_frame <- opt$frame
min_read_length <- opt$minreadlen
filter_for_frame <- opt$filter_for_frame
snapdisp <- opt$snapdisp
asite_disp_path <- opt$asite_length

print("Starting process")

# If feature_of_interest is a file then load contents. The first
# column of the file should contain the features of interest
# (codons).
if (file.exists(feature_of_interest)) {
  feature_of_interest <- read.csv(feature_of_interest)
  feature_of_interest <- feature_of_interest[, 1]
}

gff_df <- readGFFAsDf(gff)
gene_names <- unique(gff_df$Name)
yeast_codon_pos_i200 <- suppressMessages(
    readr::read_tsv(file = yeast_codon_table))

# feature_of_interest may be provided as a single feature, i.e. a
# codon, or a list of features.
# When only one feature is provided then a graph is plotted based on
# the overlayed relative count around the feature_of_interest.
if (length(feature_of_interest) == 1) {

  # Run ExpandFeatureRegionAllGenes() to get a list of occurrences of
  # feature_of_interest
  print(paste0("Finding occurences of ", feature_of_interest))
  output_feature_info <- suppressMessages(
    ExpandFeatureRegionAllGenes(
      gene_names, gff_df, asite_disp_path, dataset, hd_file,
      min_read_length, filter_for_frame, filtering_frame, snapdisp,
      yeast_codon_pos_i200, feature_of_interest, expand_width))

  # Check for the presence of the
  # feature_of_interest. Output_feature_info being
  # empty will cause problems with normalization
  if (length(output_feature_info) == 0) {
    print("No occurrences of the feature of interest")
    if (expand_width > 1) {
      print("Try script with an expand_width = 1L to check for occurences near to start or stop codon")
    }
    print("Done")
    stop()
  }

  # Run ExpandedRegionNormalization to calculate the relative number
  # of reads mapping to each position around the feature_of_interest
  print("Normalizing read counts")
  normalized_expand_list <- purrr::map(
    .x = output_feature_info,
    .f = ExpandedRegionNormalization,
    expand_width)

  # Run OverlayedTable to create an average of reads at positions at
  # and around the feature_of_interest
  print(paste0(
    "Overlaying tibbles for feature of interest and calculating the average ",
    "relative reads at each position"))
  overlayed_tibbles <- OverlayedTable(normalized_expand_list, expand_width)

  # Create a graph using ggplot
  print("Creating graph")
  overlayed_plot <- ggplot(
    overlayed_tibbles,
    mapping = aes(x = Rel_Pos, y = RelCount)) +
    geom_line() +
    theme_bw() +
    theme(text = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"),
          title = element_text(size = 14, face = "bold")) +
    labs(title = paste0("Relative read counts around feature ",
                        feature_of_interest),
         x = "Position relative to feature of interest",
         y = "Relative read count", size = 2) +
    scale_x_continuous(breaks = seq(-expand_width, expand_width, 2))
  # Save plot as PDF.
  print("Save plot as PDF")
  SavePlotPdf(overlayed_plot, feature_of_interest, dataset, output_dir)
  print("Done")
} else {
  # Use purrr::map to extract the RelCounts at position 0 of all
  # desired features of interest.
  feature_rel_use <- purrr::map_df(
    .x = feature_of_interest,
    .f = FindAllFeatures,
    yeast_codon_pos_i200 = yeast_codon_pos_i200,
    gene_names = gene_names,
    gff_df = gff_df,
    asite_disp_path = asite_disp_path,
    dataset = dataset,
    hd_file = hd_file,
    min_read_length = min_read_length,
    expand_width = expand_width,
    filter_for_frame = filter_for_frame,
    filtering_frame = filtering_frame,
    snapdisp = snapdisp)

  print("Ranking Codons based on RelCount")
  # Rearrange feature_rel_use to be in descending order, so features
  # with the highest relative use are listed at the top
  feature_rel_use <- arrange(feature_rel_use, desc(RelCount))

  # TEST::feature_rel_use should be of class "tbl"
  # TEST::feature_rel_use should have 2 columns.
  # feature_rel_use may not contain as many rows as features were
  # initially input, as if no occurrences of features are found then
  # they won't be included
  # TEST::The headings of feature_rel_use should be "Feature" and "RelCount"
  # EXAMPLE::When run on tinysim, with expand_width = 1,
  # filter_for_frame = False and using all using a list/tsv file of
  # all codons as input for feature_of_interest parameter the
  # following file is produced
  # NOTE::position 4 of MIKE edited to contain GCC to test for
  # occurrence where the feature_of_interest has reads at one position
  # but not the other
  # Feature   RelCount
  # ATC   3.00
  # ACA  1.50
  # GCC   0.75
  # AAG   0.00

  print(head(feature_rel_use, 5))
  print(tail(feature_rel_use, 5))
  print("Saving table as TSV")
  write.table(feature_rel_use,
              file = paste0("Feature_Relative_use_", dataset, ".tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  print("Done")
}
