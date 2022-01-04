# This script uses a H5, GFF and TSV file to create either a
# metafeature plot for a feature of interest, or a table comparing the
# RelCount of multiple features of interest.
# Currently designed for single codons.

## TEST::run on TinySim Dataset



suppressMessages(library(tidyverse))
suppressMessages(library(plotly))
suppressMessages(library(optparse))

if (interactive()) {
  # Use hard-coded script name and assume script is in "rscripts"
  # directory. This assumes that interactive R is being run within
  # the parent of rscripts/ but imposes no other constraints on
  # where rscripts/ or its parents are located.
  this_script <- "codon_feature_analysis.R"
  path_to_this_script <- here::here("rscripts", this_script)
  source(here::here("riboviz","rscripts", "provenance.R"))
  source(here::here("riboviz","rscripts", "read_count_functions.R"))
  source(here::here("riboviz","rscripts", "stats_figs_block_functions.R"))
  source(here::here("riboviz","rscripts", "stats_figs_block_functions.R"))
} else {
  # Deduce file name and path using reflection as before.
  this_script <- getopt::get_Rscript_filename()
  path_to_this_script <- this_script
  source(file.path(dirname(this_script), "provenance.R"))
  source(file.path(dirname(this_script), "read_count_functions.R"))
  source(file.path(dirname(this_script), "stats_figs_block_functions.R"))
}



#' FilterForFrame(): extracts A-site assigned counts, filters for a
#' reading frame for one gene
#'
#' Counts are fetched with GetGeneDatamatrix() and A-site assignment
#' is carried out with CalcAsiteFixed() (in nt positions)
#'
#' The reading frame is then filtered for and the counts are returned
#' (in codon positions)
#'
#' @param gene to get read lengths for.
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml).
#' @param snapdisp integer any additional displacement in the
#' snapping, default = 0L.
#' @param asite_displacement_length integer, lengths used for A-siete assignment.
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data
#' format, as created by readGFFAsDf() from which to extract the UTRs
#' and CDS widths.
#'
#' @return a list of numeric values (read counts) for a reading frame
#' (0/1/2).
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

GetCDSReads <- function(
  gene, left, right, dataset, hd_file, min_read_length,
  asite_displacement_length) {

  # Get the matrix of read counts.
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  # Assign reads to their A site nt, determined by read length.
  reads_asitepos <- CalcAsiteFixed(reads_pos_length,
                                   min_read_length,
                                   asite_displacement_length)
  # Extract the reads that map to the CDS.
  cds_reads <- reads_asitepos[left:right]
  
  return(cds_reads)
}


# TEST: FilterForFrame(): returns a list of numeric values = TRUE
# TEST: FilterForFrame(): length(filtered_counts) = length(cds_length)
# for gene of interest = TRUE
# gives:
# > str(filtered_counts)
# num [1:207] 811 429 488 102 994 146 173 762 13 176 ...

### Fetch and format counts from the h5 file ###

#' Calculate read counts at each codon position for a gene.
#'
#' @param gene Gene in the sample.
#' @param gff_df Data frame version of the GFF3 file contents.
#' @param asite_disp_length integer, lengths used for A-siete assignment.
#' @param dataset Name of dataset in H5 file.
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param min_read_length Minimum read length in H5 file (integer).
#' @param filter_for_frame Method to use for assigning reads (logical).
#' @param filtering_frame Frame from which to select reads from.
#' @param snapdisp Frame to filter to.
#' @return Read counts at each codon position.
GetAllCodonPosCounts1Gene <- function(
  gene,  gff_df, dataset, hd_file, min_read_length, asite_disp_length,
  filter_for_frame = NULL)
{

  # Fetch gff values for a gene, e.g. gene = YAL003W.
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df,
                                         seqnames == gene)

  # Assign start position of the CDS, e.g. for YAL003W: 251.
  left <- as.numeric(subset_gff_df_by_gene %>% 
                      dplyr::select(start))

  # Assign end position of the CDS, e.g. for YAL003W: 871.
  right <- as.numeric(subset_gff_df_by_gene %>% 
                      dplyr::select(end))
  
  cds_reads <- GetCDSReads(gene, left, right, dataset, hd_file, min_read_length,
    asite_disp_length)
  cds_length <- length(cds_reads) / 3
  # Create a tibble, assigning a frame to each nt, so the first nt in
  # each frame has the corresponding frame identity.
  cds_reads_frames <- tibble(Count = cds_reads,
                             Frame = rep(c(0, 1, 2), times = cds_length))
  # If all reads are to be kept.
  if (is.null(filter_for_frame))
  {
    codon_counts_1_gene <- RcppRoll::roll_suml(cds_reads_frames$Count, n=3L, by=3L, fill = NULL)
  } else {
    filtered_for_frame <- dplyr::filter(cds_reads_frames,
                                        Frame == filter_for_frame)
    codon_counts_1_gene<- filtered_for_frame$Count
  }
  # Make a tibble which contains "Gene", "PosCodon" and "Count".
  codon_pos_counts <- tibble(Gene = gene,
                             PosCodon = 1:length(codon_counts_1_gene),
                             Count = codon_counts_1_gene)
  return(codon_pos_counts)
}


#' GetAllCodonPosCounts(): extracts A-site assigned counts for a list
#' of genes
#'
#' This function extracts the A-site assigned counts and generates a
#' tidy data frame (tibble) which contains the counts for all genes in
#' the list of genes.
#'
#' If filter_for_frame = FALSE the GetGeneCodonPosReads1dsnap()
#' function is applied to a list of genes and generates a tibble which
#' contains the counts from all reading frames.
#'
#' If filter_for_frame = TRUE the FilterForFrame() function is applied
#' instead, this generates a tibble which contains the read counts for
#' a reading frame (snapdisp = 0L or 1L or 2L).
#'
#' @param gene from gene_names to get read lengths for.
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml).
#' @param snapdisp integer any additional displacement in the
#' snapping; Default = 0L.
#' @param asite_disp_length integer, lengths used for A-site assignment.
#' @param filter_for_frame TRUE if filtering for a reading frame,
#' FALSE if keeping and grouping all reading frames for each codon.
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data
#' format, as created by readGFFAsDf() from which to extract the UTRs
#' and CDS widths.
#'
#' @return a tibble which contains the columns "Gene", "PosCodon" and
#' "Count" for a list of genes
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
#'                        "data",
#'                        "yeast_standard_asite_disp_length.txt"),
#'                      filter_for_frame = TRUE,
#'                      gff_df)
#'
#' @export
GetAllCodonPosCounts <- function(
  gene_names, dataset, hd_file, min_read_length,
  asite_disp_length, filter_for_frame, gff_df) {
  
  # Apply GetAllCodonPosCounts1Gene() to genes contained within gene_names
  all_codon_pos_counts <- purrr::map_dfr(
    .x = gene_names,
    .f = GetAllCodonPosCounts1Gene,
    dataset = dataset,
    hd_file = hd_file,
    min_read_length = min_read_length,
    asite_disp_length = asite_disp_length,
    filter_for_frame = filter_for_frame,
    gff_df = gff_df)

  return(all_codon_pos_counts)
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
# the codon are retained.ß

#' AddCodonNamesToCodonPosCounts(): takes codon names from the
#' annotation .tsv file, joined to the codon counts table
#'
#' Uses the function GetAllCodonPosCounts().
#'
#' @param .tsv file from which to fetch the codon names associated
#' with the CDS co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_length integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame,
#' FALSE if keeping and grouping all reading frames for each codon
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data
#' format, as created by readGFFAsDf() from which to extract the UTRs
#' and CDS widths.
#'
#' @return a  tidy format data frame (tibble) which contains the genes
#' in gene_names, codon positions, counts and the codon pair.
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
#' AddCodonNamesToCodonPosCounts(codon_pos,
#'                               gene_names,
#'                               dataset = "Mok-simYAL5",
#'                               hd_file = here::here("Mok-simYAL5",
#'                                                    "output",
#'                                                    "A",
#'                                                    "A.h5"),
#'                               min_read_length = 10,
#'                               snapdisp = 0L,
#'                               asite_disp_path = here::here(
#'                                "data",
#'                                "yeast_standard_asite_disp_length.txt"),
#'                               filter_for_frame = TRUE,
#'                               gff_df)
#'
#' @export
AddCodonNamesToCodonPosCounts <- function(
  codon_pos, gene_names, dataset, hd_file, min_read_length,
  asite_disp_length, filter_for_frame, gff_df) {
  
  # Fetch read counts and positions for each gene contained within gene_names
  all_codon_pos_counts <- GetAllCodonPosCounts(
    gene_names, dataset, hd_file, min_read_length = min_read_length,
    filter_for_frame = filter_for_frame, 
    asite_disp_length = asite_disp_length, 
    gff_df = gff_df)

  # Join .tsv file with all_codon_pos_counts
  transcript_tibbles <- all_codon_pos_counts %>% left_join(
    codon_pos,
    by = c("PosCodon" = "CodonPos1", "Gene" = "Gene"),
    keep = FALSE,
    copy = TRUE)

  # Make final tibble with correct column names
  gene_pos_codon_counts <- tibble(Gene = transcript_tibbles$Gene,
                                  CodonPos1 = transcript_tibbles$PosCodon,
                                  CodonPos2 = transcript_tibbles$CodonPos2,
                                  Count = transcript_tibbles$Count,
                                  Codon_1 = transcript_tibbles$Codon_1,
                                  Codon_2 = transcript_tibbles$Codon_2,
                                  CodonPair = transcript_tibbles$CodonPair)

  return(gene_pos_codon_counts)
}
#TEST: AddCodonNamesToCodonPosCounts(): creates a tibble = TRUE
#TEST: AddCodonNamesToCodonPosCounts(): the tibble contains columns =
# TRUE
#TEST: AddCodonNamesToCodonPosCounts(): number of observations in the
# output tibble = sum of CDS (codon co-ordinates) for all genes in
# gene_names.
#TEST: AddCodonNamesToCodonPosCounts(): the column names are %in%
# c("Gene", "CodonPos1", "CodonPos2", "Count", "CodonPair")
#TEST: AddCodonNamesToCodonPosCounts(): the unique gene names in
# column "Gene" match the genes in gene_names
# (unique(all_codon_pos_counts$Gene) = gene_names) = TRUE.
# gives:
# > str(gene_pos_codon_counts)
# Classes "tbl_df", "tbl" and "data.frame":   2,749 observations of 5
# variables:
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  1 2 3 4 5 6 7 8 9 10 ...
#   $ CodonPos2: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ Count    : num  4249 825 1017 1176 1116 ...
#   $ CodonPair: chr  "ATG GCA" "GCA TCC" "TCC ACC" "ACC GAT" ...

### Slice out and expand a frame around interesting features ###
FilterForFeatureOfInterest <- function(
  gene, gene_pos_codon_counts, feature_of_interest) {

  # Check if feature of interest is single codon 
  single <- length(unlist(str_split(feature_of_interest," "))) == 1
  # Filters for the gene of interest from interesting_feature_tibble
  transcript_for_one_gene <- dplyr::filter(gene_pos_codon_counts,
                                           Gene == gene)
  
  # Filter for the occurrences of the feature_of_interest
  if (single)
  {
    transcript_feature_tibble <- dplyr::filter(
      transcript_for_one_gene,
      Codon_1 == feature_of_interest)
  } else {
	  
    transcript_feature_tibble <- dplyr::filter(
	    transcript_for_one_gene,
	    CodonPair == feature_of_interest)
    
  }

  return(transcript_feature_tibble)
}



# Expand the region around an occurrence of the feature_of_interest
ExpandRegionsAroundFeature <- function(
  feature_location, gene_pos_codon_counts, gene, gff_df, expand_width) {

  # filters gene_pos_codon_counts down to the gene being processed
  gene_pos_codon_feature <- dplyr::filter(gene_pos_codon_counts,
                                          Gene == gene)
  
  # fetch gene length for gene being processed from the gff file
  gene_length <- dplyr::filter(gff_df, Name == gene) %>%
                  dplyr::select(width)
  

  # if the window cannot be formed around the feature_of_interest
  if (feature_location <= expand_width  |
      feature_location + expand_width > (gene_length / 3)) {

    return()
  } else {
    # if the window can be formed around the feature_of_interest
    # slice and return the window around the feature_of_interest
    expand_feature_region <- tibble(
      dplyr::slice(gene_pos_codon_feature,
                   (feature_location - expand_width):(feature_location + expand_width),
                   each = FALSE),
                   RelPos =  seq(- expand_width, expand_width))

    if (dim(expand_feature_region)[1] == (2 * expand_width + 1)) {
      return(expand_feature_region)
    } else {
      return()
    }
  }
}

# Takes gene_pos_codon_counts as inputs and select for positions
# on separate genes
InterestingFeaturePerGene <- function(
  gene, feature_of_interest, gene_pos_codon_counts, 
  gff_df, expand_width) {

  features_for_one_gene <- FilterForFeatureOfInterest(gene,
                                                  gene_pos_codon_counts,
                                                  feature_of_interest)
 
  # The if statement ensures that feature positions that are
  # less/more than the expand_width value are discarded
  expand_feature_region <- purrr::map(
    .x = features_for_one_gene$CodonPos1,
    .f = ExpandRegionsAroundFeature,
    gene_pos_codon_counts,
    gene,
    gff_df,
    expand_width)
  return(expand_feature_region)
}

#' ExpandFeatureRegion(): Generates a window around each
#' feature_of_interest across all genes in gene_names
#'
#' Feature_of_interest has position 0, all adjacent codons are
#' assigned relative positions.
#'
#' Uses the function AddCodonNamesToCodonPosCounts() to get counts
#' aligned to positions and codons
#'
#' @param .tsv file from which to fetch the codon names associated
#' with the CDS co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_length integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame,
#' FALSE if keeping and grouping all reading frames for each codon
#' @param feature_of_interest character, each incidence of the feature
#' will be extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions
#' on each side of the feature_of_interest to include in the window
#'
#' @return a list of tibbles, one tibble for each occurrence of the
#' feature of interest with an expanded window around the feature
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
#' ExpandFeatureRegion(codon_pos,
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
#'                       "data",
#'                       "yeast_standard_asite_disp_length.txt"),
#'                     filter_for_frame = TRUE,
#'                     feature_of_interest = "TCC AAG",
#'                     expand_width = 5L
#'                     )
#'
#' @export
ExpandFeatureRegionForAllGenes <- function(
  gene_pos_codon_counts, gene_names, gff_df,
  feature_of_interest, expand_width) {
  # Fetch the assigned read counts for genes in gene_names
  # gene_pos_codon_counts <- AddCodonNamesToCodonPosCounts(
  #   codon_pos, gene_names, dataset, hd_file,
  #   min_read_length, asite_disp_length, filter_for_frame,
  #   gff_df)
 
  expand_feature_region <- purrr::map(
    .x = gene_names,
    .f = InterestingFeaturePerGene,
    feature_of_interest,
    gene_pos_codon_counts,
    gff_df,
    expand_width)

  expand_feature_region <- unlist(expand_feature_region,
                                  recursive = FALSE)
  
  # remove NULLS, which represent features of interest occurring
  # within one expand width of the UTRs
  expand_feature_region <- expand_feature_region[
    !sapply(expand_feature_region, is.null)]
  
  return(expand_feature_region)
}
#TEST: ExpandFeatureRegion(): output is a list of tidy format data
# frames (tibbles) = TRUE. type(expand_feature_region) = "list"
#TEST: ExpandFeatureRegion(): number of tibbles in list matches the
# number of occurrences in "feature_of_interest" list = TRUE
#TEST: ExpandFeatureRegion(): each tibble contains 6 columns = TRUE
#TEST: ExpandFeatureRegion(): the column names are %in% c("Gene",
# "Pos_Codon1", "Pos_Codon2", "Count", "CodonPair", "RelPos")
#TEST: ExpandFeatureRegion(): number of observations in each output
# tibble = "expand_width"*2+1, if "expand_width" = 5L the number of
# observations
# should be 11
#TEST: ExpandFeatureRegion(): the position from
# "interesting_feature_positions"  has "RelPos" value 0 = TRUE
#TEST: ExpandFeatureRegion(): the column "RelPos" goes from
# -"expand_width  to +"expand_width" gives:
# > str(expand_feature_region)
# List of 8 (showing 2)
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 6
# variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  2 3 4 5 6 7 8 9 10 11 ...
#   $ CodonPos2: num  3 4 5 6 7 8 9 10 11 12 ...
#   $ Count    : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ CodonPair: chr  "GCA TCC" "TCC ACC" "ACC GAT" "GAT TTC" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 6
# variables
#   $ Gene     : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1: num  52 53 54 55 56 57 58 59 60 61 ...
#   $ CodonPos2: num  53 54 55 56 57 58 59 60 61 62 ...
#   $ Count    : num  42 53 648 293 121 92 519 79 765 196 ...
#   $ CodonPair: chr  "TTC AAC" "AAC CAC" "CAC ATC" "ATC GCT" ...
#   $ RelPos   : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
# ...

### Normalisation ###

#' Change NaNs to zero.
#'
#' @param values Values.
#' @return 0 if NaN else `values`.
SetNaNToZero <- function(values) {
  if (is.nan(values)) {
    print("NaN present")
    values <- 0
  } else {
    values <- values
  }
}

#' Remove the occurrences of tibbles which contain NaN.
CheckForNaN <- function(normalized_expand_list) {
  Relcount_values <- unlist(normalized_expand_list$RelCount)
  Relcount_values <- unlist(purrr::map(.x = Relcount_values,
                                       .f = SetNaNToZero))
  dplyr::mutate(normalized_expand_list, RelCount = Relcount_values)
}

# function to normalise each tibble before they can be overlayed
Normalization <- function(expand_feature_region, expand_width) {
  dplyr::mutate(expand_feature_region,
                RelCount = Count / sum(Count) * (2 * expand_width + 1))
}

#' ExpandedRegionNormalisation(): carries out normalization within
#' each expanded frame so that they are comparable
#'
#' Normalizes the ExpandFeatureRegion() list, generates a RelCount
#' column with the normalization values
#'
#' @param .tsv file from which to fetch the codon names associated
#' with the CDS co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param gff_df data.frame or tibble; riboviz-format GFF in tidy data
#' format, as created by readGFFAsDf() from which to extract the UTRs
#' and CDS widths.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_length integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame,
#' FALSE if keeping and grouping all reading frames for each codon.
#' @param feature_of_interest character, each incidence of the feature
#' will be extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions
#' on each side of the feature_of_interest to include in the window
#'
#' @return a list of normalized tibbles, where each tibble is an
#' occurrence of the feature_of_interest
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
#' ExpandedRegionNormalisation(codon_pos,
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
#'                               "data",
#'                               "yeast_standard_asite_disp_length.txt"),
#'                             filter_for_frame = TRUE,
#'                             feature_of_interest = "TCC AAG",
#'                             expand_width = 5L)
#'
#' @export
ExpandedRegionNormalisation <- function(
  gene_pos_codon_counts, gene_names, gff_df,
  feature_of_interest, expand_width) {

  # fetch the expanded tibbles for each occurrence of the
  # feature_of_interest
  expand_feature_region <- ExpandFeatureRegionForAllGenes(
    gene_pos_codon_counts, gene_names, gff_df,
    feature_of_interest, expand_width)
  
  if (length(expand_feature_region) == 0) {
    print(paste("No occurrences of", feature_of_interest))
    if (expand_width > 1) {
      print("Use expand_width = 1L to check for occurrences near to start or stop codon")
    }
    return()
  }
  
  normalized_expand_list <- purrr::map(.x = expand_feature_region,
                                       .f = Normalization,
                                       expand_width)

  normalized_expand_list <- purrr::map(.x = normalized_expand_list,
                                       .f = CheckForNaN)
  return(normalized_expand_list)
}
#TEST: ExpandedRegionNormalisation(): creates a list of tidy format
#  data frame (tibble) (type(normalized_expand_list)) = TRUE
#TEST: ExpandedRegionNormalisation(): the tibble contains 6 columns =
# TRUE
#TEST: ExpandedRegionNormalisation(): the column names are %in%
# c("Gene", "CodonPos1", "CodonPos2", "Count", "RelPos", "RelCount")
#TEST: ExpandedRegionNormalisation(): number of observations in the
# output tibble = "expand_width"*2+1, if "expand_width" = 5L the
# number of observations should be 11
#TEST: ExpandedRegionNormalisation(): the column "RelPos" goes from
# -"expand_width to +"expand_width" gives:
# > str(normalized_expanded_feature_region)
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 6
# variables
#   $ Gene      : chr  "YAL003W" "YAL003W" "YAL003W" "YAL003W" ...
#   $ CodonPos1 : num  2 3 4 5 6 7 8 9 10 11 ...
#   $ CodonPos2 : num  3 4 5 6 7 8 9 10 11 12 ...
#   $ Count     : num  429 488 102 994 146 173 762 13 176 98 ...
#   $ RelPos    : int  -5 -4 -3 -2 -1 0 1 2 3 4 ...
#   $ RelCount  : num  1.347 1.532 0.32 3.12 0.458 ...


#' FindAllFeatures(): extracts the RelCount values for multiple
#' features of interest
#'
#' The feature contains the functions ExpandFeatureRegionAllGenes(),
#' ExpandedRegionNormalization() and OverlayedTibble(), which are
#' defined above.
#'
#' @param .tsv file from which to fetch the codon names associated
#' with the CDS co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param asite_disp_length integer, lengths used for A-site assignment
#' @param expand_width integer which provides the number of positions
#' on each side of the feature_of_interest to include in the window
#' @param filter_for_frame logical, TRUE if filtering for a reading
#' frame, FALSE if keeping and grouping all reading frames for each
#' codon
#' @param feature_of_interest character, each incidence of the feature
#' will be extracted from transcript_info_tibble
#' @param snapdisp integer any additional displacement in the snapping
#'
#' @return a tidy format data frame (tibble) consisting of the columns
#' '"Feature" and "RelCount" containing all objects listed in
#' feature_of_interest and their relative counts
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
#' FindAllFeatures(codon_pos,
#'                 gene_names,
#'                 dataset = "Mok-simYAL5",
#'                 hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"),
#'                 gff_df,
#'                 min_read_length = 10,
#'                 snapdisp = 0L,
#'                 asite_disp_path = here::here(
#'                     "data",
#'                     "yeast_standard_asite_disp_length.txt"),
#'                 filter_for_frame = TRUE,
#'                 feature_of_interest = "TCC AAG",
#'                 expand_width = 5L)
#'
#'@export
FindFeature <- function(
  feature_of_interest, gene_pos_codon_counts, gene_names,
  gff_df, expand_width) {
  
  # fetch a list of normalised tibbles consisting of each occurrence of the
  # feature_of_interest
  normalized_expand_list <- ExpandedRegionNormalisation(
    gene_pos_codon_counts, gene_names, gff_df,
    feature_of_interest, expand_width)
  
  if (length(normalized_expand_list) == 0) {
    return()
  }
  
  # the number of objects inside normalized_expand_list
  number_of_objects <- length(normalized_expand_list)
  
  # reduces normalized_expand_list to the columns RelPos and RelCount
  result <- lapply(normalized_expand_list, "[", c("RelPos", "RelCount"))
  
  # sums the columns RelPos and RelCount for all the tibbles into one
  # tibble
  joined_result <- purrr::reduce(.x = result, .f = `+`)
  
  # average the tibble by dividing by the number of objects summed
  joined_rows <- mutate(joined_result / number_of_objects)
  
  # create the final overlayed tibble with the relative positions and
  # the normalised and overlayed relative counts around the
  # feature_of_interest
  overlayed_tibbles <- tibble::tibble(
    RelPos = seq(-expand_width, expand_width),
    RelCount = joined_rows$RelCount)


  if (length(overlayed_tibbles) == 0) {
    return()
  }
  return(overlayed_tibbles)
}
#TEST: FindAllOccurences(): creates a tibble = TRUE
#TEST: FindAllOccurences(): the tibble contains 2 columns = TRUE
#TEST: FindAllOccurences(): the column names are %in% c("Feature",
# "RelCount")
#TEST: FindAllOccurences(): number of observations in the output tibble:
#      length(feature_of_interest) == length(feature_rel_use)
# gives:
# > str(feature_rel_use)
# Classes "tbl_df", "tbl" and "data.frame":   11 observations of 2
# variables
# $ Feature : chr [1:11] "CCA TGG" "AGA TGG" "GTA GTG" "TCA TAC" ...
# $ RelCount: num [1:11] 3.41 3.35 2.19 2.12 2.11 ...




#' Create metafeature plot for feature of interest.
#'
#' @param hd_file Path to H5 file holding read data for all genes.
#' @param dataset Name of dataset in H5 file.
#' @param gff Path to the GFF3 file of the organism being studied"),
#' @param codon_pos_table Path to codon positions table.
#' @param feature_of_interest Feature of interest e.g. codon pair.
#' @param output_dir Output directory.
#' @param expand_width Number of codons to take a slice of either side
#' of occurrences of the feature of interest.
#' @param filter_for_frame Frame from which to select reads from. If NULL, do not filter.
#' @param min_read_length Minimum read length in H5 file (integer).
# @param snapdisp Frame to filter to.
#' @param asite_disp_path Path to A-site file.
CodonFeatureAnalysis <- function(hd_file, dataset, gff, codon_pos_table,
  feature_of_interest, output_dir = ".", expand_width = 5,
  filter_for_frame = NULL, min_read_length = 10, asite_disp_path) {

  # If feature_of_interest is a file then load contents. The first
  # column of the file should contain the features of interest
  # (codons).
  if (file.exists(feature_of_interest)) {
    feature_of_interest <- suppressMessages(read_tsv(feature_of_interest,col_types=cols()))
    feature_of_interest <- unname(unlist(feature_of_interest[, 1]))
    names(feature_of_interest) <- feature_of_interest
  } 

  if (file.exists(asite_disp_path))
  {
    asite_displacement_length <- suppressMessages(ReadAsiteDisplacementLengthFromFile(asite_disp_path))
  } else{
    stop("ERROR: A-site offset file not found.")
  }

  gff_df <- readGFFAsDf(gff)

  gff_df <- gff_df %>% dplyr::filter(type == "CDS")

  gene_names <- rhdf5::h5ls(hd_file,
                            recursive = 1)$name

  codon_pos <- suppressMessages(read_tsv(file = codon_pos_table))

  ## This operation is relatively quick. Can make even if not doing codon pairs analysis
  codon_pos <- tibble::tibble(
    Gene = codon_pos$Gene,
    CodonPos1 = codon_pos$PosCodon,
    CodonPos2 = dplyr::lead(codon_pos$PosCodon),
    Codon_1 = codon_pos$Codon,
    Codon_2 = (dplyr::lead(codon_pos$Codon)
                       %>% str_replace_all("ATG", "NA"))
  )
  codon_pos <- codon_pos %>% mutate(CodonPair=paste(Codon_1,Codon_2))

  ## Create table with codon counts per codon per gene
  gene_pos_codon_counts <- AddCodonNamesToCodonPosCounts(
    codon_pos, gene_names, dataset, hd_file,
    min_read_length, asite_disp_length, filter_for_frame,
    gff_df)
  
  ## Map over different features of interest
  feature_rel_use <- purrr::map_df(
    .x = feature_of_interest,
    .f = FindFeature,
    gene_pos_codon_counts = gene_pos_codon_counts,
    gene_names = gene_names,
    gff_df = gff_df,
    expand_width = expand_width,
    .id="Codon")

  return(feature_rel_use)

}

FilterFeatureTableByRelativePosition <-function(feature_table,relative_pos = 0)
{
  # Create a new tibble listing the feature being studied, and the
  # RelCount at position 0, ie RelCount at the feature_of_interest
  feature_table_filter <- feature_table %>% filter(RelPos == 0) %>% dplyr::select(Codon,RelCount)
  feature_table_filter <- arrange(feature_table_filter, desc(RelCount))
  return(feature_table_filter)
}


WriteCodonFeatureAnalysis <-function(output_dir,feature_table)
{
  # Save file
  tsv_file_path <- file.path(output_dir, "codon_feature_analysis.tsv")
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    feature_table,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
}



option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              help = "Path input to h5 file"),
  make_option(c("-d", "--dataset"),
              type = "character",
              help = "Name of dataset being studied"),
  make_option(c("-g", "--gff"),
              type = "character",
              help = "Path to the GFF3 file"),
  make_option(c("-a", "--annotation"),
              type = "character",
              help = "Path to codon positions table"),
  make_option(c("--feature"),
              type = "character",
              help = "Feature of interest, e.g. codon pair"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "Path to output directory",
              default = "."),
  make_option(c("--expand_width"),
              type = "integer",
              help = "Range around feature_of_interest",
              default = 5),
  make_option(c("--minreadlen"),
              type = "integer",
              help = "Minimum read length",
              default = 10),
  make_option(c("--filter_for_frame"),
              type = "integer",
              help = "Keep all or filter for a reading frame",
              default = NULL),
  make_option(c("--asite_length"),
              type = "character",
              help = "Path to asite_disp_length")
)

# opt <- optparse::parse_args(OptionParser(option_list = option_list))
# 
# hd_file <- opt$input
# dataset <- opt$dataset
# gff <- opt$gff
# codon_pos_table <- opt$annotation
# feature_of_interest <- opt$feature
# output_dir <- opt$output
# expand_width <- opt$expand_width
# filter_for_frame <- opt$filter_for_frame
# min_read_length <- opt$minreadlen
# asite_disp_path <- opt$asite_length

hd_file <- "data/Mok-simYAL5/A.h5"
dataset <- "Mok-simYAL5"
gff <- "data/Mok-simYAL5/Scer_YAL_5genes_w_250utrs.gff3"
asite_disp_path <- "data/yeast_standard_asite_disp_length.txt"
codon_pos_table <- "data/yeast_codon_table.tsv"
feature_of_interest <- "data/codon-pairs.tsv"
output_dir <- "."
output_file <- "Feature_Relative_use_Mok-simYAL5.tsv"
min_read_length <- 10
expand_width <- 5
filter_for_frame <- 0

feature_rel_use <- CodonFeatureAnalysis(hd_file, dataset, gff, codon_pos_table,
  feature_of_interest, output_dir, expand_width, filter_for_frame,
  min_read_length, asite_disp_path)


feature_rel_use_pos0 <- FilterFeatureTableByRelativePosition(feature_rel_use)
#WriteCodonFeatureAnalysis(output_dir,feature_rel_use_pos0)

