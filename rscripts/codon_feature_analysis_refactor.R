# This script uses a H5, GFF and TSV file to create either a
# metafeature plot for a feature of interest, or a table comparing the
# RelCount of multiple features of interest.
# Currently designed for single codons.

## TEST::run on TinySim Dataset

# Example: 
# Rscript rscripts/codon_feature_analysis_refactor.R --filter_for_frame 0 -i data/Mok-simYAL5/A.h5 -d Mok-simYAL5 -g data/Mok-simYAL5/Scer_YAL_5genes_w_250utrs.gff3 -a data/yeast_codon_table.tsv --asite_length data/yeast_standard_asite_disp_length.txt -o .
# This assumes that example-datasets is in the same parent directory as riboviz


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
  source(here::here("rscripts", "provenance.R"))
  source(here::here("rscripts", "read_count_functions.R"))
  source(here::here("rscripts", "stats_figs_block_functions.R"))
} else {
  # Deduce file name and path using reflection as before.
  this_script <- getopt::get_Rscript_filename()
  path_to_this_script <- this_script
  source(file.path(dirname(this_script), "provenance.R"))
  source(file.path(dirname(this_script), "read_count_functions.R"))
  source(file.path(dirname(this_script), "stats_figs_block_functions.R"))
}



#' GetCDSReads(): extracts A-site assigned counts, filters for a
#' reading frame for one gene
#'
#' Counts are fetched with GetGeneDatamatrix() and A-site assignment
#' is carried out with CalcAsiteFixed() (in nt positions)
#'
#' The reading frame is then filtered for and the counts are returned
#' (in codon positions)
#'
#' @param gene to get read lengths for.
#' @param start Data frame version of the GFF3 file contents.
#' @param end Data frame version of the GFF3 file contents.
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml).
#' @param asite_displacement_length integer, lengths used for A-siete assignment.
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
#' GetCDSReads(gene = "YAL003W",
#'                dataset = "Mok-simYAL5",
#'                hd_file,
#'                min_read_length = 10,
#'                snapdisp = 0L,
#'                asite_disp_path = here::here(
#'                  "data", "yeast_standard_asite_disp_length.txt"))
#'
#' @export

GetCDSReads <- function(
  gene, start, end, dataset, hd_file, min_read_length,
  asite_displacement_length) {

  # Get the matrix of read counts.
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  # Assign reads to their A site nt, determined by read length.
  reads_asitepos <- CalcAsiteFixed(reads_pos_length,
                                   min_read_length,
                                   asite_displacement_length)
  # Extract the reads that map to the CDS.
  cds_reads <- reads_asitepos[start:end]
  
  return(cds_reads)
}


# TEST: GetCDSReads(): returns a list of numeric values = TRUE
# TEST: GetCDSReads(): length(filtered_counts) = length(cds_length)
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
#' @param filter_for_frame integer (0, 1, or 2), reading frame to use. If NULL, 
#' data from all reading frames are combined.
#' @param feature_type feature type to take from gff, default "CDS".
#' @return Read counts at each codon position.
GetAllCodonPosCounts1Gene <- function(
  gene,  gff_df, dataset, hd_file, min_read_length, asite_disp_length,
  filter_for_frame = NULL,
  feature_type = "CDS")
{

  # Fetch gff values for a gene, e.g. gene = YAL003W
  # Only take the CDS feature
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df,
                                         seqnames == gene,
                                         type == feature_type)

  # Assign start position of the CDS, e.g. for YAL003W: 251.
  start <- as.integer(subset_gff_df_by_gene %>% 
                      dplyr::select(start))

  # Assign end position of the CDS, e.g. for YAL003W: 871.
  end <- as.integer(subset_gff_df_by_gene %>% 
                      dplyr::select(end))
  
  cds_reads <- GetCDSReads(gene, start, end, dataset, hd_file, min_read_length,
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
#' If filter_for_frame = NULL, this generates a tibble is created which 
#' contains the counts from all reading frames.
#'
#' If filter_for_frame = 0, 1, or 2, this generates a tibble which contains 
#' the read counts for a reading frame.
#'
#' @param gene_names list of genes.
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param min_read_length integer, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml).
#' @param asite_disp_length integer, lengths used for A-site assignment.
#' @param filter_for_frame integer (0,1,or 2), reading frame to use. If NULL, 
#' data from all reading frames are combined. 
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
#'                      asite_disp_path = here::here(
#'                        "data",
#'                        "yeast_standard_asite_disp_length.txt"),
#'                      filter_for_frame = 0,
#'                      gff_df)
#'
#' @export
GetAllCodonPosCounts <- function(
  gene_names, dataset, hd_file, min_read_length,
  asite_disp_length, filter_for_frame, gff_df) {
  
  # Apply GetAllCodonPosCounts1Gene() to genes contained within gene_names
  all_codon_pos_counts <- 
    purrr::map_dfr(
    .x = gene_names,
    .f = GetAllCodonPosCounts1Gene,
    dataset = dataset,
    hd_file = hd_file,
    min_read_length = min_read_length,
    asite_disp_length = asite_disp_length,
    filter_for_frame = filter_for_frame,
    gff_df = gff_df) %>%
    dplyr::mutate(Gene = factor(Gene, levels = gene_names))

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
# the codon are retained.

WriteFeatureAnalysis <-function(feature_table,
                                filename = "feature_relpos_relcount.tsv",
                                output_dir = ".")
{
  # Save file
  tsv_file_path <- file.path(output_dir, filename)
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


### Generate plot around feature_of_interest ###

#' GeneratePositionFeaturePlot(): Generates a metafeature plot around the
#' feature_of_interest
#'
#' @param feature_rel_use tibble with columns:
#' "RelPos" (integer, position relative to feature of interest, which is at position 0)
#' "RelCount (double, relative counts at each relative position)
#' @param expand_width integer which provides the number of positions
#' on each side of the feature_of_interest to include in the window
#' 
#' @return A plot which shows the ribosomal occupancy around a feature
#' of interest
#'
#' @example
#'
#' feature_rel_use <- tibble(RelPos = seq(-5L,5L), RelCount = dnorm(RelPos))
#' GeneratePositionFeaturePlot(feature_rel_use,
#'              expand_width = 5L)
#'
#' @export
GenerateRelPosRelCountPlot <- function(feature_rel_use, expand_width) {
  ggplot(
    feature_rel_use,
    mapping = aes(x = RelPos,
                  y = RelCount)) +
    geom_line() +
    theme(text = element_text(size = 12),
          axis.title = element_text(size = 10, face = "bold"),
          title = element_text(size = 12, face = "bold")) +
    labs(x = "Relative Position",
         y = "Normalised occupancy", size = 2) +
    scale_x_continuous(breaks = seq(-expand_width, expand_width, 1)) +
    theme_bw() +
    theme(panel.grid = element_blank(), strip.background = element_blank())
}
# TEST: GenerateRelPosRelCountPlot(): produces a single plot =TRUE

CalculateCodonPairTable1Gene <- function(codon_poscodon_1gene) {
  # Currently this mostly works but
  assertthat::assert_that(
    assertthat::are_equal(codon_poscodon_1gene$PosCodon,
                          seq(min(codon_poscodon_1gene$PosCodon),
                              max(codon_poscodon_1gene$PosCodon))),
    msg = "codon_poscodon_1gene does not have PosCodons in integer sequence")
  
  codon_poscodon_1gene %>%
    dplyr::mutate(CodonNext = dplyr::lead(codon_poscodon_1gene$Codon),
                  CodonPair = paste(Codon, CodonNext, sep = "_") %>%
                    replace(nrow(codon_poscodon_1gene), NA))
}

CalculateCodonPairTable <- function(features_gene_poscodon) {
  assertthat::assert_that(
    assertthat::has_name(features_gene_poscodon, 
                         c("Gene", "PosCodon", "Codon")))
  features_gene_poscodon %>%
    dplyr::nest_by(Gene) %>%
    dplyr::summarise(CalculateCodonPairTable1Gene(data),
                     .groups = "drop")
}

#' Expand the window around a gene, poscodon of interest
#' 
#' Return an empty output if the window extends beyond either end of the gene.
#' 
ExpandWindowGenePosCodon <- function(
  gene, poscodon, gene_pos_codon_counts, expand_width) {
  
  gene_poscodon_window <- 
    dplyr::filter(gene_pos_codon_counts,
                  Gene == gene,
                  PosCodon >= poscodon - expand_width,
                  PosCodon <= poscodon + expand_width)
  
  if (nrow(gene_poscodon_window) == (2 * expand_width + 1)) {
    assertthat::assert_that(
      assertthat::are_equal(gene_poscodon_window$PosCodon,
                            seq(poscodon - expand_width,
                                poscodon + expand_width)),
      msg = "gene_pos_codon_feature does not have PosCodons in integer sequence"
    )
    return(gene_poscodon_window %>%
             dplyr::mutate(RelPos =  seq(-expand_width, expand_width))
    )
  } else {
    # window extends beyond an end of gene
    return(tibble())
  }
}

# filter for windows with minimum counts (to avoid NaNs)
# and normalize divide by mean
NormalizeFilterWindow <- function(pos_window,
                                  min_count = 1,
                                  na.rm = TRUE) {
  if(nrow(pos_window) == 0) {
    return(tibble()) 
  } else if(sum(pos_window$Count, na.rm = na.rm) <= min_count) {
    return(tibble()) 
  } else {
    return(dplyr::mutate(pos_window,
                         RelCount = Count / mean(Count, na.rm = na.rm))
    )
  }
}

#' First ExpandWindowGenePosCodon, then NormalizeFilterWindow
#'
ExpandWindowNormGenePosCodon <- 
  function(gene, poscodon, 
           gene_pos_codon_counts, expand_width,
           min_count = 1, na.rm = TRUE,
           select_relposonly = FALSE) {
    gene_poscodon_window <- 
      ExpandWindowGenePosCodon(gene = gene, 
                               poscodon = poscodon, 
                               gene_pos_codon_counts = gene_pos_codon_counts, 
                               expand_width = expand_width) %>%
      NormalizeFilterWindow(min_count = min_count, 
                            na.rm = na.rm)
    
    if(select_relposonly & nrow(gene_poscodon_window) > 0 & ncol(gene_poscodon_window) > 0) {
      # This removes Gene, PosCodon from output while keeping RelPos
      # It's designed to remove problems overwriting columns when running 
      # group_by(Gene, PosCodon) %>% summarise(ExpandWindowNormGenePosCodon,...)
      # But I don't know if this is the best way to do that.
      gene_poscodon_window <- 
        dplyr::select(gene_poscodon_window, -Gene, -PosCodon)
    }
    return(gene_poscodon_window)
}

#' For all features in a table ExpandWindowNormGenePosCodon
#' 
#' Uses rowwise() to summarise for each row (Gene, PosCodon)
#' 
#' TO DO: need additional function to count how many instances of the windows there are
#' 
ExpandWindowsNormGenePosCodon <-
  function(features_gene_poscodon,
           gene_pos_codon_counts, expand_width,
           min_count = 1, na.rm = TRUE) {
    features_gene_poscodon %>%
      rowwise() %>%
      dplyr::summarise(
        ExpandWindowNormGenePosCodon(
          gene = Gene, 
          poscodon = PosCodon, 
          gene_pos_codon_counts = gene_pos_codon_counts, 
          expand_width = expand_width,
          select_relposonly = TRUE,
          min_count = min_count, 
          na.rm = na.rm),
        .groups = "drop")
  }


CalculateSummaryCountByRelPos <- function(windows_counts_relpos, na.rm = TRUE) {
  if( nrow(windows_counts_relpos) > 0 & ncol(windows_counts_relpos) > 0) {
    windows_counts_relpos %>%
      dplyr::group_by(RelPos) %>%
      dplyr::summarise(RelCount = mean(RelCount, na.rm = na.rm),
                       TotCount = sum(Count, na.rm = na.rm),
                       .groups = "drop") %>%
      return() 
  } else {
    return(tibble())
  }
}

SummariseCountsByRelPos <- 
  function(features_gene_poscodon,
           gene_pos_codon_counts, expand_width,
           min_count = 1, na.rm = TRUE) {
    features_gene_poscodon %>%
      ExpandWindowsNormGenePosCodon(gene_pos_codon_counts = gene_pos_codon_counts, 
                                    expand_width = expand_width,
                                    min_count = min_count, 
                                    na.rm = na.rm) %>%
      CalculateSummaryCountByRelPos(na.rm = na.rm)
  }

#' Summarise counts by feature and relative position
#' 
SummariseCountsByFeatureAndRelPos <- 
  function(features_gene_poscodon,
           feature_var,
           gene_pos_codon_counts, expand_width,
           min_count = 1, na.rm = TRUE) {
    
    assertthat::assert_that(
      assertthat::has_name(features_gene_poscodon, 
                           feature_var))
    
    features_gene_poscodon %>%
      dplyr::nest_by(.data[[feature_var]]) %>%
      dplyr::summarise(
        SummariseCountsByRelPos(features_gene_poscodon = data,
                                gene_pos_codon_counts = gene_pos_codon_counts, 
                                expand_width = expand_width,
                                min_count = min_count, 
                                na.rm = na.rm),
        .groups = "keep")
  }

##### Command-line options and core script

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
              help = "Path to asite_disp_length"),
  make_option(c("--is_codon_feature"),
              type = "logical",
              help = "If true (default), features of interest are individual codons or codon pairs. If false, features of interest
              are specified as TSV with columns Gene, CodonPos, and Feature.",
              default = TRUE)
)

opt <- optparse::parse_args(OptionParser(option_list = option_list))

hd_file <- opt$input
dataset <- opt$dataset
gff <- opt$gff
codon_pos_table <- opt$annotation
feature_of_interest <- opt$feature
output_dir <- opt$output
expand_width <- opt$expand_width
filter_for_frame <- opt$filter_for_frame
min_read_length <- opt$minreadlen
asite_disp_path <- opt$asite_length
is_codon_feature_analysis <- opt$is_codon_feature

# This block of code gives example options for testing
# hd_file <- "data/Mok-simYAL5/A.h5"
# dataset <- "Mok-simYAL5"
# gff <- "data/Mok-simYAL5/Scer_YAL_5genes_w_250utrs.gff3"
# asite_disp_path <- "data/yeast_standard_asite_disp_length.txt"
# codon_pos_table <- "data/yeast_codon_table.tsv"
# feature_of_interest <- "data/codon-pairs.tsv"
# output_dir <- "."
# output_file <- "Feature_Relative_use_Mok-simYAL5.tsv"
# min_read_length <- 10
# expand_width <- 5
# filter_for_frame <- 0
# is_codon_feature_analysis <- TRUE

# Load gff as a data frame
gff_df <- readGFFAsDf(gff)

gene_names <- levels(gff_df$seqnames)

# Load codon position table
gene_pos_codon_table <- 
  readr::read_tsv(codon_pos_table, comment = "#",
                  col_types = "cic")


# Get the counts for each gene and position
all_codon_pos_counts <- GetAllCodonPosCounts(
  gene_names,
  dataset, 
  hd_file = hd_file, 
  min_read_length = min_read_length,
  filter_for_frame = filter_for_frame, 
  asite_disp_length = readr::read_tsv(asite_disp_path, 
                                      comment = "#", 
                                      show_col_types = FALSE), 
  gff_df = gff_df)


# Run for all codons
codon_rel_use <-
  gene_pos_codon_table %>%
  dplyr::filter(Gene %in% gene_names) %>%
  SummariseCountsByFeatureAndRelPos(
    feature_var = "Codon",
    gene_pos_codon_counts = all_codon_pos_counts,
    expand_width = 5,
    min_count = 1, na.rm = TRUE)

# Write the output 
WriteFeatureAnalysis(codon_rel_use, filename = "codon_relpos_relcount.tsv", output_dir = output_dir)

# Generate and save faceted plot
codon_rel_use_plot <- GenerateRelPosRelCountPlot(codon_rel_use, expand_width = 5L) + 
  facet_wrap(~Codon)

ggsave(codon_rel_use_plot,
       filename = file.path(output_dir,
                            paste0("codon_rel_use_plot_faceted.pdf")),
       width = 8, height = 8)
