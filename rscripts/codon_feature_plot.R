### Generate plot around feature_of_interest ###

#' GeneratePlot(): Generates a metafeature plot around the
#' feature_of_interest
#'
#' Fetches the overlayed tibble using the function overlayed_tibbles()
#'
#' @param .tsv file from which to fetch the codon names associated
#' with the CD co-ordinates for each gene
#' @param gene from gene_names
#' @param dataset name of dataset stored in .h5 file.
#' @param hd_file name of .h5 hdf5 file holding read data for all
#' genes, created from BAM files for dataset samples.
#' @param gff_df from which to extract the UTRs and CDS widths.
#' @param min_read_length numeric, minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param snapdisp integer any additional displacement in the snapping
#' @param asite_disp_path integer, lengths used for A-site assignment
#' @param filter_for_frame TRUE if filtering for a reading frame,
#' FALSE if keeping and grouping all reading frames for each codon
#' @param feature_of_interest character, each incidence of the feature
#' will be extracted from transcript_info_tibble
#' @param expand_width integer which provides the number of positions
#' on each side of the feature_of_interest to include in the window
#' @param size the size of the text on the metafeature plot
#'
#' @return A plot which shows the ribosomal occupancy around a feature
#' of interest
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
#' GeneratePlot(codon_pos_pair_i200,
#'              gene_names,
#'              dataset = "Mok-simYAL5",
#'              hd_file = here::here("Mok-simYAL5", "output", "A", "A.h5"),
#'              gff_df,
#'              min_read_length = 10,
#'              snapdisp = 0L,
#'              asite_disp_path = here::here(
#'                "data",
#'                "yeast_standard_asite_disp_length.txt"),
#'              filter_for_frame = TRUE,
#'              feature_of_interest = "TCC AAG",
#'              expand_width = 5L,
#'              size = 12)
#'
#' @export
GeneratePlot <- function(
  codon_pos_pair_i200, gene_names, dataset, hd_file, gff_df,
  min_read_length, snapdisp, asite_disp_path, filter_for_frame,
  feature_of_interest, expand_width, size = 12) {

  overlayed_tibbles <- OverlayedTibble(
    codon_pos_pair_i200, gene_names, dataset, hd_file, gff_df,
    min_read_length, snapdisp, asite_disp_path, filter_for_frame,
    feature_of_interest, expand_width)

  overlayed_plot <- ggplot(
    overlayed_tibbles,
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
# <feature_of_interest>
#TEST: GeneratePlot(): the x-axis is "Distance from codon pair (3nt)
# and y-axis is "Normalised ribosomal occupancy"

#' Save plot as a PDF.
#'
#' @param overlayed_plot Plot.
#' @param feature_of_interest Feature of interest.
#' @param output_dir Output directory.
SavePlotPDF <- function(
  overlayed_plot, feature_of_interest, output_dir) {
  overlayed_plot %>%
    ggsave(
      filename = file.path(output_dir,
                           paste0("Meta_feature_plot",
                                  feature_of_interest, ".pdf")),
      width = 6, height = 5
    )
}

