#' Convert BAM files to RiboViz HDF5 files.
#'
#' Given a GFF file and a BAM file, this script creates an HDF5 file
#' with information about a feature (e.g. CDS, ORF, or uORF).
#;
#' This script accepts the following command-line parameters:
#'
#' | Parameter | Description |
#' | --------- | ----------- |
#' | `bam-file` | BAM input file |
#' | `orf-gff-file` | Matched genome feature file, specifying coding sequences locations (start and stop coordinates) within the transcripts (GTF/GFF3 file) |
#' | `feature` | Feature e.g. `CDS`, `ORF`, or `uORF` |
#' | `min-read-length` | Minimum read length in H5 output |
#' | `max-read-length` | Maximum read length in H5 output |
#' | `buffer` | Length of flanking region around the feature. Used only if `is_riboviz_gff` is `FALSE`. |
#' | `primary-id` | Primary gene IDs to access the data (YAL001C, YAL003W, etc.) |
#' | `secondary-id` | Secondary gene IDs to access the data (COX1, EFB1, etc. or `NULL`) |
#' | `dataset` | Human-readable name of the dataset |
#' | `stop-in-cds` | Are stop codons part of the feature annotations in `orf-gff-file`? Used only if `is_riboviz_gff` is `FALSE`. |
#' | `is-riboviz-gff` | Does `orf-gff-file` contain 3 elements per gene - UTR5, feature, and UTR3? |
#' | `num-processes` | Number of processes to parallelize over |
#' | `hd-file` | H5 output file |
#'
#' See `BamToH5` in `bam_to_h5_functions.R` for details of behaviour.
#'
#' @export

suppressMessages(library(getopt, quietly = T))
suppressMessages(library(here))
suppressMessages(library(optparse, quietly = T))

# Load local dependencies.
if (interactive()) {
  # Use hard-coded script name and assume script is in "rscripts"
  # directory. This assumes that interactive R is being run within
  # the parent of rscripts/ but imposes no other constraints on
  # where rscripts/ or its parents are located.
  self <- "bam_to_h5.R"
  path_to_self <- here("rscripts", self)
  source(here::here("rscripts", "provenance.R"))
  source(here::here("rscripts", "bam_to_h5_functions.R"))
} else {
  # Deduce file name and path using reflection as before.
  self <- getopt::get_Rscript_filename()
  path_to_self <- self
  source(file.path(dirname(self), "provenance.R"))
  source(file.path(dirname(self), "bam_to_h5_functions.R"))
}

option_list <- list(
  make_option("--bam-file", type = "character", default = NA,
    help = "BAM input file"),
  make_option("--orf-gff-file", type = "character", default = NA,
    help = "GFF2/GFF3 Matched genome feature file, specifying coding sequences locations (start and stop coordinates) within the transcripts"),
  make_option("--feature", type = "character", default = "CDS",
    help = "Feature e.g. CDS, ORF, or uORF"),
  make_option("--min-read-length", type = "integer", default = 10,
    help = "Minimum read length in H5 output"),
  make_option("--max-read-length", type = "integer", default = 50,
    help = "Maximum read length in H5 output"),
  make_option("--buffer", type = "integer", default = 250,
    help = "Length of flanking region around the feature"),
  make_option("--primary-id", type = "character", default = "gene_id",
    help = "Primary gene IDs to access the data (YAL001C, YAL003W, etc.)"),
  make_option("--secondary-id", type = "character", default = NA,
    help = "Secondary gene IDs to access the data (COX1, EFB1, etc.)"),
  make_option("--dataset", type = "character", default = "data",
    help = "Human-readable name of the dataset"),
  make_option("--stop-in-cds", type = "logical", default = FALSE,
    help = "Are stop codons part of the feature annotations in GFF?"),
  make_option("--is-riboviz-gff", type = "logical", default = TRUE,
    help = "Does the GFF file contain 3 elements per gene - UTR5, feature, and UTR3?"),
  make_option("--hd-file", type = "character", default = "output.h5",
    help = "H5 output file"),
  make_option("--num-processes", type = "integer", default = 1,
    help = "Number of processes to parallelize over")
)

print_provenance(get_Rscript_filename())

opt <- parse_args(OptionParser(option_list = option_list),
                  convert_hyphens_to_underscores = TRUE)
print("bam_to_h5.R running with parameters:")
print(opt)

if (is.na(opt$bam_file)) {
  stop("--bam-file parameter must be provided. See usage (--help)")
}
if (is.na(opt$orf_gff_file)) {
  stop("--orf_gff_file parameter must be provided. See usage (--help)")
}

BamToH5(opt$bam_file, opt$orf_gff_file, opt$feature,
        opt$min_read_length, opt$max_read_length,
        opt$buffer, opt$primary_id, opt$secondary_id, opt$dataset,
        opt$stop_in_cds, opt$is_riboviz_gff, opt$hd_file,
        opt$num_processes)

print("bam_to_h5.R done")
