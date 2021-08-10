#' Wrapper to invoke Rmarkdown script 'rmarkdown/AnalysisOutputs.Rmd'.

suppressMessages(library(getopt, quietly = T))
suppressMessages(library(here, quietly = T))
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
  make_option("--rmd",
              type = "character", default = NA,
              help = "AnalysisOutputs.Rmd file"
  ),
  make_option("--yamlfile",
              type = "character", default = NA,
              help = "YAML configuration file"
  ),
  make_option("--sampleid",
              type = "character", default = NA,
              help = "Sample ID"
  ),
  make_option("--three-nucleotide-periodicity-data-file",
              type = "character", default = NA,
              help = "3nt periodicity TSV file"
  ),
  make_option("--gene-position-length-counts-5start-file",
              type = "character", default = NA,
              help = "Gene position length counts 5 start TSV file"
  ),
  make_option("--read-length-data-file",
              type = "character", default = NA,
              help = "Read lengths TSV file"
  ),
  make_option("--pos-sp-rpf-norm-reads-data-file",
              type = "character", default = NA,
              help = "Pos sp RPF norm reads TSV file"
  ),
  make_option("--gene-read-frames-filtered-data-file",
              type = "character", default = NA,
              help = "Gene reads filtered TSV file"
  ),
  make_option("--codon-ribodens-gathered-file",
              type = "character", default = NA,
              help = "Codons ribosome densities gathered TSV file"
  ),
  make_option("--sequence-features-file",
              type = "character", default = NA,
              help = "Sequence features TSV file"
  ),
  make_option("--intermediates-dir",
              type = "character", default = NA,
              help = "Intermediate files directory"
  ),
  make_option("--output-format",
              type = "character", default = "html_document",
              help = "Output format"
  ),
  make_option("--output-file",
              type = "character", default = NA,
              help = "Output HTML file"
  )
)

print_provenance(get_Rscript_filename())

opt <- parse_args(OptionParser(option_list = option_list),
                  convert_hyphens_to_underscores = TRUE)
print("analysis_outputs.R running with parameters:")
print(opt)

rmarkdown::render(
  opt$rmd,
  params = list(
    verbose = FALSE,
    yamlfile = opt$yamlfile,
    sampleid = opt$sampleid,
    three_nucleotide_periodicity_data_file = opt$three_nucleotide_periodicity_data_file,
    gene_position_length_counts_5start_file = opt$gene_position_length_counts_5start_file,
    read_length_data_file = opt$read_length_data_file,
    pos_sp_rpf_norm_reads_data_file = opt$pos_sp_rpf_norm_reads_data_file,
    gene_read_frames_filtered_data_file = opt$gene_read_frames_filtered_data_file,
    codon_ribodens_gathered_file = opt$codon_ribodens_gathered_file,
    sequence_features_file = opt$sequence_features_file),
  intermediates_dir = opt$intermediates_dir,
  output_format = opt$output_format,
  output_file = opt$output_file)
