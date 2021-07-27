# load getopt to allow use of get_Rscript_filename for provenance-gathering
# and sourcing read_count_functions.R
# load here for the same reason
suppressMessages(library(getopt, quietly=T))
suppressMessages(library(here))
# NOTE: other libraries loaded from read_count_functions.R

# Handle interactive session behaviours or use get_Rscript_filename():
if (interactive()) {
  # Use hard-coded script name and assume script is in "rscripts"
  # directory. This assumes that interactive R is being run within
  # the parent of rscripts/ but imposes no other constraints on
  # where rscripts/ or its parents are located.
  this_script <- "generate_stats_figs.R"
  path_to_this_script <- here("rscripts", this_script)
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

# print provenance
print_provenance(path_to_this_script)

# define input options for optparse package
option_list <- list(
  make_option("--output-dir",
              type = "character", default = "./",
              help = "Output directory"
  ),
  make_option("--orf-fasta-file",
              type = "character", default = FALSE,
              help = "FASTA file with nt seq"
  ),
  make_option("--orf-gff-file",
              type = "character", default = NA,
              help = "riboviz generated GFF2/GFF3 annotation file"
  ),
  make_option("--output-pdfs",
              type = "logical", default = TRUE,
              help = "generate .pdfs for sample-related plots"
  ),
  make_option("--num-processes",
              type = "integer", default = 1,
              help = "Number of cores for parallelization"
  ),
  make_option("--min-read-length",
              type = "integer", default = 10,
              help = "Minimum read length in H5 output"
  ),
  make_option("--max-read-length",
              type = "integer", default = 50,
              help = "Maximum read length in H5 output"
  ),
  make_option("--buffer",
              type = "integer", default = 250,
              help = "Length of flanking region around the CDS"
  ),
  make_option("--primary-id",
              type = "character", default = "gene_id",
              help = "Primary gene IDs to access the data (YAL001C, YAL003W, etc.)"
  ),
  make_option("--dataset",
              type = "character", default = "vignette",
              help = "Name of the dataset"
  ),
  make_option("--rpf",
              type = "logical", default = TRUE,
              help = "Is the dataset an RPF or mRNA dataset?"
  ),
  make_option("--features-file",
              type = "character", default = NA,
              help = "features file, columns are gene features and rows are genes"
  ),
  make_option("--do-pos-sp-nt-freq",
              type = "logical", default = TRUE,
              help = "do calculate the position-specific nucleotide frequency"
  ),
  make_option("--t-rna-file",
              type = "character", default = NA,
              help = "tRNA estimates in .tsv file"
  ),
  make_option("--codon-positions-file",
              type = "character", default = NA,
              help = "Codon positions in each gene in .Rdata file"
  ),
  make_option("--count-threshold",
              type = "integer", default = 64,
              help = "threshold for count of reads per gene to be included in plot"
  ),
  make_option("--length_threshold",
              type = "integer", default = 500,
              help = "threshold for length of genes to be included in metagene plots"
  ),
  make_option("--output-prefix",
              type = "character", default = "",
              help = "Prefix for output files"
  ),
  make_option("--hd-file",
              type = "character", default = "output.h5",
              help = "Location of H5 output file"
  ),
  make_option("--nnt-buffer",
              type = "integer", default = 25,
              help = "n nucleotides of UTR buffer to include in metagene plots"
  ),
  make_option("--nnt-gene",
              type = "integer", default = 50,
              help = "n nucleotides of gene to include in metagene plots"
  ),
  make_option("--asite-disp-length-file",
              type = "character", default = NA,
              help = "asite displacement file
    table with one displacement per read length"
  )
)

# read in commandline arguments
opt <- optparse::parse_args(OptionParser(option_list = option_list),
                            convert_hyphens_to_underscores=TRUE)

# attach opt list to be able to refer to variables in the list by names alone
# ie `height` rather than `women$height`
attach(opt)

print("generate_stats_figs.R running with parameters:")
opt

# read in positions of all exons/genes in GFF format and subset CDS locations
gene_names <- rhdf5::h5ls(hd_file, recursive = 1)$name

# read in coding sequences
coding_seqs <- Biostrings::readDNAStringSet(orf_fasta_file)

# range of read lengths between parameters set in config file
read_range <- min_read_length:max_read_length

# read in positions of all exons/genes in GFF format and convert to tibble data frame
gff_df <- readGFFAsDf(orf_gff_file)

# set ggplot2 theme for plots drawn after this; use dark on light theme
ggplot2::theme_set(theme_bw())

####

# REFACTORING NOTE: Breaking functions into large & medium chunks:

# LARGE CHUNKS:
#  3nt periodicity
#  Lengths all mapped reads
#  Biases in nucleotide composition
#  Calculate read frame for every orf
#  RPF position specific distribution of reads
#  mRNA position specific distribution of reads
#  TPMs of genes
#  TPMs correlations with features

# MEDIUM CHUNKS:
# Following this naming convention:
#  CalculateFunction - computational code
#  PlotFunction - plot creation
#  SaveFunction - PDF writeouts of plots
#  WriteFunction - .TSV file creation

# NOTE: Do not use mclapply when accessing H5 data

#
#
# START 3NT PERIODICITY
#
# 

ThreeNucleotidePeriodicity <- function(gene_names, dataset, hd_file, gff_df) {

  # check for 3nt periodicity
  print("Starting: Check for 3nt periodicity globally")

  # CalculateThreeNucleotidePeriodicity():
  three_nucleotide_periodicity_data <- CalculateThreeNucleotidePeriodicity(gene_names = gene_names, dataset = dataset, hd_file = hd_file, gff_df = gff_df)

  # CalculateGenePositionLengthCounts5Start():
  gene_poslen_counts_5start_df <- CalculateGenePositionLengthCounts5Start(gene_names = gene_names, dataset = dataset, hd_file = hd_file, gff_df = gff_df)
  
  # PlotThreeNucleotidePeriodicity()
  if(output_pdfs){
    three_nucleotide_periodicity_plot <- PlotThreeNucleotidePeriodicity(three_nucleotide_periodicity_data)
    
    # run PlotStartCodonRiboGrid()
    start_codon_ribogrid_plot <- PlotStartCodonRiboGrid(gene_poslen_counts_5start_df)
    # creates plot object
    
    # run SaveStartCodonRiboGrid():
    SaveStartCodonRiboGrid(start_codon_ribogrid_plot)
    
    # run PlotStartCodonRiboGridBar():
    start_codon_ribogrid_bar_plot <- PlotStartCodonRiboGridBar(gene_poslen_counts_5start_df)
    # creates plot object
    
    # run SaveStartCodonRiboGridBar():
    SaveStartCodonRiboGridBar(start_codon_ribogrid_bar_plot)
    
    # run SavePlotThreeNucleotidePeriodicity():
    SavePlotThreeNucleotidePeriodicity(three_nucleotide_periodicity_plot)
  } 

  # run WriteThreeNucleotidePeriodicity():
  WriteThreeNucleotidePeriodicity(three_nucleotide_periodicity_data)
  
  # WriteGenePositionLengthCounts5Start(): 
  WriteGenePositionLengthCounts5Start(gene_poslen_counts_5start_df)

  print("Completed: Check for 3nt periodicity globally")

} # end ThreeNucleotidePeriodicity() function definition

# run ThreeNucleotidePeriodicity():
ThreeNucleotidePeriodicity(gene_names, dataset, hd_file, gff_df)

#
#
# END 3NT PERIODICITY
#
# START ALL MAPPED READS
#
#

DistributionOfLengthsMappedReads <- function(gene_names, dataset, hd_file){

  # run CalculateReadLengths():
  read_length_data <- CalculateReadLengths(gene_names, dataset, hd_file)

  if(output_pdfs){
    # run PlotReadLengths():
    read_len_plot <- PlotReadLengths(read_length_data)
    # creates plot object
    
    # to run SavePlotReadLenths():
    SavePlotReadLengths(read_len_plot)
  }
  
  # to run WriteReadLengths():
  WriteReadLengths(read_length_data)

  print("Completed: Distribution of lengths of all mapped reads")

} # end of definition of function DistributionOfLengthsMappedReads()

# run DistributionOfLengthsMappedReads():
DistributionOfLengthsMappedReads(gene_names, dataset, hd_file)

#
#
# END ALL MAPPED READS
#
# START BIASES IN NUCLEOTIDE COMPOSITION
#
#

if (!do_pos_sp_nt_freq) {

  print("NOT calculating position-specific nucleotide frequency - reason: do_pos_sp_nt_freq parameter set to FALSE")

} 

# big function with arguments
BiasesInNucleotideCompositionAlongMappedReadLengths <- function(gene_names, dataset, hd_file, read_range, min_read_length) {

  print("Starting: Biases in nucleotide composition along mapped read lengths")

  ## CALCULATE Biases In Nucleotide Composition Along Mapped Read Lengths
  all_out <- CalculateBiasesInNucleotideComposition(gene_names, dataset, hd_file, read_range, min_read_length)
  
  ## WRITE DATA Biases In Nucleotide Composition Along Mapped Read Lengths
  WriteBiasesInNucleotideComposition(all_out)
  
  print("Completed: Biases in nucleotide composition along mapped read lengths")

} # end definition of function: BiasesInNucleotideCompositionAlongMappedReadLengths()

if (do_pos_sp_nt_freq) {
  
  BiasesInNucleotideCompositionAlongMappedReadLengths(gene_names, dataset, hd_file, read_range, min_read_length)

  }

#
#
# END BIASES IN NUCLEOTIDE COMPOSITION
#
# START CALCULATE READ FRAME FOR EVERY ORF
#
#

## calculate read frame for every annotated ORF

if (is.na(asite_disp_length_file)) {
  
    print("NOT checking for 3nt periodicity (frame) by gene - reason: asite_disp_length_file parameter not provided")
  
} # TODO: FLIC FIGURE OUT IF USING ELSE OR NOT # else { # if asite_disp_length_file parameter provided, calculate read frame for every ORF:

if (!is.na(asite_disp_length_file)) {

  # check frame by gene
  print("Starting: Check for 3nt periodicity (frame) by gene")

  # get a-site displacement lengths from file
  asite_displacement_length <- ReadAsiteDisplacementLengthFromFile(asite_disp_length_file)

  # run CalculateGeneReadFrames() to create data object
  gene_read_frames_data <- CalculateGeneReadFrames(dataset, hd_file, gff_df, min_read_length, asite_displacement_length)

  # filter gene_read_frames_data to remove counts over the count_threshold
  gene_read_frame_data_filtered <- FilterGeneReadFrames(gene_read_frames_data, count_threshold)
  
  if(output_pdfs){
    # run PlotGeneReadFrames():
    gene_read_frame_plot <- PlotGeneReadFrames(gene_read_frame_data_filtered)
    # creates plot object
    
    # run SaveGeneReadFrames():
    SaveGeneReadFrames(gene_read_frame_plot)
  }

  # run WriteFilteredGeneReadFrames():
  WriteFilteredGeneReadFrames(gene_read_frame_data_filtered)
  
  # run WriteGeneReadFrames():
  WriteGeneReadFrames(gene_read_frames_data)

  print("Completed: Check for 3nt periodicity (frame) by Gene")
}

#
#
# END CALCULATE READ FRAME FOR EVERY ORF
#
# START RPF POSITION SPECIFIC DISTRIBUTION OF READS
#
#

print("Starting: Position specific distribution of reads")

# For RPF datasets, generate codon-based position-specific reads
if (rpf) {
  
  pos_sp_rpf_norm_reads_data <- CalculatePositionSpecificDistributionOfReads(gene_names, dataset, hd_file, gff_df, min_read_length, count_threshold, asite_displacement_length)
  
  if(output_pdfs){
    pos_sp_rpf_norm_reads_plot <- PlotPositionSpecificDistributionOfReads(pos_sp_rpf_norm_reads_data)
    
    SavePositionSpecificDistributionOfReads(pos_sp_rpf_norm_reads_plot)
  }

  WritePositionSpecificDistributionOfReads(pos_sp_rpf_norm_reads_data)

}

print("Completed: Position specific distribution of reads - RPF method")

#
#
# END RPF POSITION SPECIFIC DISTRIBUTION OF READS
#
# START MRNA POSITION SPECIFIC DISTRIBUTION OF READS
#
#

# run mRNA dataset method for position specific distribution of reads
 # (nucleotide-based instead of codon-based as per RPF method)

if (!rpf) {

  print("Starting: Position specific distribution of reads - mRNA dataset method")

  # TODO: rewrite along the lines of rpf above
  warning("Warning: nt-based position-specific reads for non-RPF datasets not tested")

  # calculate
  pos_sp_mrna_norm_coverage <- CalculateNucleotideBasedPositionSpecificReadsMRNA(gene_names, dataset, min_read_length, read_range, buffer)
  
  if(output_pdfs){
    # plot
    PlotNucleotideBasedPositionSpecificReadsPerGeneMRNA(pos_sp_mrna_norm_coverage)
    
    # save plot out
    SaveNucleotideBasedPositionSpecificReadsPerGeneMRNA(pos_sp_mrna_norm_coverage_plot)
  }
  
  # write file out
  WriteNucleotideBasedPositionSpecificReadsPerGeneMRNA(pos_sp_mrna_norm_coverage)

  print("Completed: Position specific distribution of reads - mRNA dataset method")

}

print("Completed: Position specific distribution of reads")


#
#
# END MRNA POSITION SPECIFIC DISTRIBUTION OF READS
#
# START TPMS OF GENES
#
#

# Calculate TPMs of genes
GeneTranscriptsPerMillion <- function(gene_names, dataset, hd_file){
  
  print("Starting: Calculate TPMs of genes")
  
  tpms <- CalculateGeneTranscriptsPerMillion(gene_names, dataset, hd_file)
  
  WriteGeneTranscriptsPerMillion(tpms)
  
  print("Completed: Calculate TPMs of genes")
  
} # end GeneTranscriptsPerMillion() definition

# run GeneTranscriptsPerMillion():
GeneTranscriptsPerMillion(gene_names, dataset, hd_file)


#
#
# END TPMS OF GENES
#
# START TPMS CORRELATIONS WITH FEATURES
#
#

## Correlations between TPMs of genes with their sequence-based features

# Correlate TPMs of genes with sequence-based features, skip if missing features_file
if (!is.na(features_file)) { # do correlating

  print("Starting: Correlations between TPMs of genes with their sequence-based features")

  features <- ReadSequenceBasedFeatures(features_file)

  tpms <- CalculateGeneTranscriptsPerMillion(gene_names, dataset, hd_file) # repeated from TPMs section to make it available here

  features_plot_data <- CalculateSequenceBasedFeatures(features, tpms)

  if (output_pdfs){
    features_plot <- PlotSequenceBasedFeatures(features_plot_data)
    
    SaveSequenceBasedFeatures(features_plot)
  }

  WriteSequenceBasedFeatures(features_plot_data)

  print("Completed: Correlations between TPMs of genes with their sequence-based features")

} else { # skip

  print("Skipped: Correlations between TPMs of genes with their sequence-based features - features_file.tsv not provided")

}

## Codon-specific ribosome densities for correlations with tRNAs

# Codon-specific ribosome density for tRNA correlation; skip if missing t_rna_file & codon_positions_file
if (!is.na(t_rna_file) & !is.na(codon_positions_file)) {

  print("Starting: Codon-specific ribosome densities for correlations with tRNAs")
  # Only for RPF datasets

  if (rpf) {
    
    cod_dens_tRNA_data <- CalculateCodonSpecificRibosomeDensity(t_rna_file, codon_positions_file, gene_names, hd_file, dataset, gff_df, count_threshold, asite_displacement_length)
    
    cod_dens_tRNA_wide <- GatherCodonSpecificRibosomeDensityTRNACorrelation(cod_dens_tRNA_data)
    
    if(output_pdfs){
      cod_dens_tRNA_plot <- PlotCodonSpecificRibosomeDensityTRNACorrelation(cod_dens_tRNA_wide)
      
      SaveCodonSpecificRibosomeDensityTRNACorrelation(cod_dens_tRNA_plot)
    }
    
    WriteGatheredCodonSpecificRibosomeDensityTRNACorrelation(cod_dens_tRNA_wide)
    
    WriteCodonSpecificRibosomeDensityTRNACorrelation(cod_dens_tRNA_data)
    
  }

  print("Completed: Codon-specific ribosome densities for correlations with tRNAs")

} else {

  print("Skipped: Codon-specific ribosome densities for correlations with tRNAs - t-rna-file and/or codon-positions-file not provided")
}

#
#
# END TPMS CORRELATIONS WITH FEATURES
#
#

print("Completed generate_stats_figs.R")
