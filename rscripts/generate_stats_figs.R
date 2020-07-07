# load getopt to allow use of get_Rscript_filename for provenance-gathering
# and sourcing read_count_functions.R
# load here for the same reason
suppressMessages(library(getopt, quietly=T))
suppressMessages(library(here))
# NOTE: other libraries loaded from read_count_functions.R

# FLIC: adding testthat package for +/- temporary unit testing
suppressMessages(library(testthat))


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
  source(here::here("rscripts", "analysis_block_functions.R"))
} else {
  # Deduce file name and path using reflection as before.
  this_script <- getopt::get_Rscript_filename()
  path_to_this_script <- this_script
  source(here::here("rscripts", "provenance.R"))
  source(here::here("rscripts", "read_count_functions.R"))
  source(here::here("rscripts", "analysis_block_functions.R"))
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
coding_seqs <- readDNAStringSet(orf_fasta_file)

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

  # PlotThreeNucleotidePeriodicity()
  three_nucleotide_periodicity_plot <- PlotThreeNucleotidePeriodicity(three_nucleotide_periodicity_data)

  # NOTE: repeated from inside CalculateThreeNucleotidePeriodicity() as preferred not to return multiple objects in list (hassle :S)
  gene_poslen_counts_5start_df <- AllGenes5StartPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hd_file = hd_file, gff_df = gff_df)

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

  # run WriteThreeNucleotidePeriodicity():
  WriteThreeNucleotidePeriodicity(three_nucleotide_periodicity_data)

  print("Completed: Check for 3nt periodicity globally")

} # end ThreeNucleotidePeriodicity() function definition
# run ThreeNucleotidePeriodicity():
ThreeNucleotidePeriodicity(gene_names, dataset, hd_file, gff_df)

#26Mar:
# > ThreeNucleotidePeriodicity(gene_names, dataset, hd_file, gff_df)
# [1] "Starting: Check for 3nt periodicity globally"
# Saving 7 x 7 in image
# [1] "Completed: Check for 3nt periodicity globally"
# Warning message:
#   In write.table(three_nucleotide_periodicity_data, file = tsv_file_path,  :
#                    appending column names to file

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

  # run PlotReadLengths():
  read_len_plot <- PlotReadLengths(read_length_data)
  # creates plot object

  # to run SavePlotReadLenths():
  SavePlotReadLengths(read_len_plot)

  # to run WriteReadLengths():
  WriteReadLengths(read_length_data)

  print("Completed: Distribution of lengths of all mapped reads")

} # end of definition of function DistributionOfLengthsMappedReads()
# run DistributionOfLengthsMappedReads():
DistributionOfLengthsMappedReads(gene_names, dataset, hd_file)

#26Mar:
# > DistributionOfLengthsMappedReads(gene_names, dataset, hd_file)
# [1] "Starting: Distribution of lengths of all mapped reads"
# Saving 7 x 7 in image
# [1] "Completed: Distribution of lengths of all mapped reads"
# Warning message:
#   In write.table(read_length_data, file = tsv_file_path, append = T,  :
#                    appending column names to file


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
  # > str(all_out)
  # 'data.frame':	3690 obs. of  7 variables:
  #   $ Length  : int  10 10 10 10 10 10 10 10 10 10 ...
  # $ Position: int  1 2 3 4 5 6 7 8 9 10 ...
  # $ Frame   : int  0 0 0 0 0 0 0 0 0 0 ...
  # $ A       : num  0 0 0 0 0 0 0 0 0 0 ...
  # $ C       : num  0 0 0 0 0 0 0 0 0 0 ...
  # $ G       : num  0 0 0 0 0 0 0 0 0 0 ...
  # $ T       : num  0 0 0 0 0 0 0 0 0 0 ...

  ## WRITE DATA Biases In Nucleotide Composition Along Mapped Read Lengths
  WriteBiasesInNucleotideComposition(all_out)
  # >   WriteBiasesInNucleotideComposition(all_out)
  # Warning message:
  #   In write.table(all_out, file = tsv_file_path, append = T, sep = "\t",  :
  #                    appending column names to file

  # print("Completed nucleotide composition bias table")

  print("Completed: Biases in nucleotide composition along mapped read lengths")

} # end definition of function: BiasesInNucleotideCompositionAlongMappedReadLengths()
# to run:
#BiasesInNucleotideCompositionAlongMappedReadLengths(gene_names, dataset, hd_file, read_range, min_read_length)

if (do_pos_sp_nt_freq) {
  BiasesInNucleotideCompositionAlongMappedReadLengths(gene_names, dataset, hd_file, read_range, min_read_length)
}
# > BiasesInNucleotideCompositionAlongMappedReadLengths(gene_names, dataset, hd_file, read_range, min_read_length)
# [1] "Starting: Biases in nucleotide composition along mapped read lengths"
# [1] "finished for loop"
# [1] "finished prepping L, R, F variables"
# [1] "finished, returning all_out"
# Warning message:
#   In write.table(all_out, file = tsv_file_path, append = T, sep = "\t",  :
#                    appending column names to file


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

  # run PlotGeneReadFrames():
  gene_read_frame_plot <- PlotGeneReadFrames(gene_read_frames_data)
  # creates plot object

  # run SaveGeneReadFrames():
  SaveGeneReadFrames(gene_read_frame_plot)

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
  
  
  # create empty matrix to store position-specific read counts
  out5p <- matrix(NA, nrow = length(gene_names), ncol = 500) # 5'
  out3p <- matrix(NA, nrow = length(gene_names), ncol = 500) # 3'

  out <- lapply(gene_names, function(gene) {
    GetCodonPositionReads(
      hd_file,
      gene,
      dataset,
      left = (buffer - 15),
      right = (buffer + 11),
      min_read_length = min_read_length
    )
  }) # Get codon-based position-specific reads for each gene
  names(out) <- gene_names

  cc <- 1
  for (gene in gene_names) {
    tmp <- out[[gene]]
    # Only consider genes with at least count_threshold mapped reads along its CDS
    if (sum(tmp) >= count_threshold) {
      tmp <- tmp / mean(tmp)
      if (length(tmp) > 500) {
        out5p[cc, ] <- tmp[1:500]
        out3p[cc, ] <- rev(tmp)[1:500]
      } else {
        out5p[cc, 1:length(tmp)] <- tmp
        out3p[cc, 1:length(tmp)] <- rev(tmp)
      }
    }
    cc <- cc + 1
  }

  # Estimate position-specific mean and std error of mapped read counts
  m5p <- signif(apply(out5p, 2, mean, na.rm = T), 4)
  m3p <- signif(apply(out3p, 2, mean, na.rm = T), 4)
  s5p <- signif(apply(out5p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)
  s3p <- signif(apply(out3p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)

  # Normalize reads to last 50 codons of the 500-codon window.
  # This allows easy comparison between datasets
  s5p <- s5p / mean(m5p[450:500])
  s3p <- s3p / mean(m3p[450:500])
  m5p <- m5p / mean(m5p[450:500])
  m3p <- m3p / mean(m3p[450:500])

  # Create a dataframe to store the output for plots/analyses
  pos_sp_rpf_norm_reads <- data.frame(
    Position = c(1:500, 0:-499),
    Mean = c(m5p, m3p),
    SD = c(s5p, s3p),
    End = factor(rep(c("5'", "3'"), each = 500), levels = c("5'", "3'"))
  )
  
  # Plot
  pos_sp_rpf_norm_reads_plot <- ggplot(
    pos_sp_rpf_norm_reads,
    aes(Position, Mean, col = End)
  ) +
    geom_line() +
    facet_grid(~End, scales = "free") +
    guides(col = FALSE)

  # Save plot and file
  ggsave(pos_sp_rpf_norm_reads_plot, filename = file.path(output_dir, paste0(output_prefix, "pos_sp_rpf_norm_reads.pdf")))
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "pos_sp_rpf_norm_reads.tsv"))
  write_provenance_header(this_script, tsv_file_path)
  write.table(
    pos_sp_rpf_norm_reads,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
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
  pos_sp_mrna_norm_coverage <- CalculateNucleotideBasedPositionSpecificReads(gene, dataset, min_read_length, read_range, buffer)
  
  # plot
  PlotNucleotideBasedPositionSpecificReadsPerGene(pos_sp_mrna_norm_coverage)
  
  # save plot out
  SaveNucleotideBasedPositionSpecificReadsPerGene(pos_sp_mrna_norm_coverage_plot)
  
  # write file out
  WriteNucleotideBasedPositionSpecificReadsPerGene(pos_sp_mrna_norm_coverage)

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
GeneTranscriptsPerMillion <- function(gene, dataset, hd_file){

  print("Starting: Calculate TPMs of genes")

    tpms <- CalculateGeneTranscriptsPerMillion(gene, dataset, hd_file)

    WriteGeneTranscriptsPerMillion(tpms)

    print("Completed: Calculate TPMs of genes")

} # end GeneTranscriptsPerMillion() definition
# run GeneTranscriptsPerMillion():
GeneTranscriptsPerMillion(gene, dataset, hd_file)


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

  tpms <- CalculateGeneTranscriptsPerMillion(gene, dataset, hd_file) # repeated from TPMs section to make it available here

  features_plot_data <- CalculateSequenceBasedFeatures(features, tpms)

  features_plot <- PlotSequenceBasedFeatures(features_plot_data)

  WriteSequenceBasedFeatures(features_plot)

  print("Completed: Correlations between TPMs of genes with their sequence-based features")

} else { # skip

  print("Skipped: Correlations between TPMs of genes with their sequence-based features - features_file.tsv not provided")

}

## Codon-specific ribosome densities for correlations with tRNAs


# BIG FUNCTIONS:


# Codon-specific ribosome density for tRNA correlation; skip if missing t_rna_file & codon_positions_file
if (!is.na(t_rna_file) & !is.na(codon_positions_file)) {

  print("Starting: Codon-specific ribosome densities for correlations with tRNAs")
  # Only for RPF datasets

  if (rpf) {

    # TODO: This section needs attention. Can be refactored analogously to reads_per_codon_etc
    # Needs separate calculation of per-codon normalized reads
    # WAITING: we want new format of codon_pos from @john-s-f before editing this chunk

  #   cod_dens_tRNA_data <- CalculateCodonSpecificRibosomeDensity(t_rna_file, hd_file, gene, dataset, buffer, count_threshold)
  # 
  #   cod_dens_tRNA_plot <- PlotCodonSpecificRibosomeDensityTRNACorrelation(cod_dens_tRNA_data)
  # 
  #   SaveCodonSpecificRibosomeDensityTRNACorrelation(cod_dens_tRNA_plot)
  # 
  #   WriteCodonSpecificRibosomeDensityTRNACorrelation(cod_dens_tRNA_data)
  # 
  # }


    # This still depends on yeast-specific arguments and should be edited.
    yeast_tRNAs <- read.table(t_rna_file, h = T) # Read in yeast tRNA estimates
    load(codon_positions_file) # Position of codons in each gene (numbering ignores first 200 codons)
    # Reads in an object named "codon_pos"
    out <- lapply(gene_names, function(gene) {
      # From "Position specific distribution of reads" plot
      GetCodonPositionReads(hd_file, gene, dataset, left = (buffer - 15), right = (buffer + 11), min_read_length = min_read_length)
    }) # Get codon-based position-specific reads for each gene
    names(out) <- gene_names

    gene_len <- sapply(out, length) # Calculate gene length in codons
    out <- out[gene_len > 201] # Ignore genes with <=200 sense codons

    trim_out <- lapply(out, function(x) {
      x[201:(length(x) - 1)]
    }) # Trim first 200 codons and stop codon from each gene
    read_counts_trim <- sapply(trim_out, sum) # Calculate read counts in trimmed genes
    trim_out <- trim_out[read_counts_trim >= count_threshold] # Ignore genes with fewer than count_threshold mapped reads

    norm_out <- lapply(trim_out, function(x) {
      x / mean(x)
    }) # Normalize reads in each gene by their mean

    # Calculate codon-specific mean ribosome-densities at A/P/E sites of the mapped reads
    a_mn <- sapply(names(codon_pos), function(codon) {
      mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
        pos <- as.numeric(a[2])
        norm_out[[a[1]]][pos]
      })), na.rm = T)
    })
    p_mn <- sapply(names(codon_pos), function(codon) {
      mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
        pos <- as.numeric(a[2]) + 1
        norm_out[[a[1]]][pos]
      })), na.rm = T)
    })
    e_mn <- sapply(names(codon_pos), function(codon) {
      mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
        pos <- as.numeric(a[2]) + 2
        norm_out[[a[1]]][pos]
      })), na.rm = T)
    })

    # Sort the values
    A <- a_mn[order(names(codon_pos))]
    P <- p_mn[order(names(codon_pos))]
    E <- e_mn[order(names(codon_pos))]

    cod_dens_tRNA <- cbind(yeast_tRNAs, A, P, E)

    # Prepare data for plot
    cod_dens_tRNA_wide <- cod_dens_tRNA %>%
      gather(tRNA_type, tRNA_value, 3:6) %>%
      gather(Site, Ribodens, 3:5)

    cod_dens_tRNA_plot <- ggplot(cod_dens_tRNA_wide, aes(x = tRNA_value, y = Ribodens)) +
      geom_point(alpha = 0.3) +
      facet_grid(Site ~ tRNA_type, scales = "free_x") +
      geom_smooth(method = "lm") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Save plot and file
    ggsave(cod_dens_tRNA_plot, filename = file.path(output_dir, paste0(output_prefix, "codon_ribodens.pdf")))
    tsv_file_path <- file.path(output_dir, paste0(output_prefix, "codon_ribodens.tsv"))
    write_provenance_header(this_script, tsv_file_path)
    write.table(
      cod_dens_tRNA,
      file = tsv_file_path,
      append = T,
      sep = "\t",
      row = F,
      col = T,
      quote = F
    )
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
