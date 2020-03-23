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
  source(here("rscripts", "provenance.R"))
  source(here("rscripts", "read_count_functions.R"))
} else {
  # Deduce file name and path using reflection as before.
  this_script <- getopt::get_Rscript_filename()
  path_to_this_script <- this_script
  source(file.path(dirname(this_script), "provenance.R")) # FLIC: consider if this can be redone using here() too
  source(file.path(dirname(this_script), "read_count_functions.R")) # FLIC: consider if this can be redone using here() too
}

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

# print provenance
print_provenance(path_to_this_script)

# read in commandline arguments
opt <- optparse::parse_args(OptionParser(option_list = option_list),
                            convert_hyphens_to_underscores=TRUE)

# attach opt list to be able to refer to variables in the list by names alone
 # ie `height` rather than `women$height`
attach(opt)

print("generate_stats_figs.R running with parameters:")
opt

# prepare files, opens hdf5 file connection
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

# list of gene names taken from h5 file
gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name

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
# Following this naming convention: 
 #  CalculateFunction - computational code
 #  PlotFunction - plot creation
 #  SaveFunction - PDF writeouts of plots
 #  WriteFunction - .TSV file creation
s
#####

#
#
# START 3NT PERIODICITY
#
#

# big function with probable arguments
ThreeNucleotidePeriodicity <- function(gene_names, dataset, hdf5file, gff_df) {
  
  # check for 3nt periodicity
  print("Starting: Check for 3nt periodicity globally")
  
  # NOTE: Do not use mclapply when accessing H5 data
  
  CalculateThreeNucleotidePeriodicity <- function(gene_names, dataset, hdf5file, gff_df){  
    
    # get gene and position specific total counts for all read lengths
    gene_poslen_counts_5start_df <- AllGenes5StartPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hdf5file = hdf5file, gff_df = gff_df)
    
    gene_poslen_counts_3end_df <- AllGenes3EndPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hdf5file = hdf5file, gff_df = gff_df)
   
    # summarize by adding different read lengths
    gene_pos_counts_5start <- gene_poslen_counts_5start_df %>%
      group_by(Pos) %>%
      summarize(Counts = sum(Counts))
    # gives: 
    # > str(gene_pos_counts_5start)
    # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	75 obs. of  2 variables:
    #   $ Pos   : int  -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 ...
    #   $ Counts: int  285 318 307 386 291 347 840 330 475 355 ...
    
    gene_pos_counts_3end <- gene_poslen_counts_3end_df  %>%
      group_by(Pos) %>%
      summarize(Counts = sum(Counts))
    # gives: 
    # > str(gene_pos_counts_3end)
    # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	75 obs. of  2 variables:
    #   $ Pos   : int  -49 -48 -47 -46 -45 -44 -43 -42 -41 -40 ...
    #   $ Counts: int  19030 13023 50280 19458 12573 46012 19043 13282 36968 20053 ...
    
    three_nucleotide_periodicity_data <- bind_rows(
      gene_pos_counts_5start %>% mutate(End = "5'"),
      gene_pos_counts_3end %>% mutate(End = "3'")
    ) %>%
      mutate(End = factor(End, levels = c("5'", "3'")))
    # gives: 
    # > str(three_nucleotide_periodicity_data)
    # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	150 obs. of  3 variables:
    #   $ Pos   : int  -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 ...
    #   $ Counts: int  285 318 307 386 291 347 840 330 475 355 ...
    #   $ End   : Factor w/ 2 levels "5'","3'": 1 1 1 1 1 1 1 1 1 1 ...
    
    return(three_nucleotide_periodicity_data)
    
  } # end CalculateThreeNucleotidePeriodicity() definition
  
  # gives: 
  #   CalculateThreeNucleotidePeriodicity(gene_names = gene_names, dataset = dataset, hdf5file = hdf5file, gff_df = gff_df)
  #   # A tibble: 150 x 3
  #   Pos Counts End  
  #   <int>  <int> <fct>
  #     1   -24    285 5'   
  #     2   -23    318 5'   
  #     3   -22    307 5'   
  #     4   -21    386 5'   
  #     5   -20    291 5'   
  #     6   -19    347 5'   
  #     7   -18    840 5'   
  #     8   -17    330 5'   
  #     9   -16    475 5'   
  #    10   -15    355 5'   
  #   # … with 140 more rows
  
  # run the function:
  three_nucleotide_periodicity_data <- CalculateThreeNucleotidePeriodicity(gene_names = gene_names, dataset = dataset, hdf5file = hdf5file, gff_df = gff_df)
  
  # define PlotThreeNucleotidePeriodicity() function with reasonable arguments
  PlotThreeNucleotidePeriodicity <- function(three_nucleotide_periodicity_data){
    
    # Plot
    three_nucleotide_periodicity_plot <- ggplot(
      three_nucleotide_periodicity_data,
      aes(x = Pos, y = Counts)) +
      geom_line() +
      facet_wrap(~End, scales = "free") +
      labs(x = "Nucleotide Position", y = "Read counts")
    
    return(three_nucleotide_periodicity_plot) 
    
  } # end PlotThreeNucleotidePeriodicity() definition
  
  # run:
  #PlotThreeNucleotidePeriodicity(three_nucleotide_periodicity_data)
  # creates plot object
  
  three_nucleotide_periodicity_plot <- PlotThreeNucleotidePeriodicity(three_nucleotide_periodicity_data)
  
  # NOTE: repeated from inside CalculateThreeNucleotidePeriodicity() as preferred not to return multiple objects in list (hassle :S)
  gene_poslen_counts_5start_df <- AllGenes5StartPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hdf5file = hdf5file, gff_df = gff_df)
  
  
  # potentially replace/tweak plot_ribogrid() to follow StyleGuide
  PlotStartCodonRiboGrid <- function(gene_poslen_counts_5start_df){
  # function to do the ribogrid & ribogridbar plots?
  # ribogrid_5start
    start_codon_ribogrid_plot <- plot_ribogrid(gene_poslen_counts_5start_df)
    return(start_codon_ribogrid_plot) 
  } # end PlotStartCodonRiboGrid() definition
  
  # run: 
  start_codon_ribogrid_plot <- PlotStartCodonRiboGrid(gene_poslen_counts_5start_df)
  # creates plot object
  
  
  SaveStartCodonRiboGrid <- function(start_codon_ribogrid_plot){
    # function to do the ribogrid & ribogridbar plots?
    # ribogrid_5start
    start_codon_ribogrid_plot %>%
      ggsave(
        filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogrid.pdf")),
        width = 6, height = 3
      )
    #return() # no return as writing-out
  } # end SaveStartCodonRiboGrid() definition
  
  # run: 
  SaveStartCodonRiboGrid(start_codon_ribogrid_plot)
  
  
  PlotStartCodonRiboGridBar <- function(gene_poslen_counts_5start_df){  
  start_codon_ribogrid_bar_plot <- barplot_ribogrid(gene_poslen_counts_5start_df)
      return(start_codon_ribogrid_bar_plot) 
  } # end PlotStartCodonRiboGridBar() definition
  
  # run: 
  start_codon_ribogrid_bar_plot <- PlotStartCodonRiboGridBar(gene_poslen_counts_5start_df)
  # creates plot object
  
  
  SaveStartCodonRiboGridBar <- function(start_codon_ribogrid_bar_plot){  
    start_codon_ribogrid_bar_plot %>%  
    ggsave(
        filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogridbar.pdf")),
        width = 6, height = 5
      )
    #return() # no return as writing-out
  } # end SaveStartCodonRiboGridBar() definition
  
  # run: 
  SaveStartCodonRiboGrid(start_codon_ribogrid_bar_plot)
  
  
  SavePlotThreeNucleotidePeriodicity <- function() {
    
    # Save plot and file
    ggsave(
      three_nucleotide_periodicity_plot, 
      filename = file.path(output_dir, paste0(output_prefix, "3nt_periodicity.pdf"))
    )
    # return() # NO RETURN as writing out
    
  } # end of function definition SavePlotThreeNucleotidePeriodicity()
  
  
  WriteThreeNucleotidePeriodicity <- function() {
    
    tsv_file_path <- file.path(output_dir, paste0(output_prefix, "3nt_periodicity.tsv"))
    
    write_provenance_header(path_to_this_script, tsv_file_path)
    
    write.table(
      three_nucleotide_periodicity_data,
      file = tsv_file_path,
      append = T,
      sep = "\t",
      row = F,
      col = T,
      quote = F)
    
    # return()? NO RETURN
    
  }  # end of function definition WriteThreeNucleotidePeriodicity()
  
  print("Completed: Check for 3nt periodicity globally")
  
} # end ThreeNucleotidePeriodicity() function definition

# run
ThreeNucleotidePeriodicity(gene_names, dataset, hdf5file, gff_df)

#
#
# END 3NT PERIODICITY
#
# START ALL MAPPED READS
#
#

#####

DistributionOfLengthsMappedReads <- function(gene_names, dataset, hdf5file){
  
  # calculate function
  CalculateReadLengths <- function(gene_names, dataset, hdf5file){
    
    ## distribution of lengths of all mapped reads
    print("Starting: Distribution of lengths of all mapped reads")
    
    # read length-specific read counts stored as attributes of 'reads' in H5 file
    gene_sp_read_length <- lapply(gene_names, function(gene) {
      GetGeneReadLength(gene,dataset,hdf5file)
    })
    
    # sum reads of each length across all genes
    read_length_data <- data.frame(
      Length = read_range,
      Counts = gene_sp_read_length %>% 
        Reduce("+", .)
    )
    
    # return read length data
    return(read_length_data)
    
  } # end definition of function CalculateReadLengths()
  # run:
  read_length_data <- CalculateReadLengths(gene_names, dataset, hdf5file)
  
  # > str(read_length_data)
  # 'data.frame':	41 obs. of  2 variables:
  #   $ Length: int  10 11 12 13 14 15 16 17 18 19 ...
  #   $ Counts: num  0 0 0 0 0 ...
  
  # plot function
  PlotReadLengths <- function(read_length_data){
    
    # plot read lengths with counts
    read_len_plot <- ggplot(read_length_data, aes(x = Length, y = Counts)) +
      geom_bar(stat = "identity")
    
    return(read_len_plot) 
    
  } # end definition of function PlotReadLengths()
  # run:
  read_len_plot <- PlotReadLengths(read_length_data)
  # creates plot object
  
  # save pdf
  SavePlotReadLengths <- function(read_len_plot) {
    
    ggsave(read_len_plot, filename = file.path(output_dir, paste0(output_prefix, "read_lengths.pdf")))
    
    # return() # NO RETURN as writing out
    
  } # end definition of function SavePlotReadLenths()
  # to run:
  #SavePlotReadLengths(read_len_plot)
  
  # save read lengths plot and file
  WriteReadLengths <- function(read_length_data){
    
    tsv_file_path <- file.path(output_dir, paste0(output_prefix, "read_lengths.tsv"))
    
    write_provenance_header(path_to_this_script, tsv_file_path)
    
    write.table(
      read_length_data,
      file = tsv_file_path,
      append = T,
      sep = "\t",
      row = F,
      col = T,
      quote = F
    )
    
    # return() NO RETURN as writing out
    
  } # end definition of function WriteReadLengths()
  # to run: 
  #WriteReadLengths(read_length_data)
  
  print("Completed: Distribution of lengths of all mapped reads")
  
} # end of definition of function DistributionOfLengthsMappedReads()

# to run:
#DistributionOfLengthsMappedReads(gene_names, dataset, hdf5file)

DistributionOfLengthsMappedReads(gene_names, dataset, hdf5file)# load getopt to allow use of get_Rscript_filename for provenance-gathering
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
  source(here("rscripts", "provenance.R"))
  source(here("rscripts", "read_count_functions.R"))
} else {
  # Deduce file name and path using reflection as before.
  this_script <- getopt::get_Rscript_filename()
  path_to_this_script <- this_script
  source(file.path(dirname(this_script), "provenance.R")) # FLIC: consider if this can be redone using here() too
  source(file.path(dirname(this_script), "read_count_functions.R")) # FLIC: consider if this can be redone using here() too
}

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

# print provenance
print_provenance(path_to_this_script)

# read in commandline arguments
opt <- optparse::parse_args(OptionParser(option_list = option_list),
                            convert_hyphens_to_underscores=TRUE)

# attach opt list to be able to refer to variables in the list by names alone
# ie `height` rather than `women$height`
attach(opt)

print("generate_stats_figs.R running with parameters:")
opt

# prepare files, opens hdf5 file connection
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

# list of gene names taken from h5 file
gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name

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
# Following this naming convention: 
#  CalculateFunction - computational code
#  PlotFunction - plot creation
#  SaveFunction - PDF writeouts of plots
#  WriteFunction - .TSV file creation

#####

#
#
# START 3NT PERIODICITY
#
#

# big function with probable arguments
ThreeNucleotidePeriodicity <- function(gene_names, dataset, hdf5file, gff_df) {
  
  # check for 3nt periodicity
  print("Starting: Check for 3nt periodicity globally")
  
  # NOTE: Do not use mclapply when accessing H5 data
  
  CalculateThreeNucleotidePeriodicity <- function(gene_names, dataset, hdf5file, gff_df){  
    
    # get gene and position specific total counts for all read lengths
    gene_poslen_counts_5start_df <- AllGenes5StartPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hdf5file = hdf5file, gff_df = gff_df)
    
    gene_poslen_counts_3end_df <- AllGenes3EndPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hdf5file = hdf5file, gff_df = gff_df)
    
    # summarize by adding different read lengths
    gene_pos_counts_5start <- gene_poslen_counts_5start_df %>%
      group_by(Pos) %>%
      summarize(Counts = sum(Counts))
    # gives: 
    # > str(gene_pos_counts_5start)
    # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	75 obs. of  2 variables:
    #   $ Pos   : int  -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 ...
    #   $ Counts: int  285 318 307 386 291 347 840 330 475 355 ...
    
    gene_pos_counts_3end <- gene_poslen_counts_3end_df  %>%
      group_by(Pos) %>%
      summarize(Counts = sum(Counts))
    # gives: 
    # > str(gene_pos_counts_3end)
    # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	75 obs. of  2 variables:
    #   $ Pos   : int  -49 -48 -47 -46 -45 -44 -43 -42 -41 -40 ...
    #   $ Counts: int  19030 13023 50280 19458 12573 46012 19043 13282 36968 20053 ...
    
    three_nucleotide_periodicity_data <- bind_rows(
      gene_pos_counts_5start %>% mutate(End = "5'"),
      gene_pos_counts_3end %>% mutate(End = "3'")
    ) %>%
      mutate(End = factor(End, levels = c("5'", "3'")))
    # gives: 
    # > str(three_nucleotide_periodicity_data)
    # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	150 obs. of  3 variables:
    #   $ Pos   : int  -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 ...
    #   $ Counts: int  285 318 307 386 291 347 840 330 475 355 ...
    #   $ End   : Factor w/ 2 levels "5'","3'": 1 1 1 1 1 1 1 1 1 1 ...
    
    return(three_nucleotide_periodicity_data)
    
  } # end CalculateThreeNucleotidePeriodicity() definition
  
  # gives: 
  #   CalculateThreeNucleotidePeriodicity(gene_names = gene_names, dataset = dataset, hdf5file = hdf5file, gff_df = gff_df)
  #   # A tibble: 150 x 3
  #   Pos Counts End  
  #   <int>  <int> <fct>
  #     1   -24    285 5'   
  #     2   -23    318 5'   
  #     3   -22    307 5'   
  #     4   -21    386 5'   
  #     5   -20    291 5'   
  #     6   -19    347 5'   
  #     7   -18    840 5'   
  #     8   -17    330 5'   
  #     9   -16    475 5'   
  #    10   -15    355 5'   
  #   # … with 140 more rows
  
  # run the function:
  three_nucleotide_periodicity_data <- CalculateThreeNucleotidePeriodicity(gene_names = gene_names, dataset = dataset, hdf5file = hdf5file, gff_df = gff_df)
  
  # define PlotThreeNucleotidePeriodicity() function with reasonable arguments
  PlotThreeNucleotidePeriodicity <- function(three_nucleotide_periodicity_data){
    
    # Plot
    three_nucleotide_periodicity_plot <- ggplot(
      three_nucleotide_periodicity_data,
      aes(x = Pos, y = Counts)) +
      geom_line() +
      facet_wrap(~End, scales = "free") +
      labs(x = "Nucleotide Position", y = "Read counts")
    
    return(three_nucleotide_periodicity_plot) 
    
  } # end PlotThreeNucleotidePeriodicity() definition
  
  # run:
  #PlotThreeNucleotidePeriodicity(three_nucleotide_periodicity_data)
  # creates plot object
  
  three_nucleotide_periodicity_plot <- PlotThreeNucleotidePeriodicity(three_nucleotide_periodicity_data)
  
  # NOTE: repeated from inside CalculateThreeNucleotidePeriodicity() as preferred not to return multiple objects in list (hassle :S)
  gene_poslen_counts_5start_df <- AllGenes5StartPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hdf5file = hdf5file, gff_df = gff_df)
  
  
  # potentially replace/tweak plot_ribogrid() to follow StyleGuide
  PlotStartCodonRiboGrid <- function(gene_poslen_counts_5start_df){
    # function to do the ribogrid & ribogridbar plots?
    # ribogrid_5start
    start_codon_ribogrid_plot <- plot_ribogrid(gene_poslen_counts_5start_df)
    return(start_codon_ribogrid_plot) 
  } # end PlotStartCodonRiboGrid() definition
  
  # run: 
  start_codon_ribogrid_plot <- PlotStartCodonRiboGrid(gene_poslen_counts_5start_df)
  # creates plot object
  
  
  SaveStartCodonRiboGrid <- function(start_codon_ribogrid_plot){
    # function to do the ribogrid & ribogridbar plots?
    # ribogrid_5start
    start_codon_ribogrid_plot %>%
      ggsave(
        filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogrid.pdf")),
        width = 6, height = 3
      )
    #return() # no return as writing-out
  } # end SaveStartCodonRiboGrid() definition
  
  # run: 
  SaveStartCodonRiboGrid(start_codon_ribogrid_plot)
  
  
  PlotStartCodonRiboGridBar <- function(gene_poslen_counts_5start_df){  
    start_codon_ribogrid_bar_plot <- barplot_ribogrid(gene_poslen_counts_5start_df)
    return(start_codon_ribogrid_bar_plot) 
  } # end PlotStartCodonRiboGridBar() definition
  
  # run: 
  start_codon_ribogrid_bar_plot <- PlotStartCodonRiboGridBar(gene_poslen_counts_5start_df)
  # creates plot object
  
  
  SaveStartCodonRiboGridBar <- function(start_codon_ribogrid_bar_plot){  
    start_codon_ribogrid_bar_plot %>%  
      ggsave(
        filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogridbar.pdf")),
        width = 6, height = 5
      )
    #return() # no return as writing-out
  } # end SaveStartCodonRiboGridBar() definition
  
  # run: 
  SaveStartCodonRiboGrid(start_codon_ribogrid_bar_plot)
  
  
  SavePlotThreeNucleotidePeriodicity <- function() {
    
    # Save plot and file
    ggsave(
      three_nucleotide_periodicity_plot, 
      filename = file.path(output_dir, paste0(output_prefix, "3nt_periodicity.pdf"))
    )
    # return() # NO RETURN as writing out
    
  } # end of function definition SavePlotThreeNucleotidePeriodicity()
  
  
  WriteThreeNucleotidePeriodicity <- function() {
    
    tsv_file_path <- file.path(output_dir, paste0(output_prefix, "3nt_periodicity.tsv"))
    
    write_provenance_header(path_to_this_script, tsv_file_path)
    
    write.table(
      three_nucleotide_periodicity_data,
      file = tsv_file_path,
      append = T,
      sep = "\t",
      row = F,
      col = T,
      quote = F)
    
    # return()? NO RETURN
    
  }  # end of function definition WriteThreeNucleotidePeriodicity()
  
  print("Completed: Check for 3nt periodicity globally")
  
} # end ThreeNucleotidePeriodicity() function definition
# run
ThreeNucleotidePeriodicity(gene_names, dataset, hdf5file, gff_df)

#
#
# END 3NT PERIODICITY
#
# START ALL MAPPED READS
#
#

DistributionOfLengthsMappedReads <- function(gene_names, dataset, hdf5file){
  
  # calculate function
  CalculateReadLengths <- function(gene_names, dataset, hdf5file){
    
    ## distribution of lengths of all mapped reads
    print("Starting: Distribution of lengths of all mapped reads")
    
    # read length-specific read counts stored as attributes of 'reads' in H5 file
    gene_sp_read_length <- lapply(gene_names, function(gene) {
      GetGeneReadLength(gene,dataset,hdf5file)
    })
    
    # sum reads of each length across all genes
    read_length_data <- data.frame(
      Length = read_range,
      Counts = gene_sp_read_length %>% 
        Reduce("+", .)
    )
    
    # return read length data
    return(read_length_data)
    
  } # end definition of function CalculateReadLengths()
  # run:
  read_length_data <- CalculateReadLengths(gene_names, dataset, hdf5file)
  
  # > str(read_length_data)
  # 'data.frame':	41 obs. of  2 variables:
  #   $ Length: int  10 11 12 13 14 15 16 17 18 19 ...
  #   $ Counts: num  0 0 0 0 0 ...
  
  # plot function
  PlotReadLengths <- function(read_length_data){
    
    # plot read lengths with counts
    read_len_plot <- ggplot(read_length_data, aes(x = Length, y = Counts)) +
      geom_bar(stat = "identity")
    
    return(read_len_plot) 
    
  } # end definition of function PlotReadLengths()
  # run:
  read_len_plot <- PlotReadLengths(read_length_data)
  # creates plot object
  
  # save pdf
  SavePlotReadLengths <- function(read_len_plot) {
    
    ggsave(read_len_plot, filename = file.path(output_dir, paste0(output_prefix, "read_lengths.pdf")))
    
    # return() # NO RETURN as writing out
    
  } # end definition of function SavePlotReadLenths()
  # to run:
  #SavePlotReadLengths(read_len_plot)
  
  # save read lengths plot and file
  WriteReadLengths <- function(read_length_data){
    
    tsv_file_path <- file.path(output_dir, paste0(output_prefix, "read_lengths.tsv"))
    
    write_provenance_header(path_to_this_script, tsv_file_path)
    
    write.table(
      read_length_data,
      file = tsv_file_path,
      append = T,
      sep = "\t",
      row = F,
      col = T,
      quote = F
    )
    
    # return() NO RETURN as writing out
    
  } # end definition of function WriteReadLengths()
  # to run: 
  #WriteReadLengths(read_length_data)
  
  print("Completed: Distribution of lengths of all mapped reads")
  
} # end of definition of function DistributionOfLengthsMappedReads()

# to run:
#DistributionOfLengthsMappedReads(gene_names, dataset, hdf5file)

DistributionOfLengthsMappedReads(gene_names, dataset, hdf5file)

#
#
# END ALL MAPPED READS
#
# START BIASES IN NUCLEOTIDE COMPOSITION
#
# 


if (!do_pos_sp_nt_freq) {
  
  print("NOT calculating position-specific nucleotide frequency - reason: do_pos_sp_nt_freq parameter set to FALSE")
  
}# TODO: FLIC FIGURE OUT IF USING ELSE OR NOT # else { # if do_pos_sp_nt_freq parameter set to true, calculate position-specific nucleotide frequency:

# big function with probable arguments
BiasesInNucleotideCompositionAlongMappedReadLengths <- function(gene_names, dataset, hdf5file, read_range, min_read_length) {
  
  print("Starting: Biases in nucleotide composition along mapped read lengths")
  
  ## CALCULATE Biases In Nucleotide Composition Along Mapped Read Lengths  
  
  CalculateBiasesInNucleotideComposition <- function(gene_names, dataset, hdf5file, read_range, min_read_length){
    
    # This is in a conditional loop because it fails for some inputs
    # and has not been debugged. Needs to be rewritten in tidyverse
    all_out <- c()
    for (lid in seq_len(length(read_range))) { # TODO: WHAT IS LID
      out <- lapply(gene_names, function(x) { # TODO: fix variablename 'out' (also check below for uses)
        # For each read length convert reads to IRanges
        GetNTReadPosition(gene = as.character(x), 
                          dataset = dataset,
                          hdf5file = hdf5file, 
                          lid = lid, min_read_length = min_read_length)
      })
      names(out) <- gene_names
      
      # Get position-specific nucleotide counts for reads in each frame
      
      # # TODO: test this function
      # CalculatePositionSpecificNucleotideCountsByFrame <- function(gene_names, cframe, lid, num_processes){
      #   mcapply(gene_names, function(gene){
      #     cons_mat(gene = gene, pos_IR = out[[gene]], cframe = cframe, lid = lid)
      #   }, mc.cores = num_processes)
      # } # end of function definition CalculatePositionSpecificNucleotideCountsByFrame()
      # test_frame0 <- CalculatePositionSpecificNucleotideCountsByFrame(gene_names, cframe=0, lid, num_processes)
      
      # frame 0
      fr0 <- mclapply(gene_names, function(gene) {
        cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 0, lid = lid)
      }, mc.cores = num_processes)
      allfr0 <- do.call(rbind, fr0)
      
      # frame 1
      fr1 <- mclapply(gene_names, function(gene) {
        cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 1, lid = lid)
      }, mc.cores = num_processes)
      allfr1 <- do.call(rbind, fr1)
      
      # frame 2
      fr2 <- mclapply(gene_names, function(gene) {
        cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 2, lid = lid)
      }, mc.cores = num_processes)
      allfr2 <- do.call(rbind, fr2)
      
      # Get position-specific freq for all nucleotides
      cnt_fr0 <- signif(comb_freq(allfr0), 3)
      cnt_fr1 <- signif(comb_freq(allfr1), 3)
      cnt_fr2 <- signif(comb_freq(allfr2), 3)
      
      output <- data.frame(rbind(cnt_fr0, cnt_fr1, cnt_fr2)) # TODO: replace variablename 'output'
      all_out <- rbind(all_out, output)
    } # end of the for() loop
    
    # TODO:
    print("finished for loop")
    
    # Prepare variables for output file
    Length <- unlist(lapply(read_range, function(x) {
      rep(x, x * 3)
    }))
    Position <- unlist(lapply(read_range, function(x) {
      rep(1:x, 3)
    }))
    Frame <- unlist(lapply(read_range, function(x) {
      rep(0:2, each = x)
    }))
    
    # TODO:
    print("finished prepping L, R, F variables")
    
    all_out <- cbind(Length, Position, Frame, all_out)
    all_out[is.na(all_out)] <- 0
    
    return(all_out) # TODO: ensure this is correct; rename the returned output
    
    # TODO:
    print("finished")
    
  } # end of function definition of CalculateBiasesInNucleotideComposition()
  # run: 
  #CalculateBiasesInNucleotideComposition(gene_names, dataset, hdf5file, read_range, min_read_length)
  
  
  all_out <- CalculateBiasesInNucleotideComposition(gene_names, dataset, hdf5file, read_range, min_read_length)
  
  ## WRITE DATA Biases In Nucleotide Composition Along Mapped Read Lengths
  WriteBiasesInNucleotideComposition <- function(all_out){
    # save file
    tsv_file_path <- file.path(output_dir, paste0(output_prefix, "pos_sp_nt_freq.tsv"))
    write_provenance_header(path_to_this_script, tsv_file_path)
    write.table(all_out, file = tsv_file_path, append = T, sep = "\t", row = F, col = T, quote = F)
  } # end of function definition: WriteBiasesInNucleotideComposition()
  # run: 
  #WriteBiasesInNucleotideComposition()
  
  print("Completed nucleotide composition bias table")
  
} # end definition of function: BiasesInNucleotideCompositionAlongMappedReadLengths()

# to run:
#BiasesInNucleotideCompositionAlongMappedReadLengths(gene_names, dataset, hdf5file, read_range, min_read_length)

BiasesInNucleotideCompositionAlongMappedReadLengths(gene_names, dataset, hdf5file, read_range, min_read_length)

#
#
# END BIASES IN NUCLEOTIDE COMPOSITION
#
# START CALCULATE READ FRAME FOR EVERY ORF
#
# 

## calculate read frame for every annotated ORF

if (!is.na(asite_disp_length_file)) {
  print("Starting: Check for 3nt periodicity (frame) by Gene")
  asite_disp_length <- readr::read_tsv(asite_disp_length_file,
                                       comment = "#"
  )
  # TODO: wrap in function
  gene_read_frames <- gff_df %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::select(gene = seqnames, left = start, right = end) %>%
    purrr::pmap_dfr(GetGeneReadFrame,
                    hdf5file = hdf5file,
                    dataset = dataset,
                    min_read_length = min_read_length,
                    asite_disp_length = asite_disp_length
    )
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "3ntframe_bygene.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    gene_read_frames,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
  
  gene_read_frame_plot <- gene_read_frames %>%
    filter(Ct_fr0 + Ct_fr1 + Ct_fr2 > count_threshold) %>%
    BoxplotReadFrameProportion()
  
  # save read lengths plot and file
  ggsave(gene_read_frame_plot,
         filename = file.path(output_dir, paste0(output_prefix, "3ntframe_propbygene.pdf")),
         width = 3, height = 3
  )
  
  print("Completed: Check for 3nt periodicity (frame) by Gene")
}

#
#
# END CALCULATE READ FRAME FOR EVERY ORF
#
# START RPF POSITION SPECIFIC DISTRIBUTION OF READS
#
#

#####
## position specific distribution of reads

print("Starting: Position specific distribution of reads")

# @FlicAnderson: split this if chunk up
# For RPF datasets, generate codon-based position-specific reads
if (rpf & !is.na(asite_disp_length_file)) {
  # Get codon-based position-specific reads for each gene, in a tibble
  reads_per_codon_etc <- tibble(gene=gene_names) %>%
    mutate(CountPerCodon = map(gene, ~GetGeneCodonPosReads1dsnap(
      .,
      dataset,
      hdf5file,
      left = GetCDS5start(., gff_df),
      right = GetCDS3end(., gff_df),
      min_read_length = min_read_length,
      asite_disp_length = asite_disp_length,
    ) ),
    NormCtperCodon = map(CountPerCodon, ~NormByMean(.) ),
    SumAsiteCt = map(CountPerCodon,~sum(.)) %>% unlist,
    LengthCodons = map(CountPerCodon,~length(.)) %>% unlist
    )
  
  # metagene mean normalized read counts for 5' end positions 1 to length_threshold
  pos_sp_rpf_norm_reads_5end <- 
    reads_per_codon_etc %>%
    # select genes over count_threshold and length_threshold
    filter(SumAsiteCt > count_threshold, LengthCodons > length_threshold) %>%
    # Extract from codon positions 1 to length_threshold
    transmute(NormCtperCodon5end = 
                map(CountPerCodon,
                    ~extract(.,seq_len(length_threshold)))) %>%
    # remove NormCts as vector and (3 lines later) make these rows of a matrix
    dplyr::pull(NormCtperCodon5end) %>%
    unlist %>% 
    matrix(byrow=T,ncol=length_threshold) %>%
    # calculate column means and standard deviations in tidy data frame
    # @FlicAnderson: Why is this in {} brackets?
    {tibble(Position = seq_len(length_threshold),
            Mean = colMeans(.) %>% as.numeric,
            SD   = apply(.,2,sd),
            End  = "5'")}
  
  # metagene mean normalized read counts for 3' end stop-length_threshold to stop
  pos_sp_rpf_norm_reads_3end <- 
    reads_per_codon_etc %>%
    # select genes over count_threshold and length_threshold
    filter(SumAsiteCt > count_threshold, LengthCodons > length_threshold) %>%
    # Extract from last length_threshold codon positions
    transmute(NormCtperCodon3end = 
                map2(CountPerCodon,LengthCodons,
                     ~extract(.x, .y + 1 - seq_len(length_threshold)))) %>%              
    # remove NormCts as vector and (3 lines later) make these rows of a matrix
    dplyr::pull(NormCtperCodon3end) %>%
    unlist %>% 
    matrix(byrow=T,ncol=length_threshold)  %>% 
    # calculate column means and standard deviations in tidy data frame
    {tibble(Position = seq_len(length_threshold),
            Mean = colMeans(.),
            SD   = apply(.,2,sd),
            End  = "3'")}
  
  pos_sp_rpf_norm_reads <- bind_rows(pos_sp_rpf_norm_reads_5end,
                                     pos_sp_rpf_norm_reads_3end)
  
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
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    pos_sp_rpf_norm_reads %>%
      # @FlicAnderson: CHECK THIS! the numbers look quite different to a recent vignette package output..?
      # mutate_if() mutates if the data is numeric, rounds values to specified number digits
      # Is 4 significant figures in mutate_if() enough? Leads to (0, 1, 1.5, 0.5, 1.414)
      mutate_if(is.numeric, signif, digits=4), 
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
}

#
#
# END RPF POSITION SPECIFIC DISTRIBUTION OF READS
#
# START MRNA POSITION SPECIFIC DISTRIBUTION OF READS
#
#

# for mRNA datasets, generate nt-based position-specific reads
if (!rpf) {
  # TODO: rewrite along the lines of rpf above
  warning("Warning: nt-based position-specific reads for non-RPF datasets not tested")
  # create empty matrix to store position-specific read counts
  out5p <- matrix(NA, nrow = length(gene_names), ncol = 1500) # 5'
  out3p <- matrix(NA, nrow = length(gene_names), ncol = 1500) # 3'
  
  # TODO
  out <- lapply(gene_names, function(gene) {
    GetMRNACoverage(
      gene,
      dataset,
      left = (buffer - 49),
      right = (buffer - 3),
      min_read_length = min_read_length,
      read_range,
      buffer
    )
  })
  names(out) <- gene_names
  
  cc <- 1
  for (gene in gene_names) {
    tmp <- out[[gene]]
    # only consider genes with at least 1 mapped read along its CDS
    if (sum(tmp) > 0) {
      tmp <- tmp / mean(tmp)
      if (length(tmp) > 1500) {
        out5p[cc, ] <- tmp[1:1500]
        out3p[cc, ] <- rev(tmp)[1:1500]
      } else {
        out5p[cc, 1:length(tmp)] <- tmp
        out3p[cc, 1:length(tmp)] <- rev(tmp)
      }
    }
    cc <- cc + 1
  }
  
  # estimate position-specific mean and std error of mapped read counts
  m5p <- signif(apply(out5p, 2, mean, na.rm = T), 4)
  m3p <- signif(apply(out3p, 2, mean, na.rm = T), 4)
  s5p <- signif(apply(out5p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)
  s3p <- signif(apply(out3p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)
  
  # normalize reads to last 150 nts of the 1500-nt window.
  # this allows easy comparison between datasets
  s5p <- s5p / mean(m5p[1350:1500])
  s3p <- s3p / mean(m3p[1350:1500])
  m5p <- m5p / mean(m5p[1350:1500])
  m3p <- m3p / mean(m3p[1350:1500])
  
  # create a dataframe to store the output for plots/analyses
  pos_sp_mrna_norm_coverage <- data.frame(
    Position = c(1:1500, 0:-1499),
    Mean = c(m5p, m3p),
    SD = c(s5p, s3p),
    End = factor(rep(c("5'", "3'"), each = 1500), levels = c("5'", "3'"))
  )
  
  # plot
  pos_sp_mrna_norm_coverage_plot <- ggplot(pos_sp_mrna_norm_coverage, aes(Position, Mean, col = End)) +
    geom_line() +
    facet_grid(~End, scales = "free") +
    guides(col = FALSE)
  
  # Save plot and file
  ggsave(pos_sp_mrna_norm_coverage_plot, filename = file.path(output_dir, paste0(output_prefix, "pos_sp_mrna_norm_coverage.pdf")))
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "pos_sp_mrns_norm_coverage.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    pos_sp_mrna_norm_coverage,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
}

print("Completed: Position specific distribution of reads")


#
#
# END MRNA POSITION SPECIFIC DISTRIBUTION OF READS
#
# START TPMS OF GENES
#
#


## Calculate TPMs of genes

print("Starting: Calculate TPMs of genes")

# calculate transcripts per million (TPM)
# TODO: make this a single function
gene_sp_reads <- sapply(gene_names, GetGeneReadsTotal, dataset, hdf5file)
reads_per_b <- sapply(gene_names, GetGeneReadDensity, dataset, hdf5file)

tpms <- data.frame(
  ORF = gene_names,
  readcount = gene_sp_reads,
  rpb = reads_per_b,
  tpm = reads_per_b * 1e6 / sum(reads_per_b)
)

# write out to *_tpms.tsv
tsv_file_path <- file.path(output_dir, paste0(output_prefix, "tpms.tsv"))
write_provenance_header(path_to_this_script, tsv_file_path)
write.table(
  tpms,
  file = tsv_file_path,
  append = T,
  sep = "\t",
  row = F,
  col = T,
  quote = F
)

#
#
# END TPMS OF GENES
#
# START TPMS CORRELATIONS WITH FEATURES
#
#

## Correlations between TPMs of genes with their sequence-based features

# Correlate TPMs of genes with sequence-based features, skip if missing features_file
if (!is.na(features_file)) {
  print("Starting: Correlations between TPMs of genes with their sequence-based features")
  
  features <- read.table(features_file, h = T)
  
  # Prepare data for plot
  # Consider only genes with at least count_threshold mapped reads
  features_plot_data <- merge(features, tpms, by = "ORF") %>%
    filter(readcount >= count_threshold, !is.na(ORF)) %>%
    select(-readcount, -rpb) %>%
    gather(Feature, Value, -ORF, -tpm)
  
  features_plot <- ggplot(features_plot_data, aes(x = tpm, y = Value)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~Feature, scales = "free") +
    scale_x_log10() +
    geom_smooth(method = "lm") +
    xlab("TPM (transcripts per million)")
  
  # Save plot and file
  ggsave(features_plot, filename = file.path(output_dir, paste0(output_prefix, "features.pdf")))
  
  print("Completed: Correlations between TPMs of genes with their sequence-based features")
} else {
  print("Skipped: Correlations between TPMs of genes with their sequence-based features - features_file.tsv not provided")
}

## Codon-specific ribosome densities for correlations with tRNAs

# Codon-specific ribosome density for tRNA correlation; skip if missing t_rna_file & codon_positions_file
if (!is.na(t_rna_file) & !is.na(codon_positions_file)) {
  print("Starting: Codon-specific ribosome densities for correlations with tRNAs")
  # Only for RPF datasets
  if (rpf) {
    # TODO: This section needs attention. Can be refactored analogously to reads_per_codon_etc
    # Needs separate calculation of per-codon normalized reads
    # WAITING: we want new format of codon_pos from @john-s-f before editing this chunk
    t_rna_df <- read.table(t_rna_file, h = T) # Read in yeast tRNA estimates
    load(codon_positions_file) # Position of codons in each gene (numbering ignores first 200 codons)
    # Reads in an object named "codon_pos"
    out <- lapply(gene_names, function(gene) {
      # From "Position specific distribution of reads" plot
      GetCodonPositionReads(gene, dataset, hdf5file, 
                            left = (buffer - 15), right = (buffer + 11), 
                            min_read_length = min_read_length)
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
    
    
    # TODO: figure this out
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
    
    # TODO: this is a misnomer. Can calculate A/P/E-site norm density without t_rna_df.
    # Should replace t_rna_df with an argument "codon_features_file"
    # then plot against features in that, analogosly to (gene) to "features_file" above
    cod_dens_tRNA <- cbind(t_rna_df, A, P, E)
    
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
    write_provenance_header(path_to_this_script, tsv_file_path)
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

#
#
# END ALL MAPPED READS
#
# START BIASES IN NUCLEOTIDE COMPOSITION
#
# 

#####
## biases in nucleotide composition along mapped read lengths

## Calculate Biases In Nucleotide Composition Along Mapped Read

print("Starting: Biases in nucleotide composition along mapped read lengths")

if (do_pos_sp_nt_freq) {
  # This is in a conditional loop because it fails for some inputs
  # and has not been debugged. Needs to be rewritten in tidyverse
  all_out <- c()
  for (lid in seq_len(length(read_range))) {
    out <- lapply(gene_names, function(x) {
      # For each read length convert reads to IRanges
      GetNTReadPosition(gene = as.character(x), 
                        dataset = dataset,
                        hdf5file = hdf5file, 
                        lid = lid, min_read_length = min_read_length)
    })
    names(out) <- gene_names

    # Get position-specific nucleotide counts for reads in each frame
    fr0 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 0, lid = lid)
    }, mc.cores = num_processes)
    allfr0 <- do.call(rbind, fr0)
    fr1 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 1, lid = lid)
    }, mc.cores = num_processes)
    allfr1 <- do.call(rbind, fr1)
    fr2 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 2, lid = lid)
    }, mc.cores = num_processes)
    allfr2 <- do.call(rbind, fr2)

    # Get position-specific freq for all nts
    cnt_fr0 <- signif(comb_freq(allfr0), 3)
    cnt_fr1 <- signif(comb_freq(allfr1), 3)
    cnt_fr2 <- signif(comb_freq(allfr2), 3)

    output <- data.frame(rbind(cnt_fr0, cnt_fr1, cnt_fr2))
    all_out <- rbind(all_out, output)
  }

  # Prepare variables for output file
  Length <- unlist(lapply(read_range, function(x) {
    rep(x, x * 3)
  }))
  Position <- unlist(lapply(read_range, function(x) {
    rep(1:x, 3)
  }))
  Frame <- unlist(lapply(read_range, function(x) {
    rep(0:2, each = x)
  }))

  all_out <- cbind(Length, Position, Frame, all_out)
  all_out[is.na(all_out)] <- 0

  # save file
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "pos_sp_nt_freq.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(all_out, file = tsv_file_path, append = T, sep = "\t", row = F, col = T, quote = F)

  print("Completed nucleotide composition bias table")
}

#
#
# END BIASES IN NUCLEOTIDE COMPOSITION
#
# START CALCULATE READ FRAME FOR EVERY ORF
#
# 

## calculate read frame for every annotated ORF

if (!is.na(asite_disp_length_file)) {
  print("Starting: Check for 3nt periodicity (frame) by Gene")
  asite_disp_length <- readr::read_tsv(asite_disp_length_file,
    comment = "#"
  )
  # TODO: wrap in function
  gene_read_frames <- gff_df %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::select(gene = seqnames, left = start, right = end) %>%
    purrr::pmap_dfr(GetGeneReadFrame,
      hdf5file = hdf5file,
      dataset = dataset,
      min_read_length = min_read_length,
      asite_disp_length = asite_disp_length
    )
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "3ntframe_bygene.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    gene_read_frames,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )

  gene_read_frame_plot <- gene_read_frames %>%
    filter(Ct_fr0 + Ct_fr1 + Ct_fr2 > count_threshold) %>%
    BoxplotReadFrameProportion()

  # save read lengths plot and file
  ggsave(gene_read_frame_plot,
    filename = file.path(output_dir, paste0(output_prefix, "3ntframe_propbygene.pdf")),
    width = 3, height = 3
  )

  print("Completed: Check for 3nt periodicity (frame) by Gene")
}

#
#
# END CALCULATE READ FRAME FOR EVERY ORF
#
# START RPF POSITION SPECIFIC DISTRIBUTION OF READS
#
#

#####
## position specific distribution of reads

print("Starting: Position specific distribution of reads")

# @FlicAnderson: split this if chunk up
# For RPF datasets, generate codon-based position-specific reads
if (rpf & !is.na(asite_disp_length_file)) {
  # Get codon-based position-specific reads for each gene, in a tibble
  reads_per_codon_etc <- tibble(gene=gene_names) %>%
	mutate(CountPerCodon = map(gene, ~GetGeneCodonPosReads1dsnap(
		.,
		dataset,
		hdf5file,
		left = GetCDS5start(., gff_df),
		right = GetCDS3end(., gff_df),
		min_read_length = min_read_length,
		asite_disp_length = asite_disp_length,
		) ),
		NormCtperCodon = map(CountPerCodon, ~NormByMean(.) ),
		SumAsiteCt = map(CountPerCodon,~sum(.)) %>% unlist,
		LengthCodons = map(CountPerCodon,~length(.)) %>% unlist
	)
  
  # metagene mean normalized read counts for 5' end positions 1 to length_threshold
  pos_sp_rpf_norm_reads_5end <- 
    reads_per_codon_etc %>%
    # select genes over count_threshold and length_threshold
    filter(SumAsiteCt > count_threshold, LengthCodons > length_threshold) %>%
    # Extract from codon positions 1 to length_threshold
    transmute(NormCtperCodon5end = 
                map(CountPerCodon,
                    ~extract(.,seq_len(length_threshold)))) %>%
    # remove NormCts as vector and (3 lines later) make these rows of a matrix
    dplyr::pull(NormCtperCodon5end) %>%
    unlist %>% 
    matrix(byrow=T,ncol=length_threshold) %>%
    # calculate column means and standard deviations in tidy data frame
    # @FlicAnderson: Why is this in {} brackets?
    {tibble(Position = seq_len(length_threshold),
            Mean = colMeans(.) %>% as.numeric,
            SD   = apply(.,2,sd),
            End  = "5'")}
  
  # metagene mean normalized read counts for 3' end stop-length_threshold to stop
  pos_sp_rpf_norm_reads_3end <- 
    reads_per_codon_etc %>%
        # select genes over count_threshold and length_threshold
    filter(SumAsiteCt > count_threshold, LengthCodons > length_threshold) %>%
    # Extract from last length_threshold codon positions
    transmute(NormCtperCodon3end = 
              map2(CountPerCodon,LengthCodons,
                   ~extract(.x, .y + 1 - seq_len(length_threshold)))) %>%              
    # remove NormCts as vector and (3 lines later) make these rows of a matrix
    dplyr::pull(NormCtperCodon3end) %>%
    unlist %>% 
    matrix(byrow=T,ncol=length_threshold)  %>% 
    # calculate column means and standard deviations in tidy data frame
    {tibble(Position = seq_len(length_threshold),
            Mean = colMeans(.),
            SD   = apply(.,2,sd),
            End  = "3'")}
  
  pos_sp_rpf_norm_reads <- bind_rows(pos_sp_rpf_norm_reads_5end,
                                     pos_sp_rpf_norm_reads_3end)
  
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
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    pos_sp_rpf_norm_reads %>%
      # @FlicAnderson: CHECK THIS! the numbers look quite different to a recent vignette package output..?
      # mutate_if() mutates if the data is numeric, rounds values to specified number digits
      # Is 4 significant figures in mutate_if() enough? Leads to (0, 1, 1.5, 0.5, 1.414)
      mutate_if(is.numeric, signif, digits=4), 
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
}

#
#
# END RPF POSITION SPECIFIC DISTRIBUTION OF READS
#
# START MRNA POSITION SPECIFIC DISTRIBUTION OF READS
#
#

# for mRNA datasets, generate nt-based position-specific reads
if (!rpf) {
  # TODO: rewrite along the lines of rpf above
  warning("Warning: nt-based position-specific reads for non-RPF datasets not tested")
  # create empty matrix to store position-specific read counts
  out5p <- matrix(NA, nrow = length(gene_names), ncol = 1500) # 5'
  out3p <- matrix(NA, nrow = length(gene_names), ncol = 1500) # 3'
  
  # TODO
  out <- lapply(gene_names, function(gene) {
    GetMRNACoverage(
      gene,
      dataset,
      left = (buffer - 49),
      right = (buffer - 3),
      min_read_length = min_read_length,
      read_range,
      buffer
      )
  })
  names(out) <- gene_names

  cc <- 1
  for (gene in gene_names) {
    tmp <- out[[gene]]
    # only consider genes with at least 1 mapped read along its CDS
    if (sum(tmp) > 0) {
      tmp <- tmp / mean(tmp)
      if (length(tmp) > 1500) {
        out5p[cc, ] <- tmp[1:1500]
        out3p[cc, ] <- rev(tmp)[1:1500]
      } else {
        out5p[cc, 1:length(tmp)] <- tmp
        out3p[cc, 1:length(tmp)] <- rev(tmp)
      }
    }
    cc <- cc + 1
  }

  # estimate position-specific mean and std error of mapped read counts
  m5p <- signif(apply(out5p, 2, mean, na.rm = T), 4)
  m3p <- signif(apply(out3p, 2, mean, na.rm = T), 4)
  s5p <- signif(apply(out5p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)
  s3p <- signif(apply(out3p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)

  # normalize reads to last 150 nts of the 1500-nt window.
  # this allows easy comparison between datasets
  s5p <- s5p / mean(m5p[1350:1500])
  s3p <- s3p / mean(m3p[1350:1500])
  m5p <- m5p / mean(m5p[1350:1500])
  m3p <- m3p / mean(m3p[1350:1500])

  # create a dataframe to store the output for plots/analyses
  pos_sp_mrna_norm_coverage <- data.frame(
    Position = c(1:1500, 0:-1499),
    Mean = c(m5p, m3p),
    SD = c(s5p, s3p),
    End = factor(rep(c("5'", "3'"), each = 1500), levels = c("5'", "3'"))
  )

  # plot
  pos_sp_mrna_norm_coverage_plot <- ggplot(pos_sp_mrna_norm_coverage, aes(Position, Mean, col = End)) +
    geom_line() +
    facet_grid(~End, scales = "free") +
    guides(col = FALSE)

  # Save plot and file
  ggsave(pos_sp_mrna_norm_coverage_plot, filename = file.path(output_dir, paste0(output_prefix, "pos_sp_mrna_norm_coverage.pdf")))
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "pos_sp_mrns_norm_coverage.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    pos_sp_mrna_norm_coverage,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
}

print("Completed: Position specific distribution of reads")


#
#
# END MRNA POSITION SPECIFIC DISTRIBUTION OF READS
#
# START TPMS OF GENES
#
#


## Calculate TPMs of genes

print("Starting: Calculate TPMs of genes")

# calculate transcripts per million (TPM)
# TODO: make this a single function
gene_sp_reads <- sapply(gene_names, GetGeneReadsTotal, dataset, hdf5file)
reads_per_b <- sapply(gene_names, GetGeneReadDensity, dataset, hdf5file)

tpms <- data.frame(
  ORF = gene_names,
  readcount = gene_sp_reads,
  rpb = reads_per_b,
  tpm = reads_per_b * 1e6 / sum(reads_per_b)
)

# write out to *_tpms.tsv
tsv_file_path <- file.path(output_dir, paste0(output_prefix, "tpms.tsv"))
write_provenance_header(path_to_this_script, tsv_file_path)
write.table(
  tpms,
  file = tsv_file_path,
  append = T,
  sep = "\t",
  row = F,
  col = T,
  quote = F
)

#
#
# END TPMS OF GENES
#
# START TPMS CORRELATIONS WITH FEATURES
#
#

## Correlations between TPMs of genes with their sequence-based features

# Correlate TPMs of genes with sequence-based features, skip if missing features_file
if (!is.na(features_file)) {
  print("Starting: Correlations between TPMs of genes with their sequence-based features")

  features <- read.table(features_file, h = T)

  # Prepare data for plot
  # Consider only genes with at least count_threshold mapped reads
  features_plot_data <- merge(features, tpms, by = "ORF") %>%
    filter(readcount >= count_threshold, !is.na(ORF)) %>%
    select(-readcount, -rpb) %>%
    gather(Feature, Value, -ORF, -tpm)

  features_plot <- ggplot(features_plot_data, aes(x = tpm, y = Value)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~Feature, scales = "free") +
    scale_x_log10() +
    geom_smooth(method = "lm") +
    xlab("TPM (transcripts per million)")

  # Save plot and file
  ggsave(features_plot, filename = file.path(output_dir, paste0(output_prefix, "features.pdf")))

  print("Completed: Correlations between TPMs of genes with their sequence-based features")
} else {
  print("Skipped: Correlations between TPMs of genes with their sequence-based features - features_file.tsv not provided")
}

## Codon-specific ribosome densities for correlations with tRNAs

# Codon-specific ribosome density for tRNA correlation; skip if missing t_rna_file & codon_positions_file
if (!is.na(t_rna_file) & !is.na(codon_positions_file)) {
  print("Starting: Codon-specific ribosome densities for correlations with tRNAs")
  # Only for RPF datasets
  if (rpf) {
    # TODO: This section needs attention. Can be refactored analogously to reads_per_codon_etc
    # Needs separate calculation of per-codon normalized reads
    # WAITING: we want new format of codon_pos from @john-s-f before editing this chunk
    t_rna_df <- read.table(t_rna_file, h = T) # Read in yeast tRNA estimates
    load(codon_positions_file) # Position of codons in each gene (numbering ignores first 200 codons)
    # Reads in an object named "codon_pos"
    out <- lapply(gene_names, function(gene) {
      # From "Position specific distribution of reads" plot
      GetCodonPositionReads(gene, dataset, hdf5file, 
                            left = (buffer - 15), right = (buffer + 11), 
                            min_read_length = min_read_length)
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
    
    
    # TODO: figure this out
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
    
    # TODO: this is a misnomer. Can calculate A/P/E-site norm density without t_rna_df.
    # Should replace t_rna_df with an argument "codon_features_file"
    # then plot against features in that, analogosly to (gene) to "features_file" above
    cod_dens_tRNA <- cbind(t_rna_df, A, P, E)

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
    write_provenance_header(path_to_this_script, tsv_file_path)
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