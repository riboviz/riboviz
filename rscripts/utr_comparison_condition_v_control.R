# create a scatter plot comparing the how the UTR:CDS use changes between a condition and a control

print('Starting')

suppressMessages(library(here))
suppressMessages(library(optparse))
source(here::here('rscripts','read_count_functions.R'))

option_list <- list(make_option(c("-i", "--input"),
                                type = "character",
                                help = "Path input to h5 file"),
                    make_option(c("-c", "--compare"),
                                type = "character",
                                help = "path to H5 to be compared with input file. Ideally a control file"),
                    make_option(c("-d", "--dataset"),
                                type = "character",
                                help = "Name of the dataset being studied"),
                    make_option(c("-g", "--gff"),
                                type = "character",
                                help = "Path t the GFF3 file of the organism being studied"),
                    make_option(c("-o", "--output"),
                                type = "character",
                                help = "Path to output directory",
                                default = '.'),
                    make_option("--gene",
                                type = "character",
                                help = "Gene of interest to be highlighted in plots",
                                default = NULL),
                    make_option(c("-t", "--treatment"),
                                type = "character",
                                help ='Treatment used on the cells, for labelling of plots',
                                default = NULL),
                    make_option("--read_threshold",
                                type = "numeric",
                                help = "The minimum reads per base value needed in both samples for a gene to be included in the plot",
                                default = 0.02)
)

opt <- optparse::parse_args(OptionParser(option_list = option_list))



h5_file_path <- opt$input
comparison_file <- opt$compare
gff_in <- opt$gff
Gene_of_interest <- opt$gene
dataset <- opt$dataset
output_dir <- opt$output
treatment <- opt$treatment
read_threshold <- opt$read_threshold

## functions 


# GeneFeatureTotalCountsPerBase takes a gene, splits it into features (5UTR and CDS),
# and returns the reads per base in the features of that gene in the form of a tibble

GeneFeatureTotalCountsPerBase <- function(gene = gene, dataset = dataset, hd_file = hd_file, startpos = 1){
  
  #get the start and stop codons
  start <- rhdf5::h5readAttributes(hd_file, paste0('/', gene,'/', dataset, '/reads'))[['start_codon_pos']]
  stop <- rhdf5::h5readAttributes(hd_file, paste0('/', gene,'/', dataset, '/reads'))[['stop_codon_pos']]
  
  #make a matric to be processed for each gene
  # turn matric into a tidy datamatrix 
  tmp_tidy <- TidyDatamatrix(
    GetGeneDatamatrix(gene = gene, dataset = dataset, hd_file = hd_file), 
    startpos = 1
  )
  
  # filter the tmp_tidy matric into two new matrixes, on with 5UTR positions
  # and one with the CDS positions and reads, had to do CDS in two steps
  
  FiveUTR <- filter(tmp_tidy, tmp_tidy$Pos < min(start))
  tmp_CDS <- filter(tmp_tidy, tmp_tidy$Pos >= min(start)) 
  CDS <- filter(tmp_CDS,tmp_CDS$Pos <= max(stop))
  
  # for the premade tibble for the h5 file, add the number of counts per base for 5UTRs and CDS for each .x 
  # doing counts per base allows comparison even if 5UTRs are super short/CDS is long
  
  gene_feature_reads_per_base <- tibble(
    Gene = gene, 
    fiveUTR_reads_per_base = sum(FiveUTR$Counts)/length(unique(FiveUTR$Pos)),
    CDS_reads_per_base = sum(CDS$Counts)/length(unique(CDS$Pos))
  ) %>%
    mutate(ratio = fiveUTR_reads_per_base/CDS_reads_per_base) %>%
    mutate(Sample = basename(hd_file))
  
  return(gene_feature_reads_per_base)
}

#test <- GeneFeatureTotalCountsPerBase(gene = gene_names[1], dataset = dataset, hd_file = h5_file_path[1], startpos = 1, condition_df = condition_df)

# > test
# # A tibble: 1 x 3
# Gene          fiveUTR_reads_per_base CDS_reads_per_base
# <chr>                          <dbl>              <dbl>
#   1 SPAC1002.01.1                      0              0.202


# Create function of applying GeneFeatureTotalCountsPerBase to genes,
# so it can be used in a nested map() function to apple to all datasets listed

make_sample_tables <- function(gene_names, hd_file, dataset){
  
  sample_table <- purrr::map_dfr(.x = gene_names, .f = GeneFeatureTotalCountsPerBase, hd_file = hd_file, dataset, startpos = 1)
  
  return(sample_table)
}


# filter to remove those with infinite values or less than 0.02 reads per base. 
# reduces chance that the observed change is due to random variation as a higher number is needed to be plotted

filter_tables <- function(all_genes_comparison, all_genes_one_sample, read_threshold){
  
  all_genes_comparison <- all_genes_comparison %>% filter(all_genes_comparison$fiveUTR_reads_per_base >= read_threshold & 
                                                          all_genes_comparison$CDS_reads_per_base >= read_threshold)
  
  all_genes_one_sample <- all_genes_one_sample %>% filter(all_genes_one_sample$fiveUTR_reads_per_base >= read_threshold &
                                                          all_genes_one_sample$CDS_reads_per_base >= read_threshold & 
                                                          all_genes_one_sample$Gene %in% all_genes_comparison$Gene)
  
  all_genes_comparison <- all_genes_comparison %>% filter(all_genes_comparison$Gene %in% all_genes_one_sample$Gene)
  
  ratio_change <- tibble(Gene = all_genes_one_sample$Gene,
                           control = all_genes_comparison$ratio,
                           Treated = all_genes_one_sample$ratio)
  return(ratio_change)
}


# plot the ratio of 5UTR/CDS reads per base in AT v 5UTR/CDS reads per base in no AT 
# x = y line. any points on this line show little change between conditions
# above x = y line, with.At 5UTR/CDS Ratio is greater than no.AT 5UTR/CDS Ratio, increased relative 5UTR use when starved
# below x = y line, with.At 5UTR/CDS Ratio is less than no.AT 5UTR/CDS Ratio, decreased relative 5UTR use when starved 

plot_scatter_comparison <- function(ratio_change, highlighted_genes, treatment){
  
  scatter_plot <- ggplot(ratio_change,aes(x = control, y = Treated))+
    theme_bw()+
    geom_point(color = 'sky blue')+
    geom_abline(slope = 1, intercept = 0)+
    geom_point(data = highlighted_genes, colour = 'black', size = 2)+ 
    geom_point(data = highlighted_genes, aes(colour = Gene))+
    scale_y_log10(limits = c(0.001, 100),
                  breaks = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0), 
                  labels = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0))+
    scale_x_log10(limits = c(0.001, 100),
                  breaks = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0), 
                  labels = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0))+
    coord_equal(ratio = 1)+
    labs(title = paste0(treatment, 'change in
5UTR rpb:CDS rpb value'),
         x = '5UTR rpb/CDS rpb - control',
         y = '5UTR rpb/CDS rpb - treated')+
    theme(text=element_text(size=14),
          axis.title=element_text(size=14, face='bold'),
          title = element_text(size = 14, face='bold'))
  
  return(scatter_plot)
}

# function to save plot as PDF in the desired location

save_plot_pdf <- function(scatter_plot, output_dir){
  scatter_plot %>%
  ggsave(
    filename = file.path(output_dir,"sample_v_control_relative_UTR_use.pdf"),
    width = 6, height = 5
  )
}


## check files exist 

if (!file.exists(h5_file_path)){
  stop("H5 file specified by --input argument not found.")
}
if (!file.exists(comparison_file)){
  stop("H5 file specified by --comparison argument not found.")
}
if (!file.exists(gff_in)){
  stop("GFF file specified by --gff argument not found.")
}

# run code 


print('Getting gene names')

  GFF <- readGFFAsDf(gff_in)
  gene_names <- levels(GFF$seqnames)
  
print('Creating sample tables')

  all_genes_one_sample <- make_sample_tables(gene_names, hd_file = h5_file_path, dataset)
  all_genes_comparison <- make_sample_tables(gene_names,hd_file = comparison_file, dataset)
  
print('removing genes with few or no reads in features')

  ratio_change <- filter_tables(all_genes_comparison, all_genes_one_sample, read_threshold)
  highlighted_genes <- ratio_change %>% filter(ratio_change$Gene %in% Gene_of_interest)
  
print('Creating UTR use comparison scatter plot')

  scatter_plot <- plot_scatter_comparison(ratio_change, highlighted_genes, treatment)
  
print('Saving Plots as PDF')

  save_plot_pdf(scatter_plot, output_dir)
  
print('Done')