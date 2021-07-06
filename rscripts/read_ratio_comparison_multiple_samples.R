

#make summary statistics. Compare how the ratio of reads in the UTRs:CDS changes
# between conditions for gene of interest 
# compare this to the global ratio change

# take all the reads before start codon and add them together
# take all the reads in the CDS and add them together 
# Compare. set CDS to 1 by dividing by CDS


# input: gff or fasta for names, h5 file
# create a dataframe with 6 columns; gene names, number of reads 5UTR per base, number of reads CDS per base, ratio of 5UTR to CDS in control, Treatment, Sample
# to do that need: for each condition dataframe of gene, reads CDS, read 5'UTR 
# make that; for gene in genenames, get all of the reads up to position buffer left, add to dataframe
# get all reads in CDS, add to dataframe 
# do reads per base to account for different lenghts of CDS and UTRs 

# have multiple arguments for multiple files being added ie have one argumet for a file and another for the control 
print('Process starting')

suppressMessages(library(here))
suppressMessages(library(argparse))
source(here::here('rscripts','read_count_functions.R'))

parser <- ArgumentParser()
parser$add_argument('-i', '--input',
                    help="input H5 file, Path to the H5 file containing the data to be studied. Able to take multiple files, with file paths separated by a space",
                    type="character", nargs='+')
parser$add_argument('-g', '--gff', help="Path to GFF3 file corresponding to the species being studies", type="character")
parser$add_argument('--gene', help="Gene of interest, name used should match the name listed in H5 file so check format", type="character")
parser$add_argument('-d', '--dataset', help='Name of the dataset being run, ie D-Sp_2018, should match the dataset listed in H5 file', type="character")
parser$add_argument('-o', '--output', help='Output directory for plots',default=".")

args <- parser$parse_args()

h5_file_path <- args$input
gff_in <- args$gff
Gene_of_interest <- args$gene
dataset <- args$dataset
output_dir <- args$output

# GeneFeatureTotalCountsPerBase takes a gene, splits it into features (5UTR and CDS),
# and returns the reads per base in the features of that gene in the form of a tibble
GeneFeatureTotalCountsPerBase <- function(gene = gene, dataset = dataset, hd_file = h5_file_path, startpos = 1){
  
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

ApplyGeneFeatureTotalCountsPerBaseToSamples <- function(gene = gene_names, dataset = dataset, hd_file = .x, startpos = 1){
  
  all_genes_one_sample <- purrr::map_dfr(.x = gene_names, .f = GeneFeatureTotalCountsPerBase, hd_file = hd_file, dataset = dataset, startpos = 1)
   
  return(all_genes_one_sample)
  
} 


#trying
# test <- purrr::map_dfr(.x = h5_file_path, .f =ApplyGeneFeatureTotalCountsPerBaseToSamples, gene = gene_names_less, dataset = dataset, startpos = 1)
# test
# # A tibble: 6 x 4
# Gene           fiveUTR_reads_per_base CDS_reads_per_base  ratio
# <chr>                           <dbl>              <dbl>  <dbl>
#   1 SPAC1002.01.1                   0                  0.202 0     
# 2 SPAC1002.02.1                   0.212              1.13  0.187 
# 3 SPAC1002.03c.1                  0.140              1.76  0.0794
# 4 SPAC1002.01.1                   0                  0.131 0     
# 5 SPAC1002.02.1                   0.121              0.599 0.203 
# 6 SPAC1002.03c.1                  0.115              0.967 0.119 


# filter to remove those with infinite values or less than 0.02 reads per base. 
# reduces chance that the observed change is due to random variation as a higher number is needed to be plotted

filter_tables <- function(all_samples_all_genes){
  all_samples_all_genes <- all_samples_all_genes %>% filter(all_samples_all_genes$fiveUTR_reads_per_base >= 0.02 &
                                                              all_samples_all_genes$CDS_reads_per_base >= 0.02)
  return(all_samples_all_genes)
}



# create box plots comparing the spread of 5UTR:CDS ratios between samples, highlighting desired genes 
# create a box plot 
# the higher the position on the boxplot, the larger the 5UTR portion of the ration, thus the more reads in the 5UTR
# If the number of reads in the 5UTR increases between conditions, then use of 5UTR (and degredation) increases compared to that in the CDS
# if 5UTR use compared to CDS increases between conditions then the point for that gene will be higher on the graph. 
# boxplot seems to automatically remove any values where one of the values is 0, but may need to do this earlier. 
  
plot_boxplot <- function(all_samples_all_genes, highlighted){
  
  box_plot <- ggplot(all_samples_all_genes, aes(x = Sample, y = ratio))+
    geom_jitter(col='sky blue', alpha = 0.1)+
    geom_boxplot()+
    theme_bw()+
    geom_point(data = highlighted, size = 2.7, colour = 'black')+
    geom_point(data= highlighted, aes(colour = Gene), size = 2)+
    scale_y_log10(limits = c(0.001, 100),
                  breaks = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0), 
                  labels = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0))+
    theme(text=element_text(size=14),
          axis.title=element_text(size=14, face='bold'),
          title = element_text(size = 16, face='bold'))+
    labs(title = '5UTR:CDS usage in different samples',
         y = '5UTR rpb/CDS rpb', size = 2)
  return(box_plot)
}


# function to save plot as PDF in the desired location

save_plot_pdf <- function(box_plot, output_dir){
  box_plot %>%
    ggsave(
      filename = file.path(output_dir,"UTR_use_in_different_samples.pdf"),
      width = 6, height = 5
    )
}

## actual code

print('Getting gene names')
  GFF <- readGFFAsDf(gff_in)
  gene_names <- levels(GFF$seqnames)
print('Creating table of samples')
  all_samples_all_genes <- purrr::map_dfr(.x = h5_file_path,
                                          .f =ApplyGeneFeatureTotalCountsPerBaseToSamples, 
                                          gene = gene_names,
                                          dataset = dataset,
                                          startpos = 1)
print('Removing genes with few or no reads')
  all_samples_all_genes <- filter_tables(all_samples_all_genes)
  highlighted <- all_samples_all_genes %>% 
    filter(all_samples_all_genes$Gene %in% Gene_of_interest)
print('Plotting boxplot')
  box_plot <- plot_boxplot(all_samples_all_genes, highlighted)
print('Saving boxplot as PDF')
  save_plot_pdf(box_plot, output_dir)
print('Done')
