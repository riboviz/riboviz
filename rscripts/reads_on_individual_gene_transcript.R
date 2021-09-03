# This script takes a h5 file and extracts the positions of reads for a gene of interest. 
print('Loading required files')

suppressMessages(library(getopt, quietly=T))
suppressMessages(library(here))
suppressMessages(library(optparse))

source(here::here("rscripts", "read_count_functions.R"))
source(here::here("rscripts", "stats_figs_block_functions.R"))



option_list <- list(make_option(c("-i", "--input"),
                                type = "character",
                                help = "Path input to h5 file"),
                    make_option(c("-d", "--dataset"),
                                type = "character",
                                help = "Name of the dataset being studied"),
                    make_option(c("-g", "--gff"),
                                type = "character",
                                help = "Path t the GFF3 file of the organism being studied"),
                    make_option(c("-o", "--output"),
                                type = "character",
                                help = "Path to output directory"),
                    make_option("--gene",
                                type = "character",
                                help = "Gene of interest, name used should match the name listed in H5 file so check format")
)


opt <- optparse::parse_args(OptionParser(option_list = option_list))

file_url <- opt$input
gff_in <- opt$gff
Gene_of_interest <- opt$gene
dataset <- opt$dataset
output_dir <- opt$output


## Actual code

GFF <- readGFFAsDf(gff_in)
gene_names <- levels(GFF$seqnames)
read_range <- c(10:50)


print('reads_on_individual_gene_transcript.R is running')


# Create the gene data matrix 
gene_data_matrix <- GetGeneDatamatrix(gene= Gene_of_interest,
                                      dataset = dataset,
                                      hd_file = file_url)


# extract the buffer so can identify stop codon later
buffer_left <- h5readAttributes(file_url, base::paste('/',Gene_of_interest,'/',dataset,'/reads', sep = ''))[['buffer_left']]
buffer_right <- h5readAttributes(file_url, base::paste('/',Gene_of_interest,'/',dataset,'/reads', sep = ''))[['buffer_right']]
nnt_buffer <- buffer_left

# Create a Tidy data matrix, using the the gene data matrix. set start as -nnt_buffer + 1 so the actual start codon lies on 1
 
gene_of_interest_tidy_matrix <- TidyDatamatrix(data_mat = gene_data_matrix, startpos = -nnt_buffer +1 )


# make tidy data matrix, with all the counts of different read lengths at different positions
# aim: make a tidyDataMatrix with two columns: position and count
# input: Fil1TidyDataMatrix 
# process: for position i, take the sum of the counts at that position across read lengths
#          and return a new row for the new matrix with the position and the counts 
# for each position, take all the counts of reads of different sizes, and combine them into one 
# 'total read' column for each position in the gene
# output: TidyDataMatrix

# create an empty tibble with two columns; Pos and Counts.

gene_Total_reads_at_position <- tibble(Pos =integer(),
                                            Counts = integer())

# Loop through tidy data matrix, and for each position add up all of the counts for different 
# read lengths, storing them in the new tidy data matrix
# To do: edit to use apply function

for(i in gene_of_interest_tidy_matrix[1,]$Pos:max(gene_of_interest_tidy_matrix$Pos)){
  tmp_row <- gene_of_interest_tidy_matrix %>% 
    filter(gene_of_interest_tidy_matrix$Pos ==i)
  new_row <- tibble(Pos = i, Counts = sum(tmp_row$Counts))
  gene_Total_reads_at_position <- gene_Total_reads_at_position %>% 
    bind_rows(new_row)
}


print('Creating reads on transcript plot')


# plot the data, so positions are along the x axis and number of counts is along the y axis
reads_on_transcript_plot <- ggplot(gene_Total_reads_at_position,aes(x = Pos, y = Counts))+
   theme_bw()+
   geom_col( width = 1, color = 'red')+
   scale_y_continuous(limits = c(0, max(gene_Total_reads_at_position)))+
   labs(title = paste(dataset, strsplit(basename(file_url), '.h5'), ' - ', Gene_of_interest), 
        x = 'Position relative to start codon', y = 'Number of reads') +
   theme(axis.text=element_text(size=14),
         axis.title=element_text(size=14, face='bold'),
         title = element_text(size = 16, face='bold'))

# save reads_on_transcript_plot as pdf

reads_on_transcript_plot %>%
  ggsave(
    filename = file.path(output_dir,"reads_on_transcript_plot.pdf"),
    width = 6, height = 5
  )

## reads per million plot 
# calculate total number of reads by getting sum of total reads per length 

# get the number of reads per length
read_length_data <- CalculateReadLengths(gene_names, dataset, hd_file = file_url)
# get total number of reads for dataset
total_reads <- sum(read_length_data$Counts)

#calculate reads per million for gene
reads_per_million <- gene_Total_reads_at_position %>% mutate(Counts = (Counts/total_reads)*1e6)


print('Creating reads per million plot')

reads_per_million_plot <- ggplot(reads_per_million,aes(x = Pos, y = Counts))+
  theme_bw()+
  geom_col( width = 1, color = 'red')+
  scale_y_continuous(limits = c(0, max(gene_Total_reads_at_position)))+
  labs(title = paste(dataset, strsplit(basename(file_url), '.h5'), ' - ', Gene_of_interest), 
       x = 'Position relative to start codon', y = 'Reads per million 
  reads (RPFs)') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14, face='bold'),
        title = element_text(size = 14, face='bold'))

# save as reads_per_million_plot as PDF

reads_per_million_plot %>%
  ggsave(
    filename = file.path(output_dir,"reads_per_million_plot.pdf"),
    width = 6, height = 5
  )


print('Done')