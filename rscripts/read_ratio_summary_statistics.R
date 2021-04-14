suppressMessages(library(here))

#make summary statistics. Compare how the ratio of reads in the UTRs:CDS changes
# between conditions for Fil1
# compare this to the global ratio change

# take all the reads before start codon and add them together
# take all the reads in the CDS and add them together 
# Compare. set CDS to 1 by dividing by CDS


# make a loop to do that for all the reads before start codons in the h5, or in CDS

# input: gff or fasta for names, h5 file
# create a dataframe with 6 columns; gene names, number of reads 5UTR per base, number of reads CDS per base, ratio of 5UTR to CDS in control, Treatment, Sample
# to do that need: for each condition dataframe of gene, reads CDS, read 5'UTR 
# make that; for gene in genenames, get all of the reads up to position buffer left, add to dataframe
# get all reads in CDS, add to dataframe 
# do reads per base to account for different lenghts of CDS and UTRs 

source(here::here('rscripts','read_count_functions.R'))

# get inputs, gene names, and paths to h5 files 
# expecting example-datasets to be on same level as main riboviz repo, currently using branch Pombe_annotation_full_UTR-70
gff_in <- here::here('..','example-datasets','fungi', 'schizosaccharomyces','annotation','Schizosaccharomyces_pombe_full_UTR_or_50nt_buffer.gff3')
GFF <- readGFFAsDf(gff_in)
# get the name of genes 
gene_names <- levels(GFF$seqnames)
# file is a list of paths to h5 files, files needed to be read
h5_file_path <- c(here::here('D-Sp_2018','output','wt.AT.ribo.4_s','wt.AT.ribo.4_s.h5'),
          #'wt.AT.ribo.3_s/8c5f643d83fe7ce08246d386e8303f/wt.AT.ribo.3_s.h5',
          #'wt.AT.CHX.ribo.11_s/wt.AT.CHX.ribo.11_s.h5',
          #'wt.AT.CHX.ribo.13/ec18ef21336806402ba541b8f12717/wt.AT.CHX.ribo.13_s.h5',
          #'wt.AT.noCHX.ribo.11_s/wt.AT.noCHX.ribo.11_s.h5',
          #'wt.AT.noCHX.ribo.13_s/6c237de3fe552baf04b30fdb08bbcc/wt.AT.noCHX.ribo.13_s.h5',
          here::here('D-Sp_2018','output','wt.noAT.ribo.4_s','wt.noAT.ribo.4_s.h5')
          #'wt.noAT.ribo.3_s/dc0ac21aa1b8dfab77ec9366797b4e/wt.noAT.ribo.3_s.h5',
          #'wt.noAT.CHX.ribo.11_s/wt.noAT.CHX.ribo.11_s.h5',
          #'wt.noAT.CHX.ribo.13_s/wt.noAT.CHX.ribo.13_s.h5',
          #'wt.noAT.noCHX.ribo.11_s/wt.noAT.noCHX.ribo.11_s.h5',
          #'wt.noAT.noCHX.ribo.13_s/wt.noAT.noCHX.ribo.13_s.h5'
          )
# could use yaml parameter in the future to be integrated 
dataset <- 'D-Sp_2018'


# create a tibble for each sample, marked with sample name, and get total reads per base 5UTR and CDS
# reads per base is used to allow comparison no matter the length of the 5UTR or CDS 

for(hd_file in h5_file_path){
  #for each h5 file listed, make an  
   condition_df <- tibble(Gene = character(),
          fiveUTR_reads_per_base = integer(),
          CDS_reads_per_base = integer())

  for(gene in gene_names){
    #make a matric to be processed for each gene
    tmp_mat <- GetGeneDatamatrix(gene, dataset = dataset, hd_file = hd_file)
    #get the start and stop codons
    start <- h5readAttributes(hd_file,base::paste('/',gene,'/',dataset,'/reads', sep = ''))[['start_codon_pos']]
    stop <- h5readAttributes(hd_file, base::paste('/',gene,'/',dataset,'/reads', sep = ''))[['stop_codon_pos']]
    # turn matric into a tidy datamatrix 
    tmp_tidy <- TidyDatamatrix(tmp_mat, 
                             startpos = 1 )
    # filter the tmp_tidy matric into two new matrixes, on with 5UTR positions
    # and one with the CDS positions and reads, had to do CDS in two steps 
    FiveUTR <- filter(tmp_tidy, tmp_tidy$Pos < min(start))
    tmp_CDS <- filter(tmp_tidy, tmp_tidy$Pos >= min(start)) 
    CDS <- filter(tmp_CDS,tmp_CDS$Pos <= max(stop))
    # for the premade tibble for the 5 file, add the number of counts per base for 5UTRs and CDS for each gene 
    # doing counts per base allows comparison even if 5UTRs are super short/CDS is long
    new_row <- tibble(Gene = gene, 
                    fiveUTR_reads_per_base = sum(FiveUTR$Counts)/length(unique(FiveUTR$Pos)),
                    CDS_reads_per_base = sum(CDS$Counts)/length(unique(CDS$Pos)))
    # add row to the dataframe for the condition  
    condition_df<- condition_df %>% bind_rows(new_row)
     
  }
   # save the condition_df into an object with the condition name 
   # create a column 'ratio' by dividing the 5UTR counts per base by CDS counts per base
   # create new column 'Treatment' by exracting the treatment
  # create column Sample to get the sample origin of the data 
   condition_df <- condition_df %>% mutate(ratio = fiveUTR_reads_per_base/CDS_reads_per_base,
            Treatment = paste(unlist(strsplit(basename(hd_file), split = "[.]"))[2:3], collapse = '.'),
            Sample = basename(hd_file),
            .keep = 'all')
   
   assign(paste('gene_reads_5UTR_CDS.',basename(hd_file), sep = ''),condition_df) 
     
}



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
    mutate(ratio = fiveUTR_reads_per_base/CDS_reads_per_base)
    
  return(gene_feature_reads_per_base)
}

#test <- GeneFeatureTotalCountsPerBase(gene = gene_names[1], dataset = dataset, hd_file = h5_file_path[1], startpos = 1, condition_df = condition_df)

# > test
# # A tibble: 1 x 3
# Gene          fiveUTR_reads_per_base CDS_reads_per_base
# <chr>                          <dbl>              <dbl>
#   1 SPAC1002.01.1                      0              0.202

gene_names_less <- gene_names[1:3]

All_genes_sample1 <- purrr::map_dfr(.x = gene_names_less, .f = GeneFeatureTotalCountsPerBase, dataset = dataset, hd_file = h5_file_path[1], startpos = 1) %>%
  mutate(Sample = basename(h5_file_path[1]))
All_genes_sample2 <- purrr::map_dfr(.x = gene_names_less, .f = GeneFeatureTotalCountsPerBase, dataset = dataset, hd_file = h5_file_path[2], startpos = 1) %>%
  mutate(Sample = basename(h5_file_path[2]))

All_samples_all_genes <- bind_rows(All_genes_sample1, All_genes_sample2)
#ApplyGeneFeatureTotalCountsPerBaseToSamples <- function(.x = gene_names_less, dataset = dataset, hd_file = .y, startpos = 1){
#   
#   all_genes_one_sample <- purrr::map_dfr(.x = gene_names_less, .f = GeneFeatureTotalCountsPerBase, hd_file = hd_file, dataset = dataset, startpos = 1)
#   
#   return(all_genes_one_sample)
#   
# } 
# trying
#test <- purrr::map_dfr(.x = gene_names_less, .y = h5_file_path, .f =ApplyGeneFeatureTotalCountsPerBaseToSamples, dataset = dataset, startpos = 1)


# combine all of the dataframes into one, so I can make a box plot 
# needs to be automated
All_samples_one_df <- rbind(gene_reads_5UTR_CDS.wt.AT.CHX.ribo.11_s.h5,
                            gene_reads_5UTR_CDS.wt.AT.CHX.ribo.13_s.h5,
                            gene_reads_5UTR_CDS.wt.AT.noCHX.ribo.11_s.h5,
                            gene_reads_5UTR_CDS.wt.AT.noCHX.ribo.13_s.h5,
                            gene_reads_5UTR_CDS.wt.AT.ribo.3_s.h5, 
                            gene_reads_5UTR_CDS.wt.AT.ribo.4_s.h5, 
                            gene_reads_5UTR_CDS.wt.noAT.CHX.ribo.11_s.h5,
                            gene_reads_5UTR_CDS.wt.noAT.CHX.ribo.13_s.h5,
                            gene_reads_5UTR_CDS.wt.noAT.noCHX.ribo.11_s.h5,
                            gene_reads_5UTR_CDS.wt.noAT.noCHX.ribo.13_s.h5,
                            gene_reads_5UTR_CDS.wt.noAT.ribo.3_s.h5,
                            gene_reads_5UTR_CDS.wt.noAT.ribo.4_s.h5)

# get gene of interest information, so can highlight Fil1, or any other 
genes_of_interest <- All_samples_one_df %>% filter(All_samples_one_df$Gene %in% c('SPCC1393.08.1',))

# create box plots comparing overall ratios of no AT and AT treated sample 
# create a box plot 
# the higher the position on the boxplot, the larger the 5UTR portion of the ration, thus the more reads in the 5UTR
# If the number of reads in the 5UTR increases between conditions, then use of 5UTR (and degredation) increases compared to that in the CDS
# if 5UTR use compared to CDS increases between conditions then the point for that gene will be higher on the graph. 
# boxplot seems to automatically remove any values where one of the values is 0, but may need to do this earlier. 
ggplot(All_samples_one_df, aes(x = Treatment, y = ratio))+
  geom_boxplot()+
  scale_y_log10()+
  geom_point(data = genes_of_interest,aes(col = Gene) )+
  theme(legend.position = 'right')+
  labs(title = 'difference in 5UTR to CDS reads per base ratio between treatments',
       y = '5UTR_reads_per_base/CDS_reads_per_base')
       


# make a jitter plot of ratio change 
ratio_change <- tibble(Gene = gene_reads_5UTR_CDS.wt.AT.ribo.4_s.h5$Gene,
                       no_AT_Ratio = gene_reads_5UTR_CDS.wt.noAT.ribo.4_s.h5$ratio,
                       AT_Ratio = gene_reads_5UTR_CDS.wt.AT.ribo.4_s.h5$ratio)

highlighted_genes <- ratio_change %>% filter(ratio_change$Gene %in% c('SPCC1393.08.1'))
# plot the ratio of 5UTR/CDS reads per base in AT v 5UTR/CDS reads per base in no AT 
# AT presence == starvation 
# x = y line. any points on this line are unchanged
# above x = y line, with.At 5UTR/CDS Ratio is greater than no.AT 5UTR/CDS Ratio, increased relative 5UTR use when starved
# below x = y line, with.At 5UTR/CDS Ratio is less than no.AT 5UTR/CDS Ratio, decreased relative 5UTR use when starved 
ggplot(ratio_change,aes(x = no_AT_Ratio, y = AT_Ratio))+
  geom_jitter(color = 'red')+
  geom_point(data=highlighted_genes, colour="green") +
  geom_text(data=highlighted_genes, label= highlighted_genes$Gene, vjust=1.5, size = 2.5)+
  geom_abline(slope = 1, intercept = 0)+
  scale_y_log10()+
  scale_x_log10()+
  labs(title = 'Change in 5UTR to CDS reads per base ratio',
       x = '5UTR_reads_per_base/CDS_reads_per_base - no.AT',
       y = '5UTR_reads_per_base/CDS_reads_per_base - with.AT')
