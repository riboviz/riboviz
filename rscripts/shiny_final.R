#!/usr/bin/env Rscript

# Uncomment the following 6 lines to take an argument from command line
# args = commandArgs(trailingOnly=TRUE)
# 
# # Check for a single file argument from the command line
# if (length(args) == 0) {
#   stop("Provide the path to the YAML file used in your riboviz2 run", call.=FALSE)
# }

# load required packages
library(shiny)
library(yaml)
library(tidyverse)
library(scales)

### Get information from the YAML ###
# read in the yaml, the path here would get replaced with "args" if you're grabbing it 
# from command line
yaml <- read_yaml("/data2/john/projects/riboviz/riboviz/vignette/simdata_multiplex_config.yaml")

# use this function to pull the sample names and directories from the yaml
# the sample names and dirs will be used to find data files.
find_sample_names <- function(yaml_file){
 # if there are no entries in fq_files, data is multiplexed
  if (is.null(yaml_file$fq_files)){
    # get the names from the barcodes file, this is the location of the sample sheet
    sample_sheet_loc <- paste0("../", yaml_file$dir_in, yaml_file$sample_sheet)
    
    # get the sample names from it
    sample_names <- read_tsv(sample_sheet_loc)$SampleID
  } else {
    # it is not multiplexed, read the names from the yaml_file
    sample_names <- names(yaml_file$fq_files)
  }
  
  # construct the paths to the dirs
  sample_dir_paths <- paste0("../", yaml_file$dir_out, "/", sample_names)
  
  names(sample_dir_paths) <- sample_names
  
  # in the case where reads are not demultiplexed, as in Tag3 in the vignette
  # remove that from the vector
  do_dirs_exist <- unlist(lapply(sample_dir_paths, dir.exists))
  
  # retain only sample names for dirs that exist
  these_dirs_exist <- names(do_dirs_exist[do_dirs_exist == TRUE])
  
  sample_dir_paths_that_exist <- sample_dir_paths[names(sample_dir_paths) %in% these_dirs_exist]
  
  return(sample_dir_paths_that_exist)
}

# use the function
sample_names <- find_sample_names(yaml)

###############################################################################
# Data acquisition
###############################################################################

### Read counts bar plot
# rename the categories for the graph
graph_cats <- data.frame(
  Description = c("Demultiplexed reads",
                  "Reads after removal of sequencing library adapters",
                  "Reads with rRNA and other contaminating reads removed by alignment to rRNA index files",
                  "Reads after trimming of 5' mismatches and removal of those with more than 2 mismatches",
                  "Reads aligned to ORFs index files",
                  "Deduplicated reads"),
  short_name = c("Reads after demultiplexing",
                 "Reads after adapter removal",
                 "Reads after contaminant removal",
                 "Reads after trimming of 5' mismatches",
                 "Reads aligned to ORFs",
                 "Reads after deduplication"))

# Load the data frame and add the short names
read_counts_df <- read_tsv(paste0("../", yaml$dir_out, "/read_counts.tsv"), skip = 5) %>% 
  left_join(graph_cats, by = "Description") %>% 
  filter(Description %in% graph_cats$Description & !is.na(SampleName) & SampleName != "Unassigned")

# find which categories actually occur
occurring_cats <- unique(read_counts_df$short_name)

read_counts_df <- read_counts_df %>% 
  mutate(short_name = factor(short_name, levels = graph_cats$short_name[graph_cats$short_name %in% occurring_cats]))

### Read length distributions
read_length_df <- lapply(sample_names, function(x){
  # path to the read length file
  file_loc <- paste0(x, "/read_lengths.tsv")
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>% 
  bind_rows(.id = "samplez")

### Periodicity line plot
periodicity_df <- lapply(sample_names, function(x){
  # path to the read length file
  file_loc <- paste0(x, "/3nt_periodicity.tsv")
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>% 
  bind_rows(.id = "samplez")

### Read counts by frame bar plot
frame_bar_df <- lapply(sample_names, function(x){
  # path to the read length file
  file_loc <- paste0(x, "/3ntframe_bygene.tsv")
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>% 
  bind_rows(.id = "samplez") %>%
  group_by(samplez) %>%
  summarise(Frame_0 = sum(Ct_fr0),
            Frame_1 = sum(Ct_fr1),
            Frame_2 = sum(Ct_fr2)) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "Frame", values_to = "Count")

### Read counts by frame boxplot
frame_box_df <- lapply(sample_names, function(x){
  # path to the read length file
  file_loc <- paste0(x, "/3ntframe_bygene.tsv")
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>% 
  bind_rows(.id = "samplez") %>%
  mutate(Frame_0 = Ct_fr0/(Ct_fr0 + Ct_fr1 + Ct_fr2),
         Frame_1 = Ct_fr1/(Ct_fr0 + Ct_fr1 + Ct_fr2),
         Frame_2 = Ct_fr2/(Ct_fr0 + Ct_fr1 + Ct_fr2)) %>%
  select(samplez, gene, Frame_0, Frame_1, Frame_2) %>%
  pivot_longer(names_to = "Frame", cols = 3:5)

### Position specific normalized reads 
pos_sp_df <- lapply(sample_names, function(x){
  # path to the read length file
  file_loc <- paste0(x, "/pos_sp_rpf_norm_reads.tsv")
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>% 
  bind_rows(.id = "samplez")

### Ribogrid 
ribogrid_df <- lapply(sample_names, function(x){
  # path to the read length file
  file_loc <- paste0(x, "/gene_position_length_counts_5start.tsv")
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>% 
  bind_rows(.id = "samplez")

### Sample correlations heatmap
collated_tpms_df <- read_tsv(paste0("../", yaml$dir_out, "/TPMs_collated.tsv"), skip = 4)

###############################################################################
# Plotting
###############################################################################

### define the universal plot theme
plot_theme <- theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 17))

### Start the server
server <- function(input, output, session) {
  ### Read counts summary
  # Read count summary at the various steps as a bar plot
  # Requires read_counts.tsv
  output$read_counts_bar_plot <- renderPlot({
    read_counts_df %>% 
      filter(SampleName %in% input$sample) %>% 
      ggplot(., aes(SampleName, NumReads, fill = short_name))+
      geom_col(position = position_dodge())+
      plot_theme+
      scale_fill_discrete(name = NULL)+
      labs(x = NULL,
           y = "Read count",
           title = "Read counts per step")+
      guides(fill = guide_legend(ncol = 1))
  })
  ### Sample to sample cors
  # load the TPM df
  cor.df <- read_tsv("/data2/john/projects/riboviz/riboviz/vignette/simdata_multiplex_output/TPMs_collated.tsv", skip = 4) %>% 
    column_to_rownames("ORF")
  
  # sort the columns to ensure proper display of heatmap
  cor.cols <- sort(names(cor.df))
  
  # create correlations
  cor.df <- cor.df %>% 
    select(all_of(cor.cols)) %>% 
    as.matrix() %>%
    cor(use = "pairwise.complete.obs")
  
  # remove lower triangle and autocorrelations
  cor.df[lower.tri(cor.df)] <- NA
  
  cor.df[cor.df == 1] <- NA
  
  # convert to a df and reshape
  cor.df2 <- as_tibble(cor.df, rownames = "samp1") %>% 
    pivot_longer(cols = where(is.numeric), names_to = "samp2", values_to = "R")
    
  # plot it
  output$sample_cors_plot <- renderPlot({
    ggplot(cor.df2, aes(samp1, samp2, fill = R, label = signif(R,2)))+
      geom_raster()+
      scale_fill_viridis_c(option = "B", na.value = NA)+
      theme(panel.background = element_blank(),
            text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.grid.major.x = element_line(color = "grey50", linetype = 3),
            axis.ticks = element_blank())+
      labs(x = NULL, y = NULL,
           subtitle = "Pairwise correlations based on TPM")+
      geom_text()
  })
  
  # QUALITY CONTROL VISUALIZATIONS
  
  # Read length distribution plot
  # Requires read_lengths.tsv
  output$read_length_dist_plot <- renderPlot({
    read_length_df %>%
      filter(samplez %in% input$sample) %>%
      ggplot(., aes(Length, Counts, fill = samplez))+
      geom_col(width = 1)+
      facet_wrap(~samplez, scales = "free_y")+
      scale_fill_discrete(name = NULL, guide = FALSE)+
      plot_theme+
      labs(title = "Read length distributions",
           x = "Read length")+
      scale_x_continuous(breaks = seq(10, 50, 5))+
      scale_y_continuous(breaks = breaks_pretty(n = 4))
  })
  
  # Periodicity metagene plot for both 5' and 3' end as a line plot
  # Requires 3nt_periodicity.tsv
  output$periodicity_line_plot <- renderPlot({
    periodicity_df %>%
      filter(samplez %in% input$sample) %>%
      ggplot(., aes(x = Pos, y = Counts, color = samplez)) + 
      geom_line() +
      scale_color_discrete(name = NULL, guide = FALSE) +
      facet_grid(cols = vars(samplez), rows = vars(End), scales = "free_y") +
      labs(title = "3-nucleotide periodicity")+
      plot_theme
  })
  
  # Read counts per frame as a bar plot
  # Requires 3ntframe_bygene.tsv
  output$frame_counts_bar_plot <- renderPlot({
    frame_bar_df %>% 
      filter(samplez %in% input$sample) %>% 
      ggplot(., aes(x = Frame, y = Count, fill=samplez)) + 
      geom_bar(stat="identity", pos = "dodge") + 
      scale_fill_discrete(name = NULL)+
      labs(title = "Read counts per frame")+
      plot_theme+
      scale_color_discrete(guide = FALSE)
  })
  
  # Proportion of reads per frame for each gene as a boxplot
  # Requires 3ntframe_bygene.tsv
  output$frame_proportions_boxplot <- renderPlot({
    frame_box_df %>%
      filter(samplez %in% input$sample) %>% 
      ggplot(., aes(x=Frame, y=value, fill=samplez)) +
      labs(y="Proportion", title = "Proportion of reads per frame per gene")+
      scale_fill_discrete(name = NULL)+
      geom_boxplot()+
      plot_theme
  })

  # Position specific normalized reads, barplot
  # Requires pos_sp_rpf_norm_reads.tsv
  output$pos_sp_plot3 <- renderPlot({
    pos_sp_df %>%
      filter(End == "3'") %>% 
      ggplot(., aes(Position, Mean, color = samplez))+
      geom_line()+
      facet_grid(End ~ samplez, scales = "free")+
      plot_theme+
      scale_color_discrete(guide = FALSE)
  })
  
  output$pos_sp_plot5 <- renderPlot({
    pos_sp_df %>%
      filter(End == "5'") %>% 
      ggplot(., aes(Position, Mean, color = samplez))+
      geom_line()+
      facet_grid(End ~ samplez, scales = "free")+
      plot_theme+
      scale_color_discrete(guide = FALSE)
  })
  
  # ANALYSIS VISUALIZATIONS
  
  # Abundance of specific gene(s) compared to all, density plot
  # Requires TPMs_collated.tsv
  output$tpm_density_plot <- renderPlot({
    collated_tpms_df %>%
      pivot_longer(where(is.numeric), names_to = "samplez", values_to = "tpm") %>%
      filter(samplez %in% input$sample) %>%
      mutate(labz = ifelse(ORF %in% input$gene, ORF, NA),
             xint = ifelse(ORF %in% input$gene, tpm, NA)) %>%
      ggplot(., aes(x = tpm, fill = samplez)) +
      geom_density(alpha = 0.5) +
      geom_vline(aes(xintercept = xint), size = 1) +
      scale_x_log10()+
      facet_wrap(~samplez)+
      scale_fill_discrete(guide = FALSE)+
      plot_theme+
      labs(title="TPM density", x="TPM", y="Density")
  })
  
  # Ribogrid plot
  # Need to find plot_ribogrid() in Riboviz code
  output$ribogrid_plot <- renderPlot({
    ribogrid_df %>%
      ggplot(., aes(Pos, ReadLen, fill=Counts))+
               geom_raster()+
               scale_fill_gradient(low="white", high="navy")+
               plot_theme+
               labs(title="Ribogrid", x="Read position", y="Read length")
  })
  
  # Ribogrid bar plot
  output$ribogridbar_plot <- renderPlot({
    ribogrid_df %>%
      filter(ReadLen >= min(input$ribogrid_len_range) & ReadLen <= max(input$ribogrid_len_range)) %>%
      filter(Pos >= min(input$ribogrid_pos_range) & Pos <= max(input$ribogrid_pos_range)) %>%
      ggplot(., aes(Pos, ReadLen, Counts))+
      geom_col()+
      facet_wrap(~ReadLen, ncol=1, strip.position = "right")+
      plot_theme+
      labs(x="Read position", y="Read length", title="Ribogrid (bar)")
  })
  
  
  # TPM vs other features, scatterplots
  # Requires tpms.tsv and other features
  # TODO
  
}


#
# END SERVER FUNCTION
#

#
# START UI FUNCTION
#

ui <- fluidPage(
  
  # App title
  headerPanel("Riboviz: selective data visualization"),
  
  # Sidebar panel for inputs
  fluidRow(
    column(12,
           wellPanel(
             checkboxGroupInput(inputId = "sample",
                                label = "Sample:",
                                choiceNames = names(sample_names),
                                choiceValues = names(sample_names),
                                selected = names(sample_names)[1], 
                                inline = TRUE,
                                width = "100%"),
             conditionalPanel(condition="input.panelId == 'Analysis'", 
                              textInput(inputId = "gene",
                                        label = "Gene:",
                                        width = "100%"),
                              sliderInput("ribogrid_len_range", "Ribogrid read length range:",
                                          min = 0, max = 50,
                                          value = c(20,25)),
                              sliderInput("ribogrid_pos_range", "Ribogrid read position range:",
                                          min = -24, max = 50,
                                          value = c(-10,10)),
             )
             
           )
    )
  ),
  
  # Main panel for displaying outputs
  mainPanel(
    width = 12,
    fluidRow(
      tabsetPanel(
        tabPanel("Read count summary",
                 headerPanel(""),
                 plotOutput("read_counts_bar_plot"),
                 headerPanel(""),
                 plotOutput("sample_cors_plot")
        ),
        
        tabPanel("Read length distributions",
                 headerPanel(""),
                 plotOutput("read_length_dist_plot")
        ),
        
        tabPanel("3nt periodicity",
                 headerPanel(""),
                 plotOutput("periodicity_line_plot"),
                 headerPanel(""),
                 plotOutput("frame_counts_bar_plot"),
                 headerPanel(""),
                 plotOutput("frame_proportions_boxplot")
        ),
        
        tabPanel("3nt periodicity",
                 headerPanel(""),
                 plotOutput("periodicity_line_plot"),
                 headerPanel(""),
                 plotOutput("frame_counts_bar_plot"),
                 headerPanel(""),
                 plotOutput("frame_proportions_boxplot")
        ),
        
        tabPanel("Ribogrids",
                 headerPanel(""),
                 plotOutput("ribogrid_plot"),
                 headerPanel(""),
                 plotOutput("ribogridbar_plot", height = "1000px")
        
        tabPanel("Position Specific read densities",
                 headerPanel(""),
                 plotOutput("pos_sp_plot3"),
                 headerPanel(""),
                 plotOutput("pos_sp_plot5")
        ),
        
        tabPanel("Gene specific plots",
                 headerPanel(""),
                 plotOutput("tpm_density_plot"),
                 
        ),
        id = "panelId"
      )
    )
  )
)

#
# END UI FUNCTION
#

# I cast spell, ~run shiny app~
shinyApp(ui, server)