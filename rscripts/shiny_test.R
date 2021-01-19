library(shiny)
library(yaml)
library(tidyverse)
library(scales)

# define the universal plots theme
plot_theme <- theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 14))

# this should be some way to get the yaml location, it's the only thing that needs to change
# unless we want to get the sample names themselves from nextflow too
# probably it will be better to get necessary files straight from nextflow because the ribogrid,
# among other plots, will require that. 
yaml_loc <- "../vignette/simdata_multiplex_config.yaml"

# the find sample names function
find_sample_names <- function(path_to_yaml){
  # read in the yaml
  yaml <- read_yaml(path_to_yaml)
  
  # if fq_files is NULL
  if (is.null(yaml$fq_files)){
    # then it is multiplexed, get the names from the barcodes file
    # this is the location of the sample sheet
    sample_sheet_loc <- paste0("../", yaml$dir_in, yaml$sample_sheet)
    
    # get the sample names from it
    sample_names <- read_tsv(sample_sheet_loc)$SampleID
  } else {
    # it is not multiplexed, read the names from the yaml
    sample_names <- names(yaml$fq_files)
  }
  
  # construct the paths to the dirs
  sample_dir_paths <- paste0("../", yaml$dir_out, "/", sample_names)
  
  names(sample_dir_paths) <- sample_names
  
  # in the case where reads are not demultiplexed, as in Tag3
  # remove that from the vector
  do_dirs_exist <- unlist(lapply(sample_dir_paths, dir.exists))
  
  # retain only sample names for dirs that exist
  these_dirs_exist <- names(do_dirs_exist[do_dirs_exist == TRUE])
  
  sample_dir_paths_that_exist <- sample_dir_paths[names(sample_dir_paths) %in% these_dirs_exist]
  
  # load read_counts.tsv into it's own global variable
  read_counts_df <<- read_tsv(paste0("../", yaml$dir_out, "/read_counts.tsv"), skip = 5)
  
  return(sample_dir_paths_that_exist)
}

sample_names <- find_sample_names(yaml_loc)

# get read length dists into one data frame
read_length_df <- lapply(sample_names, function(x){
  # path to the read length file
  file_loc <- paste0(x, "/read_lengths.tsv")
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>% 
  bind_rows(.id = "samplez")

# Define UI
ui <- pageWithSidebar(
  
  # App title
  headerPanel("Riboviz"),
  
  # Sidebar panel for inputs
  sidebarPanel(
    checkboxGroupInput(inputId = "sample",
                       label = "Sample:",
                       choiceNames = names(sample_names),
                       choiceValues = names(sample_names),
                       selected = names(sample_names)[1]),
    width = 2
  ),
  
  # Main panel for displaying outputs
  mainPanel(
    plotOutput("read_length_dist_plot")
  )
)

# Define server logic to plot various things
server <- function(input, output, session) {
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
}

shinyApp(ui, server)