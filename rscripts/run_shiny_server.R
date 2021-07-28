### Get information from the YAML ###
library(yaml)

# load required packages
library(shiny)
library(tidyverse)
library(scales)
library(here)
library(rhdf5)

args <- commandArgs(TRUE)
if (!is.na(args[1])) {
  yaml <- read_yaml(args[1])
} else {
  return("Provide a path to a yaml file used in a riboviz run.")
}

# use this function to pull the sample names and directories from the yaml
# the sample names and dirs will be used to find data files.
find_sample_names <- function(yaml_file) {
  # if there are no entries in fq_files, data is multiplexed
  if (is.null(yaml_file$fq_files)) {
    # get the names from the barcodes file, this is the location of the sample sheet
    sample_sheet_loc <-
      normalizePath(here(yaml_file$dir_in, yaml_file$sample_sheet))
    
    # get the sample names from it
    sample_names <- read_tsv(sample_sheet_loc)$SampleID
  } else {
    # it is not multiplexed, read the names from the yaml_file
    sample_names <- names(yaml_file$fq_files)
  }
  
  # construct the paths to the dirs
  sample_dir_paths <-
    normalizePath(here(yaml_file$dir_out, sample_names))
  
  names(sample_dir_paths) <- sample_names
  
  # in the case where reads are not demultiplexed, as in Tag3 in the vignette
  # remove that from the vector
  do_dirs_exist <- unlist(lapply(sample_dir_paths, dir.exists))
  
  # retain only sample names for dirs that exist
  these_dirs_exist <- names(do_dirs_exist[do_dirs_exist == TRUE])
  
  sample_dir_paths_that_exist <-
    sample_dir_paths[names(sample_dir_paths) %in% these_dirs_exist]
  
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
  Description = c(
    "Demultiplexed reads",
    "Reads after removal of sequencing library adapters",
    "Reads with rRNA and other contaminating reads removed by alignment to rRNA index files",
    "Reads after trimming of 5' mismatches and removal of those with more than 2 mismatches",
    "Reads aligned to ORFs index files",
    "Deduplicated reads"
  ),
  short_name = c(
    "Reads after demultiplexing",
    "Reads after adapter removal",
    "Reads after contaminant removal",
    "Reads after trimming of 5' mismatches",
    "Reads aligned to ORFs",
    "Reads after deduplication"
  )
)

# Load the data frame and add the short names
read_counts_df <-
  read_tsv(normalizePath(here(yaml$dir_out, "read_counts.tsv")), skip = 5) %>%
  left_join(graph_cats, by = "Description") %>%
  filter(
    Description %in% graph_cats$Description &
      !is.na(SampleName) & SampleName != "Unassigned"
  )

# find which categories actually occur
occurring_cats <- unique(read_counts_df$short_name)

read_counts_df <- read_counts_df %>%
  mutate(short_name = factor(short_name, levels = graph_cats$short_name[graph_cats$short_name %in% occurring_cats]))

### Read length distributions
read_length_df <- lapply(sample_names, function(x) {
  # path to the read length file
  file_loc <- normalizePath(here(x, "read_lengths.tsv"))
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>%
  bind_rows(.id = "samplez")

### Periodicity line plot
periodicity_df <- lapply(sample_names, function(x) {
  # path to the read length file
  file_loc <- normalizePath(here(x, "3nt_periodicity.tsv"))
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>%
  bind_rows(.id = "samplez")

### Read counts by frame bar plot
frame_bar_df <- lapply(sample_names, function(x) {
  # path to the read length file
  file_loc <- normalizePath(here(x, "3ntframe_bygene.tsv"))
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>%
  bind_rows(.id = "samplez") %>%
  group_by(samplez) %>%
  summarise(
    Frame_0 = sum(Ct_fr0),
    Frame_1 = sum(Ct_fr1),
    Frame_2 = sum(Ct_fr2)
  ) %>%
  pivot_longer(cols = where(is.numeric),
               names_to = "Frame",
               values_to = "Count")

### Read counts by frame boxplot
frame_box_df <- lapply(sample_names, function(x) {
  # path to the read length file
  file_loc <- normalizePath(here(x, "3ntframe_bygene.tsv"))
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>%
  bind_rows(.id = "samplez") %>%
  mutate(
    Frame_0 = Ct_fr0 / (Ct_fr0 + Ct_fr1 + Ct_fr2),
    Frame_1 = Ct_fr1 / (Ct_fr0 + Ct_fr1 + Ct_fr2),
    Frame_2 = Ct_fr2 / (Ct_fr0 + Ct_fr1 + Ct_fr2)
  ) %>%
  select(samplez, gene, Frame_0, Frame_1, Frame_2) %>%
  pivot_longer(names_to = "Frame", cols = 3:5)

### Position specific normalized reads
pos_sp_df <- lapply(sample_names, function(x) {
  # path to the read length file
  file_loc <- normalizePath(here(x, "pos_sp_rpf_norm_reads.tsv"))
  # read it in
  return(read_tsv(file_loc, skip = 4))
}) %>%
  bind_rows(.id = "samplez")

### collated TPMs
collated_tpms_df <-
  read_tsv(normalizePath(here(yaml$dir_out, "TPMs_collated.tsv")), skip = 4)

### Ribogrid
# ribogrid_df <- lapply(sample_names, function(x) {
#   # path to the read length file
#   file_loc <- paste0(x, "/gene_position_length_counts_5start.tsv")
#   # read it in
#   return(read_tsv(file_loc, skip = 4))
# }) %>%
#   bind_rows(.id = "samplez")


### Gene features
# this is commented out because we're not sure what the format of the sequence features file is
# presumably it's to match the one in the vignette but I'm not sure if it's described anywhere.
if (any(names(yaml) == "features_file")) {
  if (length(yaml$features_file) > 0) {
    features_df <- lapply(sample_names, function(x) {
      file_loc <- paste0(x, "/sequence_features.tsv")

      return(read_tsv(file_loc, skip = 4))
    }) %>%
      bind_rows(.id = "samplez")

    possible_features <- unique(features_df$Feature)
  }
}

###############################################################################
# Plotting
###############################################################################

### define the universal plot theme
plot_theme <- theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    text = element_text(size = 15)
  )

### Start the server
server <- function(input, output, session) {
  ######################################
  # Read count summary
  ######################################
  
  ### bar plot
  output$read_counts_bar_plot <- renderPlot({
    read_counts_df %>%
      filter(SampleName %in% input$sample) %>%
      ggplot(., aes(SampleName, NumReads, fill = short_name)) +
      geom_col(position = position_dodge()) +
      plot_theme +
      scale_fill_discrete(name = NULL) +
      labs(x = NULL,
           y = "Read count",
           title = "Read counts per step") +
      guides(fill = guide_legend(ncol = 1))+
      scale_y_continuous(labels = label_number_si())
  })
  
  ######################################
  # TPM summary
  ######################################
  
  ### Sample to sample cors
  # load the TPM df
  cor.df <- collated_tpms_df %>%
    column_to_rownames("ORF")
  
  # get a fixed column order, alphabetical in this case
  cor.cols <- sort(names(cor.df))
  
  # create correlations, fix col order
  cor.df <- cor.df %>%
    select(all_of(cor.cols)) %>%
    apply(., 2, log10) %>%
    na_if(., -Inf) %>%
    cor(use = "pairwise.complete.obs")
  
  # remove lower triangle
  cor.df[lower.tri(cor.df)] <- NA
  
  # convert to a df and reshape
  cor.df2 <- as_tibble(cor.df, rownames = "samp1") %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "samp2",
      values_to = "R"
    ) %>%
    filter(samp1 != samp2) %>%
    mutate(
      samp1 = factor(samp1, levels = cor.cols),
      samp2 = factor(samp2, levels = cor.cols)
    )
  
  # find the floor of min correlation
  min.cor <- floor(min(cor.df2$R, na.rm = TRUE) * 10) / 10
  
  # plot it
  output$sample_cors_plot <- renderPlot({
    ggplot(cor.df2, aes(samp1, samp2, fill = R, label = signif(R, 2))) +
      geom_raster() +
      scale_fill_viridis_c(
        option = "B",
        na.value = NA,
        limits = c(min.cor, 1)
      ) +
      theme(
        panel.background = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1
        ),
        panel.grid.major.x = element_line(color = "grey50", linetype = 3),
        axis.ticks = element_blank()
      ) +
      labs(x = NULL,
           y = NULL,
           title = expression(paste(
             "Pairwise correlations based on ", log[10], "(TPM)"
           ))) +
      geom_text()
  })
  
  ### Abundance of specific gene compared to all
  output$tpm_density_plot <- renderPlot({
    collated_tpms_df %>%
      pivot_longer(where(is.numeric),
                   names_to = "samplez",
                   values_to = "tpm") %>%
      filter(samplez %in% input$sample) %>%
      mutate(
        labz = ifelse(ORF %in% input$gene, ORF, NA),
        xint = ifelse(ORF %in% input$gene, tpm, NA)
      ) %>%
      ggplot(., aes(x = tpm, fill = samplez)) +
      geom_density() +
      geom_vline(aes(xintercept = xint), size = 1) +
      scale_x_log10(labels = label_number_si()) +
      facet_wrap(~ samplez) +
      scale_fill_discrete(guide = FALSE) +
      plot_theme +
      labs(title = "TPM distributions",
           x = expression(paste(log[10], "(TPM)")),
           y = "Density")
  })
  
  ######################################
  # Read length distributions
  ######################################
  
  ### Read length distribution plot
  output$read_length_dist_plot <- renderPlot({
    read_length_df %>%
      filter(samplez %in% input$sample) %>%
      ggplot(., aes(Length, Counts, fill = samplez)) +
      geom_col(width = 1) +
      facet_wrap(~ samplez, scales = "free_y") +
      scale_fill_discrete(name = NULL, guide = FALSE) +
      plot_theme +
      labs(title = "Read length distributions",
           x = "Read length",
           y = "Read counts") +
      scale_x_continuous(breaks = seq(10, 50, 5)) +
      scale_y_continuous(breaks = breaks_pretty(n = 4),
                         labels = label_number_si())
  })
  
  ######################################
  # 3nt periodicity
  ######################################
  
  ### Periodicity metagene plot for both 5' and 3' end as a line plot
  output$periodicity_line_plot <- renderPlot({
    periodicity_df %>%
      filter(samplez %in% input$sample) %>%
      ggplot(., aes(x = Pos, y = Counts, color = samplez)) +
      geom_line() +
      scale_color_discrete(name = NULL, guide = FALSE) +
      facet_grid(cols = vars(samplez),
                 rows = vars(End),
                 scales = "free_y") +
      labs(title = "3-nucleotide periodicity",
           x = "Codon position",
           y = "Read counts") +
      plot_theme+
      scale_y_continuous(labels = label_number_si())
  })
  
  ### Read counts per frame as a bar plot
  output$frame_counts_bar_plot <- renderPlot({
    frame_bar_df %>%
      filter(samplez %in% input$sample) %>%
      ggplot(., aes(x = Frame, y = Count, fill = samplez)) +
      geom_bar(stat = "identity", pos = "dodge") +
      scale_fill_discrete(name = NULL) +
      labs(title = "Read counts per frame",
           x = NULL,
           y = "Read counts") +
      plot_theme +
      scale_color_discrete(guide = FALSE)+
      scale_y_continuous(labels = label_number_si())
  })
  
  ### Proportion of reads per frame for each gene as a boxplot
  output$frame_proportions_boxplot <- renderPlot({
    frame_box_df %>%
      filter(samplez %in% input$sample) %>%
      ggplot(., aes(x = Frame, y = value, fill = samplez)) +
      labs(y = "Proportion",
           title = "Proportion of reads per frame per gene",
           x = NULL) +
      scale_fill_discrete(name = NULL) +
      geom_boxplot() +
      plot_theme
  })
  
  ######################################
  # Position specific normalized reads
  ######################################
  
  ### 3' end
  output$pos_sp_plot3 <- renderPlot({
    pos_sp_df %>%
      filter(End == "3'" & samplez %in% input$sample) %>%
      ggplot(., aes(Position, Mean, color = samplez)) +
      geom_line() +
      facet_grid(End ~ samplez, scales = "free") +
      plot_theme +
      scale_color_discrete(guide = FALSE)
  })
  
  ### 5' end
  output$pos_sp_plot5 <- renderPlot({
    pos_sp_df %>%
      filter(End == "5'" & samplez %in% input$sample) %>%
      ggplot(., aes(Position, Mean, color = samplez)) +
      geom_line() +
      facet_grid(End ~ samplez, scales = "free") +
      plot_theme +
      scale_color_discrete(guide = FALSE)
  })
  
  ######################################
  # Ribogrids
  ######################################
  
  # ### Ribogrid plot
  # output$ribogrid_plot <- renderPlot({
  #   ribogrid_df %>%
  #     filter(
  #       samplez %in% input$sample &
  #         ReadLen >= min(input$ribogrid_len_range) &
  #         ReadLen <= max(input$ribogrid_len_range) &
  #         Pos >= min(input$ribogrid_pos_range) &
  #         Pos <= max(input$ribogrid_pos_range)
  #     ) %>%
  #     ggplot(., aes(Pos, ReadLen, fill = Counts)) +
  #     geom_raster() +
  #     scale_fill_gradient(low = "white", high = "navy") +
  #     plot_theme +
  #     labs(title = "Ribogrid", x = "Position of 5' end of read", y = "Read length") +
  #     facet_wrap( ~ samplez)
  # })
  #
  # ### Ribogrid bar plot
  # output$ribogridbar_plot <- renderPlot({
  #   ribogrid_df %>%
  #     filter(
  #       samplez %in% input$sample &
  #         ReadLen >= min(input$ribogrid_len_range) &
  #         ReadLen <= max(input$ribogrid_len_range) &
  #         Pos >= min(input$ribogrid_pos_range) &
  #         Pos <= max(input$ribogrid_pos_range)
  #     ) %>%
  #     ggplot(., aes(Pos, Counts)) +
  #     geom_col() +
  #     facet_grid(ReadLen ~ samplez) +
  #     plot_theme +
  #     labs(x = "Position of 5' end of read", y = "Read count", title = "Ribogrid (bar)")
  # })
  
  ######################################
  # Single-gene ribogrids
  ######################################
  
  ### Single-gene ribogrid
  # Go through h5 file for each sample and in a specific gene
  sgribogrid_df <- reactive({
    lapply(sample_names, function(x) {
      only_name <- str_split(x, "/")
      
      nm <- only_name[[length(only_name)]]
      
      file_loc <-
        normalizePath(here(yaml$dir_out, paste0(nm[[length(nm)]], "/", nm[[length(nm)]], ".h5")))
      
      ret <-
        h5read(file_loc, name = file.path(input$gene, yaml$dataset, "reads/data"))$data %>%
        as_tibble() %>%
        rowid_to_column("read_length") %>%
        pivot_longer(cols = starts_with("V"), names_to = "position") %>%
        mutate(position = as.numeric(str_remove(position, "V")))
      
      h5closeAll()
      
      return(ret)
    }) %>%
      bind_rows(.id = "samplez")
  })
  
  output$sgribogrid_plot <- renderPlot({
    sgribogrid_df %>%
      filter(
        samplez %in% input$sample &
          read_length >= min(input$ribogrid_len_range) &
          read_length <= max(input$ribogrid_len_range) &
          position >= min(input$ribogrid_pos_range) &
          position <= max(input$ribogrid_pos_range)
      ) %>%
      ggplot(., aes(position, read_length, fill = value)) +
      geom_raster() +
      facet_wrap(~samplez, ncol = 1) +
      scale_fill_gradient(low = "white",
                          high = "navy",
                          name = "Read count") +
      theme_bw() +
      theme(
        text = element_text(size = 14),
        panel.grid = element_blank()
      ) +
      labs(x = "Position",
           y = "Read length",
           title = input$gene)+
      scale_x_continuous(breaks = breaks_pretty(n = 8))
  })
  
  
  ######################################
  # Features
  ######################################
  
  # Features are an optional file if it exists it's here
  # but the format isn't described so I'll hold off on this for now
  
  ### TPM vs other features
  if (any(names(yaml) == "features_file") &&
      length(yaml$features_file) > 0) {
    output$features_plot <- renderPlot({
      if (input$gene2 == "") {
        features_df %>%
          filter(samplez %in% input$sample &
                   Feature %in% input$feature) %>%
          ggplot(., aes(tpm, Value)) +
          geom_point() +
          geom_smooth(method = "lm") +
          facet_wrap(samplez ~ Feature, scales = "free") +
          scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                                     .x))) +
          plot_theme +
          labs(x = expression(paste(log[10], "(TPM)"))) +
          scale_size_manual(values = c(4, 1), guide = FALSE)
      } else {
        features_df %>%
          filter(samplez %in% input$sample &
                   Feature %in% input$feature) %>%
          mutate(labz = ifelse(ORF == input$gene2, "a_label", "no_lab")) %>%
          arrange(labz) %>%
          ggplot(., aes(tpm, Value)) +
          geom_point(aes(color = labz, size = labz)) +
          geom_smooth(method = "lm") +
          facet_wrap(samplez ~ Feature, scales = "free") +
          scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                                     .x))) +
          plot_theme +
          scale_color_manual(guide = FALSE,
                             values = c("firebrick3", "grey70")) +
          labs(x = expression(paste(log[10], "(TPM)"))) +
          scale_size_manual(values = c(4, 1), guide = FALSE)
      }
    })
  } else {
    output$features_plot <- renderPlot({
      # this is simply a dummy plot to plot nothing
      ggplot(mtcars, aes(x = wt, y = mpg)) + 
        geom_blank() +
        theme_void()
    })
  }
  
}

#
# END SERVER FUNCTION
#

###############################################################################
# UI design
###############################################################################

ui <- fluidPage(# App title
  headerPanel("riboviz2"),
  
  # Sidebar panel for inputs
  fluidRow(column(
    12,
    wellPanel(
      checkboxGroupInput(
        inputId = "sample",
        label = "Sample:",
        choiceNames = names(sample_names),
        choiceValues = names(sample_names),
        selected = names(sample_names)[1],
        inline = TRUE,
        width = "100%"
      ),
      submitButton(text = "Apply Changes", icon = NULL, width = NULL)
    )
    
  )),
  
  # Main panel for displaying outputs
  mainPanel(width = 12,
            fluidRow(
              tabsetPanel(
                tabPanel(
                  "Read count summary",
                  headerPanel(""),
                  plotOutput("read_counts_bar_plot")
                ),
                
                tabPanel(
                  "TPM summary",
                  wellPanel(textInput(
                    inputId = "gene",
                    label = "Gene:",
                    width = "50%"
                  )),
                  headerPanel(""),
                  plotOutput("sample_cors_plot"),
                  headerPanel(""),
                  plotOutput("tpm_density_plot")
                ),
                
                tabPanel(
                  "Read length distributions",
                  headerPanel(""),
                  plotOutput("read_length_dist_plot")
                ),
                
                tabPanel(
                  "3nt periodicity",
                  headerPanel(""),
                  plotOutput("periodicity_line_plot"),
                  headerPanel(""),
                  plotOutput("frame_counts_bar_plot"),
                  headerPanel(""),
                  plotOutput("frame_proportions_boxplot")
                ),
                
                tabPanel(
                  "Normalized reads by position",
                  headerPanel(""),
                  plotOutput("pos_sp_plot5"),
                  headerPanel(""),
                  plotOutput("pos_sp_plot3")
                ) ,
                
                tabPanel(
                  "Single-gene ribogrid",
                  wellPanel(textInput(
                    inputId = "gene",
                    label = "Gene:",
                    width = "50%"
                  )),
                  headerPanel(""),
                  sliderInput(
                    "ribogrid_len_range",
                    "Ribogrid read length range:",
                    min = 0,
                    max = 50,
                    value = c(20, 35),
                    width = "50%"
                  ),
                  sliderInput(
                    "ribogrid_pos_range",
                    "Ribogrid read position range:",
                    min = 0,
                    max = 100,
                    value = c(-10, 10),
                    width = "50%"
                  ),
                  headerPanel(""),
                  plotOutput("sgribogrid_plot")
                ),
                
                # tabPanel(
                #   "Ribogrid",
                #   sliderInput(
                #     "ribogrid_len_range",
                #     "Ribogrid read length range:",
                #     min = 0,
                #     max = 50,
                #     value = c(20, 35),
                #     width = "50%"
                #   ),
                #   sliderInput(
                #     "ribogrid_pos_range",
                #     "Ribogrid read position range:",
                #     min = -24,
                #     max = 50,
                #     value = c(-10, 10),
                #     width = "50%"
                #   ),
                #   headerPanel(""),
                #   plotOutput("ribogrid_plot"),
                #   headerPanel(""),
                #   plotOutput("ribogridbar_plot", height = "1000px")
                # ) #,
                
                tabPanel(
                  "Features",
                  wellPanel(textInput(
                    inputId = "gene2",
                    label = "Gene:",
                    width = "50%"
                  ),
                  if (any(names(yaml) == "features_file") &&
                      length(yaml$features_file) > 0) {
                    checkboxGroupInput(
                      inputId = "feature",
                      label = "Feature:",
                      choiceNames = possible_features,
                      choiceValues = possible_features,
                      selected = possible_features[1],
                      inline = TRUE,
                      width = "100%"
                    )
                  }),
                  headerPanel(""),
                  plotOutput("features_plot", height = "1000px")
                )
              )
            )))

# I cast spell, ~run shiny app~
shinyApp(ui, server, options = list(port = 4254))