library(shiny)
library(rhdf5)
library(ggplot2)
library(reshape2)
library(plotly)
library(RColorBrewer)
library(data.table) # To use rbindlist()
library(tidyr) # To use gather()



load("./riboViz.RData")

shinyUI <- fixedPage(
  fixedRow(
  # Where selection for databases and genes are made
    column(2,
         selectInput("slct_gene", "Genes",
                    names(d2)[6:length(names(d2))], 
                    selected = 'YAL001C')),
     column(4,    
          selectizeInput("slct_db", "Databases",
                         paste(df$year, df$author, df$database, sep = " "), 
                         multiple = TRUE, 
                         selected = "2016 Weinberg unselected_total_RNA",
                         options = list(maxItems=9))),

     column(5, 
       br(),
       htmlOutput("txt_legends"))              
    ),

  fixedRow(
    column(12,
           htmlOutput("txt_test"),
           textOutput("txt_errmsg"),
           tags$head(tags$style("#txt_errmsg{color: red}"))
    )
  ),
  
  ###################################### Conditional Panel
    conditionalPanel(
      condition="input.slct_db != null", # Plots only show up when database selection is not empty
    fixedRow(            # Where the 1st plot appears
      column(4,
             div(h3("Reads by Length"),style = "color:grey"),
             div(p("The distribution of reads of specific lengths along the ORF."),
                 style = "color:grey")
             
      ),
      column(8,
             plotlyOutput("plotrbl"))#,width = 600, height = 350))
      ),
  
    br(),
    br(),

    fixedRow(
      column(4,offset = 4,
      selectInput("slct_readlen", "Select the read length of interest",
                  c("All read lengths",15:50), 
                  selected = "All read lengths")),
      column(4, 
      radioButtons(inputId = 'btn_rbpl',label = "Select read type", 
                   choices = list("Regular reads" = 1, "Normalized reads" = 2), 
                   selected = 1, inline=TRUE))),

    fixedRow(  
    column(4,
           div(h3("Reads by Position and Length"),style = "color:grey"),
           div(p("Gene/read length/position-specific distribution of mapped reads along the 5' UTR, ORF and 3' UTR lengths.
                 The ORF is indicated by the shaded grey area. The normalized reads are expressed as the proportion of actual
                 reads and the mean within each ORF."),
               style = "color:grey")
    ),
    column(8,
           plotlyOutput("plotrbpl"))
    ),
  
    br(),
    br(),

    fixedRow(
      column(4, offset=8,
             radioButtons(inputId = 'btn_rbc',label = "Select read type", 
                          choices = list("Regular reads" = 1, "Normalized reads" = 2), 
                          selected = 1, inline=TRUE))),
    
    fixedRow(  
    column(4,
           div(h3("Reads by Codon"),style = "color:grey"),
           div(p("Codon specific mapped reads. The normalized reads are expressed as the proportion of total 
                 reads for corresponding nucleotides, and the mean for all codons."), style = "color:grey")
    ),
    column(8,
           plotlyOutput("plotrbc"))
    ),
    
    br(),
    br(),
    fixedRow(
      column(8, offset=4,
             sliderInput(inputId = 'bin_rpkm',label='Choose bins for the histogram',min=5,max=30, value = 10))
    ),
  
    br(),
    br(),
  
    fixedRow( 
      column(4,
             div(h3("RPKM"),style = "color:grey"),
             div(p("RPKM stands for Reads Per Kilobase of transcript per Million mapped reads. This plot presents
                   the overall abundance of specific gene relative to its abundance in a curated set of wild-type datasets."),
                 style = "color:grey")
      ),
      column(8,
             plotlyOutput("plotrpkm"))),
    
    
    fixedRow( 
      column(4,
             div(h3("RPKM Density"),style = "color:grey"),
             div(p("The abundance of specific gene relative to the abundances of all other genes in the same data set."),
                 style = "color:grey")
      ),
      column(8,
             plotlyOutput("plotrpkmdens"))),
    
    br(),
    br(),
    
    fixedRow(
      column(4,
            div(h3("Download"),style = "color:grey"), 
            br()),
            
      column(8,
             br(),
             downloadButton(outputId = 'dwld_data',label='Download all data')
            )
      ),
    br()
    )
  )
