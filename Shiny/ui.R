library(shiny)
library(rhdf5)
library(ggplot2)
library(reshape2)
library(cowplot)

#For future database update, you need an updated file listing all year,author and database name like 'data_unq.tsv',

shinyUI<-fixedPage(
  tags$head(
    tags$style(HTML('#clickadd{background-color:orange}'))
  ),
  options = list(scrollX = TRUE), #scroll bar for the page
  #titlePanel('Explore the RPF and mRNA databases'),
  
  fixedRow(
    column(4, #id='sidebar',
           #select input with the list of datasets
           selectInput(inputId = 'selectyear',label='Select the year of publication',choices=c(2016,2015,2014,2013,2012,2009)),
           
           uiOutput('vauthor'), #3 dynamic drop down menus
           uiOutput('vdb'),
           
           actionButton(inputId = 'clickadd',label='Add to selected database (max 3)'),
           actionButton(inputId = 'go',label='Go'), 
           br(),
           br(),
           actionButton(inputId = 'clickempty',label='Reset database selections'),
           br(),
           br(),
           
           textOutput("dbselect_text"),
           tags$head(tags$style("#dbselect_text{color: red;   #manage the font color and style of the output text
                                #font-size: 20px;
                                #font-style: italic;
                                }"
        )
           ), 
        htmlOutput('sdb'),
        br(),
        
        #------------above is to select databases,below is to select plot type 
        conditionalPanel(
          condition='input.go!=0 & input.clickadd!=0',
          
          textInput(inputId = 'txt',label='Write the gene name of interest',value = 'YAL001C', placeholder = 'e.g. YAL001C'),
          selectInput(inputId = 'select',label='Select a type of plot',
                      choices = c('Reads by length'='rbl','Reads by position'='rbp','Reads by position and length'='rbpl','Reads by codon'='rbc','RPKM comparison with WT'='rpkm')),
          conditionalPanel(
            condition = "input.select == 'rbpl'",
            textInput(inputId = 'txt2',label='Write the read length of interest',placeholder='ranging from 15-50')
          ),
          conditionalPanel(
            condition = "input.select == 'rpkm'",
            sliderInput(inputId = 'bins',label='Choose bins for the histogram',min=5,max=30, value = 10)),
          conditionalPanel(
            condition="input.select=='rbp' | input.select=='rbpl' | input.select=='rbc'",
            radioButtons(inputId = 'radioreads',label = "Select read type", choices = list("Regular reads" = 1, "Normalized reads" = 2), selected = 1)
          ),
          actionButton(inputId = 'clickplot',label='Click to plot'),
          downloadButton(outputId = 'down',label='Download data')
        )
           ),
    column(8,
           textOutput('errormsg'),
           tags$head(tags$style("#errormsg{color: red}")), 
           
           #---------below is legend for the plots 
           textOutput('l1'),
           tags$head(tags$style("#l1{color: rgb(67,162,202);font-size: 18px}")),
           textOutput('l2'),
           tags$head(tags$style("#l2{color: rgb(254,178,76);font-size: 18px}")),
           textOutput('l3'),
           tags$head(tags$style("#l3{color: rgb(136,86,167);font-size: 18px;}")),
           
           plotOutput('plot11',
                      dblclick = "plot11_dblclick",
                      brush = brushOpts(
                        id = "plot11_brush",
                        resetOnNew = TRUE),
                      inline = TRUE),
           #width = 1000,height=500),
           
           plotOutput('plot22',inline = TRUE)
           #width = 1000,height=500)
           
           #uiOutput("plots",width = 1000,height=500)
    )
           )
  )