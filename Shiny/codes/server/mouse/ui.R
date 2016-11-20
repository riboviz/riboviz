library(shiny)
library(rhdf5)

#-------------------------------------------- part2:define webpage components and layout

shinyUI<-fluidPage(
  titlePanel('Explore mouse RPF databases'),
  sidebarLayout(
    sidebarPanel(
	  #select input with the list of datasets
      selectInput(inputId = 'selecttissue',label='Select the tissue',choices=c('Choose','Liver','Muscle','Retina','nonChx-Liver')),
      
      uiOutput('vsample'), #3 dynamic drop down menus
      uiOutput('vsampletype'),
      
      actionButton(inputId = 'clickadd',label='Add to selected database (max 5)'),
      br(),
      br(),
      actionButton(inputId = 'clickempty',label='Reset database selections'),
      br(),
      br(),
      #------------above is to select databases,below is to select plot type     
      textInput(inputId = 'txt',label='Write the gene name of interest',value='ENSMUST00000000049.5', placeholder = 'e.g. ENSMUST00000000049.5 or 10006'),
      selectInput(inputId = 'select',label='Select a type of plot',
                  choices = c('Choose','Reads by length'='rbl','Reads by position'='rbp','Reads by position and length'='rbpl','Reads by codon'='rbc')),
      conditionalPanel(
        condition = "input.select == 'rbpl'",
        textInput(inputId = 'txt2',label='Write the read length of interest',placeholder='ranging from 15-50')),
      actionButton(inputId = 'clickplot',label='Click to plot'),
      downloadButton(outputId = 'down',label='Download data')
    ),
    mainPanel(
      htmlOutput('sdb'),
      plotOutput("plot1")
    )
  )
)
