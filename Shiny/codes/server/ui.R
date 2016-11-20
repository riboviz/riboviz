library(shiny)
library(rhdf5)
#-------------------------------------------- part2:define webpage components and layout
shinyUI<-fluidPage(
#	df<-read.table('Supp/data_unq.tsv'); #stores the year/author/database information
#	names(df) <- c('year','author','journal','NA','database','NA','NA','NA','dbtype','condition','geoID')
#	year_vector<-as.vector(unique(df[,1])); #vector containing all years
#	author_vector<-as.vector(unique(df[,2])); #vector containing all authors from all years
#	max_db <- 5; #set the maximum databases to be selected 

  titlePanel('Explore the RPF and mRNA databases'),
  sidebarLayout(
    sidebarPanel(
      
	  #select input with the list of datasets
      selectInput(inputId = 'selectyear',label='Select the year of publication',choices=c('Choose',2009:2016)),
      
      uiOutput('vauthor'), #3 dynamic drop down menus
      uiOutput('vrpf'),
      uiOutput('vmrna'),
      
      actionButton(inputId = 'clickadd',label='Add to selected database (max 5)'),
      actionButton(inputId = 'clickempty',label='Reset database selections'),
      
#------------above is to select databases,below is to select plot type     
      textInput(inputId = 'txt',label='Write the gene name of interest',placeholder = 'e.g. YAL001C'),
      selectInput(inputId = 'select',label='Select a type of plot',
                  choices = c('Choose','Reads by length'='rbl','Reads by position'='rbp','Reads by position and length'='rbpl','Reads by codon'='rbc')),
      conditionalPanel(
        condition = "input.select == 'rbpl'",
        textInput(inputId = 'txt2',label='Write the read length of interest',placeholder='ranging from 15-50')),
      actionButton(inputId = 'clickplot',label='Click to plot')
    ),
    mainPanel(
      textOutput('sdb'),
      plotOutput("plot1")
    )
  )
)

