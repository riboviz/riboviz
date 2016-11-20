
library(shiny)
library(ggplot2)

path<-'/home/txing/shinyProject/data/RPKM/'
d1<-read.delim(paste0(path,'F8_RPKMs_modified.tsv'))
d2<-read.delim(paste0(path,'bgdb.tsv'))

year_vector<-as.vector(sort(unique(d1$Year),decreasing=TRUE)) #vector containing all years
author_vector<-as.vector(unique(d1$Author)) #vector containing all authors from all years

#------------------------------------------------
myfun1<-function(str_year){ #returns the author selections for UI year
  as.character(unique(d1[d1$Year==str_year,2])) 
}

myfun2<-function(str_author,str_year,str_dbtype){ #returns the database selections for UI year/author/ for either mRNA or RPF databases
  as.character(d1[(d1$Year==str_year)&(d1$Author==str_author)&(d1$Type==str_dbtype),3])
}

#------------------------------------------------
shinyUI<-fluidPage(
  titlePanel('Compare RPKM from experiments using WT and other conditions'),
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = 'gene',label='Write the gene name of interest',value='YAL001C',placeholder = 'e.g.YAL001C'),
      # actionButton(inputId = 'gorpf',label='See RPF background'),
      # actionButton(inputId = 'gomrna',label='See mRNA background'),
      # br(),
      # br(),
      selectInput(inputId = 'selectyear',label='Select the year of publication',choices=c('Choose',year_vector)),
      uiOutput('vauthor'), #3 dynamic drop down menus
      uiOutput('vrpf'),
      uiOutput('vmrna'),
      sliderInput(inputId = 'bins',label='Choose bins for the histogram',min=5,max=30, value = 10),
      actionButton(inputId = 'goplot',label='Click to plot the designated database')
    ),
    mainPanel(
      textOutput('errormsg'),
      plotOutput('plot')
    )
  )
)
#-------------------------------------------- part3: shinyServer
shinyServer<-function(input,output){
  var1<-reactive({ #returns a list of author name choices corresponding to UI year
    i <- which(year_vector==input$selectyear)
    myfun1(year_vector[i])})
  
  var2<-reactive({ #returns a list of RPF database names corresponding toUI year and UI author
    i<-which(author_vector==input$var_author)
    myfun2(author_vector[i],input$selectyear,'RPF')})
  
  var3<-reactive({ #returns a list of mRNA database names corresponding to UI year and UI author
    i<-which(author_vector==input$var_author)
    myfun2(author_vector[i],input$selectyear,'mRNA')})
  
  output$vauthor<-renderUI({ #create the dropdown menu for author selection
    selectInput('var_author','Select the author',choices=c('Choose',var1()))})
  
  output$vrpf<-renderUI({ #create the dropdown menu for RPF databse selection
    selectInput('var_rpf','Select the RPF database',choices=c('Choose',var2()))})
  
  output$vmrna<-renderUI({ #create the dropdown menu for mRNA databse selection
    selectInput('var_mrna','Select the mRNA database',choices=c('Choose',var3()))})
  
  #-------------------------plot the RPKM histogram for background database
  observeEvent(
    input$goplot,{
      output$plot<-renderPlot({
        input$goplot
        isolate({
          i<-which(names(d2)==input$gene) #indexes for same gene in d1 and d2 are the same
          if (input$var_rpf=='Choose' & input$var_mrna!='Choose'){
            vline<-d1[d1$Year==input$selectyear & d1$Author==input$var_author & d1$Dataset==input$var_mrna,i]
            dat<-d2[d2$Type=='mRNA',c(5,i)]
            ggplot(dat, aes_string(x=input$gene)) +
              geom_histogram(bins = input$bins, position='dodge',alpha=0.5,color='blue',fill='blue') +
              theme(legend.position='top') +
              geom_vline(xintercept = vline,color='red')
          }
          else if (input$var_rpf!='Choose' & input$var_mrna=='Choose'){
            vline<-d1[d1$Year==input$selectyear & d1$Author==input$var_author & d1$Dataset==input$var_rpf,i]
            dat<-d2[d2$Type=='RPF',c(5,i)]
            ggplot(dat, aes_string(x=input$gene)) +
              geom_histogram(bins = input$bins, position='dodge',alpha=0.5,color='blue',fill='blue') +
              theme(legend.position='top') +
              geom_vline(xintercept = vline,color='red')
          }
          else if (input$var_rpf!='Choose' & input$var_mrna!='Choose'){
            output$errormsg<-renderText({print('!!!Please select either a mRNA or RPF database and de-select the other by selecting "Choose"')})
          }
        })  
      })
    })
}
shinyApp(ui = shinyUI, server = shinyServer)