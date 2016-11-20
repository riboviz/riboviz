library(shiny)
library(rhdf5)
#For future database update, you need an updated file listing all year,author and database name like 'data_unq.tsv',

path<-'/data/riboseq/'
df<-read.table('Supp/data_unq.tsv',sep='') #stores the year/author/database information
names(df)=c('year','author','journal','NA','database','NA','NA','NA','dbtype','condition','geoID')
year_vector<-as.vector(unique(df$year)) #vector containing all years
author_vector<-as.vector(unique(df$author)) #vector containing all authors from all years
max_plots <- 5

#-------------------------------------------- part1
myfun1<-function(str_year){
  as.character(unique(df[df$year==str_year,2]))
}

myfun2<-function(str_author,str_year,str_dbtype){
  as.character(df[(df$year==str_year)&(df$author==str_author)&(df$dbtype==str_dbtype),5])
}
#-------------------------------------------- part2
shinyUI<-fluidPage(
  titlePanel('Explore the RPF databases'),
  sidebarLayout(
    sidebarPanel(
      
      #select input with the list of datasets
      selectInput(inputId = 'selectyear',label='Select the year of publication',choices=c('Choose',year_vector)),
      
      uiOutput('vauthor'),
      uiOutput('vrpf'),
      uiOutput('vmrna'),
      
      actionButton(inputId = 'clickadd',label='Add to selected database (max 5)'),
      actionButton(inputId = 'clickempty',label='Reset database selections'),
      #------------above is to select databases,below is to select plot type     
      textInput(inputId = 'txt',label='Write the gene name of interest',value='YAL001C',placeholder = 'e.g. YAL001C'),
      selectInput(inputId = 'select',label='Select a type of plot',
                  choices = c('Choose','Reads by length'='rbl','Reads by position'='rbp','Reads by position and length'='rbpl','Reads by codon'='rbc')),
      conditionalPanel(
        condition = "input.select == 'rbpl'",
        textInput(inputId = 'txt2',label='Write the read length of interest',value='15',placeholder='ranging from 15-50'))
      #actionButton(inputId = 'clickplot',label='Click to plot')
    ),
    mainPanel(
      textOutput('sdb'),
      uiOutput("plots")
    )
  )
)
#-------------------------------------------- part3
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
  
  output$vauthor<-renderUI({
    selectInput('var_author','Select the author',choices=c('Choose',var1()))})
  
  output$vrpf<-renderUI({
    selectInput('var_rpf','Select the RPF database',choices=c('Choose',var2()))})
  
  output$vmrna<-renderUI({
    selectInput('var_mrna','Select the mRNA database',choices=c('Choose',var3()))})
  
  values<-reactiveValues()
  values$selectdb<-vector()#create a reactive value to store the UI selectInputs
  a<-eventReactive(input$clickadd,{if((input$var_rpf!='Choose')&(input$selectyear!='Choose')&(input$var_author!='Choose')){paste(input$selectyear,input$var_author,input$var_rpf,sep='/')}})
  b<-eventReactive(input$clickadd,{if((input$var_mrna!='Choose')&(input$selectyear!='Choose')&(input$var_author!='Choose')){paste(input$selectyear,input$var_author,input$var_mrna,sep='/')}})
  
  
  observe({
    values$selectdb<-rbind(isolate(values$selectdb),a(),b())
  })
  
  observe({ #when reset button(input$clickempty) is pressed, values$selectdb is reset to empty
    if (input$clickempty==0) return()
    values$selectdb<-NULL
  })
  
  temp<-reactive({tail(unique(values$selectdb),5)}) #stores the last 5 UI databases
  output$sdb<-renderPrint(if(length(values$selectdb)>0){print(temp())})
  
#-----------------------------------below is to plot
  # Insert the right number of plot output objects into the web page
  output$plots <- renderUI({
    plot_output_list <- lapply(1:length(temp()), function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname, height = 320, width = 700)
    })
    
  do.call(tagList, plot_output_list)# Convert the list to a tagList - this is necessary for the list of items to display properly.
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they are visible on the web page.
  for (i in 1:max_plots) {
    # Need local so that each item gets its own number. Without it, the value of i in the renderPlot() will be the same across all instances, because of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("plot", my_i, sep="")
      
      output[[plotname]] <- renderPlot({
        tempi<-reactive({temp()[my_i]})
        geoid<-reactive({as.character(df[(df$year==unlist(strsplit(tempi(),'/'))[1])&(df$author==unlist(strsplit(tempi(),'/'))[2])&(df$database==unlist(strsplit(tempi(),'/'))[3]),11])})
        inf<-reactive({paste(path,tempi(),'/',geoid(),'.h5',sep='')}) #e.g. temp()[1]='2016/Weinberg/RPF'
        transtemp<-reactive({paste(unlist(strsplit(tempi(),'/'))[1],unlist(strsplit(tempi(),'/'))[2],unlist(strsplit(tempi(),'/'))[3],sep='_')}) #e.g. transtemp='2016_Weinberg_RPF'
        attrloc<-reactive({h5readAttributes(file=inf(),name=paste(input$txt,'/',transtemp(),'/reads',sep=''))}) #e.g. name='YAL001C/2016/Weinberg/RPF/reads'
        group<-reactive({paste(input$txt,'/',transtemp(),'/reads/data',sep='')}) #e.g. name='YAL001C/2016/Weinberg/RPF/reads/data'
        coord<-reactive({dim(h5read(inf(),group()))[2]}) #nucleotide coordinate e.g.3980 for YAL001C
        endcoord<-reactive({coord()-247}) #end coordinate for the last nucleotide in the ORF
        d<-reactive({h5read(inf(),group())}) #e.g. for YAL001C is a 36x3980 matrix
        sumd<-reactive({lapply(data.frame(d()),sum)}) #stores reads for all length, all positions
        
        d1<-reactive({d()[14:15,(251-15):(endcoord()-15)]}) #matrix 2 rows (ribosome protected length=28,29),columns=endcoord-251+1
        d2<-reactive({d()[16,(251-16):(endcoord()-16)]}) #matrix 1 row (ribosome protected length=30)
        cd<-reactive({rbind(d1(),d2())}) #3X7482, matrix to calculate the ribosome density
        num_of_codons<-reactive({(endcoord()-251+1)/3}) #number of codons= integer
        cd2<-reactive({matrix(cd(),nrow=nrow(cd())*3,ncol=num_of_codons())})
        newcd<-reactive({lapply(data.frame(cd2()),sum)})
        
        if(input$select=='rbl'){
          # plot({attrloc()$lengths},{attrloc()$reads_by_len},type='l',main='Read counts by Length',xlab='Read length',ylab='Read count')}
          barplot(attrloc()$reads_by_len,main='Read counts by Length',xlab='Read length',ylab='Read count',names.arg = c(15:50),col='blue')
          }
        
        else if(input$select=='rbp'){
          plot(1:coord(),sumd(),type='l',main='Read counts by Position',xlab='Read position',ylab='Read count')}
        
        else if(input$select=='rbpl'){
          plot(1:coord(),d()[(as.numeric(input$txt2)-14),],type='l',main='Read counts by Position and length',xlab='Read position',ylab='Read count')}
        
        else if(input$select=='rbc'){
          plot(1:num_of_codons(),newcd(),type='l',main='Read counts by codon position',xlab='Codon position',ylab='Read count')}
      })
    })
  }
  
  reactive({H5close()})

}
  

shinyApp(ui = shinyUI, server = shinyServer)