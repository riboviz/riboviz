library(shiny)
library(rhdf5)
#For future database update, you need an updated file listing all year,author and database name like 'data_unq.tsv',

setwd('/data/riboseq/')
df<-read.table('Supp/data_unq.tsv') #stores the year/author/database information
names(df)=c('year','author','journal','NA','database','NA','NA','NA','dbtype','condition','geoID')
year_vector<-as.vector(unique(df$year)) #vector containing all years
author_vector<-as.vector(unique(df$author)) #vector containing all authors from all years
max_db <- 5 #set the maximum databases to be selected 

#-------------------------------------------- part1:functions
myfun1<-function(str_year){ #returns the author selections for UI year
  as.character(unique(df[df$year==str_year,2])) 
}

myfun2<-function(str_author,str_year,str_dbtype){ #returns the database selections for UI year/author/ for either mRNA or RPF databases
  as.character(df[(df$year==str_year)&(df$author==str_author)&(df$dbtype==str_dbtype),5])
}
#-------------------------------------------- part2:define webpage components and layout
shinyUI<-fluidPage(
  titlePanel('Explore the RPF and mRNA databases'),
  sidebarLayout(
    sidebarPanel(
      
      #select input with the list of datasets
      selectInput(inputId = 'selectyear',label='Select the year of publication',choices=c('Choose',year_vector)),
      
      uiOutput('vauthor'), #3 dynamic drop down menus
      uiOutput('vrpf'),
      uiOutput('vmrna'),
      
      actionButton(inputId = 'clickadd',label='Add to selected database (max 5)'),
      br(),
      br(),
      actionButton(inputId = 'clickempty',label='Reset database selections'),
      br(),
      br(),
#------------above is to select databases,below is to select plot type     
      textInput(inputId = 'txt',label='Write the gene name of interest',placeholder = 'e.g. YAL001C'),
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
  
  values<-reactiveValues() 
  values$selectdb<-vector()#create a reactive value to store the UI combination of year/author/database
  a<-eventReactive(input$clickadd,{if((input$var_rpf!='Choose')&(input$selectyear!='Choose')&(input$var_author!='Choose')){paste(input$selectyear,input$var_author,input$var_rpf,sep='/')}})
  b<-eventReactive(input$clickadd,{if((input$var_mrna!='Choose')&(input$selectyear!='Choose')&(input$var_author!='Choose')){paste(input$selectyear,input$var_author,input$var_mrna,sep='/')}})
  
  observe({
    values$selectdb<-rbind(isolate(values$selectdb),a(),b())
  })
  
  observe({ #when reset button(input$clickempty) is pressed, values$selectdb is reset to empty
    if (input$clickempty==0) return()
    values$selectdb<-NULL
  })
  
  temp<-reactive({tail(unique(values$selectdb),max_db)}) #stores the last 5 user selected databases
  temp1<-reactive({gsub('/','\t',temp())})
  output$sdb<-renderUI(if(length(values$selectdb)>0){HTML(paste(temp1(),collapse = '<br/>'))}) #print on the webpage what these <=5 databases are
  
#-----------------------------------below is to plot
   observeEvent(
    input$clickplot,{ #store all the intermediate reactive values needed for making the plot
       state<-reactiveValues()
       state$temp<-temp() #use this reactive value because otherwise the temp() can't be isolated in the renderPlot function,leading to error
       state$geoid<-vector()
       state$inf<-vector()
       state$transtemp<-vector()
       state$attrdat<-vector()
       
       state$coord<-vector()
       state$d<-list()
       state$sumd<-vector()
       state$cd<-list()
       state$num_of_codons<-vector()
       state$newcd<-vector()
       
      myfun3<-function(x,y,lb,fname){ #download handler
       #x and y are the ones used by plotting, lb is the 1st column name for the *.csv file, fname is the file name
       data<-reactive({
         d<-matrix(c(x,y),ncol=(length(state$temp))+1)
         colnames(d) <- c(lb,state$temp)
         d
       })
       output$down <- downloadHandler(
         filename = function() {
           paste(fname, sep='')
         },
         content = function(file) {
           write.csv(data(), file)})
     }
     
      for (i in 1:length(temp())){
        state$geoid[i]<-as.character(df[(df$year==unlist(strsplit(temp()[i],'/'))[1])&(df$author==unlist(strsplit(temp()[i],'/'))[2])&(df$database==unlist(strsplit(temp()[i],'/'))[3]),11])
        state$inf[i]<-paste(temp()[i],'/',state$geoid[i],'.h5',sep='')
        state$transtemp[i]<-paste(unlist(strsplit(temp()[i],'/'))[1],unlist(strsplit(temp()[i],'/'))[2],unlist(strsplit(temp()[i],'/'))[3],sep='_')
        state$attrdat<-c(state$attrdat,h5readAttributes(file=state$inf[i],name=paste(input$txt,'/',state$transtemp[i],'/reads',sep=''))$reads_by_len)
        
        state$coord<-dim(h5read(state$inf[1],paste(input$txt,'/',state$transtemp[1],'/reads/data',sep='')))[2] #nucleotide coordinate e.g.3980 for YAL001C
        state$d[[i]]<-h5read(state$inf[i],paste(input$txt,'/',state$transtemp[i],'/reads/data',sep='')) #e.g. for YAL001C is a 36x3980 matrix
        state$sumd<-c(state$sumd,lapply(data.frame(state$d[[i]]),sum)) #stores reads for all length, all positions
        
        state$cd[[i]]<-rbind(state$d[[i]][14:15,(251-15):(state$coord-247-15)],state$d[[i]][16,(251-16):(state$coord-247-16)]) #3X3483, matrix to calculate the ribosome density
        state$num_of_codons<-(state$coord-247-251+1)/3 #number of codons= integer
        state$newcd<-c(state$newcd,lapply(data.frame(matrix(state$cd[[i]],nrow=nrow(state$cd[[i]])*3,ncol=state$num_of_codons)),sum))
        }
       
      state$ddat<-list()
      output$plot1<-renderPlot({ #isolate in the renderPlot() so that it will only update when the input$clickplot button is fired
        input$clickplot
        isolate({
        if(input$select=='rbl'){
          matplot(15:50,matrix(state$attrdat,ncol=length(state$temp)),type='l',main='Read counts by Length',xlab='Read length',ylab='Read count',pch=2,col =1:length(state$temp))
          legend("topleft", legend = state$temp, col=1:length(state$temp),pch=1)
          myfun3(15:50,state$attrdat,'readLength','yeast_reads_based_on_location.csv')}
        
        else if(input$select=='rbp'){
          matplot(1:state$coord,matrix(state$sumd,ncol=length(state$temp)),type='l',main='Read counts by Position',xlab='Read position',ylab='Read count',pch=2,col =1:length(state$temp))
          legend("topleft", legend = state$temp, col=1:length(state$temp),pch=1)
          myfun3(1:state$coord,state$sumd,'readPosition','yeast_reads_based_on_position.csv')}
        
        else if(input$select=='rbpl'){
          for (i in 1:length(state$temp)){
            #print(i) always print two sets of i, strange!!!
            state$ddat[[i]]<-state$d[[i]][isolate((as.numeric(input$txt2)-14)),]}
            matplot(1:state$coord,matrix(unlist(state$ddat),ncol=length(state$temp)),type='l',main='Read counts by Position and length',xlab='Read position',ylab='Read count',pch=2,col =1:length(state$temp))
            legend("topleft", legend = state$temp, col=1:length(state$temp),pch=1)
            myfun3(1:state$coord,unlist(state$ddat),paste('readPosition_length',input$txt2,sep=''),'yeast_reads_based_on_position_location.csv')}
        
        else if(input$select=='rbc'){
          matplot(1:state$num_of_codons,matrix(unlist(state$newcd),ncol=length(state$temp)),type='l',main='Read counts by codon position',xlab='Codon position',ylab='Read count',pch=2,col =1:length(state$temp))
          legend("topleft", legend = state$temp, col=1:length(state$temp),pch=1)
          myfun3(1:state$num_of_codons,unlist(state$newcd),'codonPosition','yeast_reads_based_on_codons.csv')}
        })
      })
    })

  reactive({H5close()})
}
  
shinyApp(ui = shinyUI, server = shinyServer)