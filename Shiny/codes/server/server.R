library(shiny)
library(rhdf5)
#For future database update, you need an updated file listing all year,author and database name like 'data_unq.tsv',

#-------------------------------------------- part3: shinyServer
shinyServer<-function(input,output){
	df<-read.table('Supp/data_unq.tsv') #stores the year/author/database information
	names(df)=c('year','author','journal','NA','database','NA','NA','NA','dbtype','condition','geoID')
	year_vector<-as.vector(unique(df$year)) #vector containing all years
	author_vector<-as.vector(unique(df$author)) #vector containing all authors from all years
	max_db <- 5 #set the maximum databases to be selected 

	#-------------------------------------------- part1:self defined functions
	myfun1<-function(str_year){ #returns the author selections for UI year
	  as.character(unique(df[df$year==str_year,2])) 
	}

	myfun2<-function(str_author,str_year,str_dbtype){ #returns the database selections for UI year/author/ for either mRNA or RPF databases
	  as.character(df[(df$year==str_year)&(df$author==str_author)&(df$dbtype==str_dbtype),5])
	}

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
  output$sdb<-renderPrint(if(length(values$selectdb)>0){print(temp())}) #print on the webpage what these <=5 databases are
  
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
          legend("topleft", legend = state$temp, col=1:length(state$temp),pch=1)}

          
        else if(input$select=='rbp'){
          matplot(1:state$coord,matrix(state$sumd,ncol=length(state$temp)),type='l',main='Read counts by Position',xlab='Read position',ylab='Read count',pch=2,col =1:length(state$temp))
          legend("topleft", legend = state$temp, col=1:length(state$temp),pch=1)}

        else if(input$select=='rbpl'){
          for (i in 1:length(state$temp)){
            #print(i) always print two sets of i, strange!!!
            state$ddat[[i]]<-state$d[[i]][isolate((as.numeric(input$txt2)-14)),]}
            matplot(1:state$coord,matrix(unlist(state$ddat),ncol=length(state$temp)),type='l',main='Read counts by Position and length',xlab='Read position',ylab='Read count',pch=2,col =1:length(state$temp))
            legend("topleft", legend = state$temp, col=1:length(state$temp),pch=1)}
        
        else if(input$select=='rbc'){
          matplot(1:state$num_of_codons,matrix(unlist(state$newcd),ncol=length(state$temp)),type='l',main='Read counts by codon position',xlab='Codon position',ylab='Read count',pch=2,col =1:length(state$temp))
          legend("topleft", legend = state$temp, col=1:length(state$temp),pch=1)}
        })
      })
    })

  reactive({H5close()})
}
  
#shinyApp(ui = shinyUI, server = shinyServer)
