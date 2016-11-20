library(shiny)
library(rhdf5)

#-------------------------------------------- part3: shinyServer
shinyServer<-function(input,output){
  path<-'mouse/'
  datapath<-'mouse/'
  df<-read.table(paste0(path,'samplePath.csv')) #stores the tissue/*/sample information
  df2<-read.delim(paste0(path,'mouse_genes.txt'))
  names(df)=c('tissue','rpf','sample','sampletype')
  tissue_vector<-as.vector(unique(df$tissue)) #vector containing all years
  sample_vector<-as.vector(unique(df$sample)) #vector containing all authors from all years
  max_db <- 5 #set the maximum databases to be selected 
  
  #-------------------------------------------- part1:self defined functions
  myfun1<-function(str_tissue){ #returns the tissue selections for UI year
    as.character(unique(df[df$tissue==str_tissue,3])) 
  }
  
  myfun2<-function(str_tissue,str_samplename){ #returns the database selections for UI tissue/sample
    as.character(df[(df$tissue==str_tissue)&(df$sample==str_samplename),4])
  }
  
  #-------------------------------------------- 
  var1<-reactive({ #returns a list of sample name choices corresponding to UI tissue
    i <- which(tissue_vector==input$selecttissue)
    myfun1(tissue_vector[i])})
  
  var2<-reactive({ #returns a list of sample type corresponding to UI tissue and UI sample
    i<-which(sample_vector==input$var_sample)
    myfun2(input$selecttissue,sample_vector[i])})
  
  output$vsample<-renderUI({ #create the dropdown menu for sample selection
    selectInput('var_sample','Select the sample name',choices=c('Choose',var1()))})
  
  output$vsampletype<-renderUI({ #create the dropdown menu for sample type selection
    selectInput('var_sampletype','Select the sample type',choices=c('Choose',var2()))})
  
  values<-reactiveValues()
  values$selectdb<-vector()#create a reactive value to store the UI combination of tissue/sample/sampletype
  values$selectdb1<-vector()# create a reactive value to store the UI combination of tissue/sample/sampletype to print on the screen
  a<-eventReactive(input$clickadd,{if((input$selecttissue!='Choose')&(input$var_sample!='Choose')&(input$var_sampletype!='Choose')){paste(input$selecttissue,input$var_sample,input$var_sampletype,sep='/')}})
  b<-eventReactive(input$clickadd,{if((input$selecttissue!='Choose')&(input$var_sample!='Choose')&(input$var_sampletype!='Choose')){
    ifelse (tail(unlist(strsplit(input$var_sampletype,'_')),1)=='lncrna.h5',
            paste(input$selecttissue,input$var_sample,'Long non-coding RNA',sep='\t'), 
            paste(input$selecttissue,input$var_sample,'ORF',sep='\t'))
  }})
  
  observe({
    values$selectdb<-rbind(isolate(values$selectdb),a())
    values$selectdb1<-rbind(isolate(values$selectdb1),b())
  })
  
  observe({ #when reset button(input$clickempty) is pressed, values$selectdb is reset to empty
    if (input$clickempty==0) return()
    values$selectdb<-NULL
    values$selectdb1<-NULL
  })
  
  temp<-reactive({tail(unique(values$selectdb),max_db)}) #stores the last 5 user selected databases for downstream plotting
  temp1<-reactive({tail(unique(values$selectdb1),max_db)}) #stores the last 5 user selected databases for printing
  output$sdb<-renderUI(if(length(values$selectdb1)>0){HTML(paste(temp1(),collapse = '<br/>'))}) #print on the webpage what these <=5 databases are
  
  #-----------------------------------below is to plot
  observeEvent(
    input$clickplot,{ #store all the intermediate reactive values needed for making the plot
      state<-reactiveValues()
      state$temp<-temp() #use this reactive value because otherwise the temp() can't be isolated in the renderPlot function,leading to error
      state$rpfid<-vector()
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
        state$rpfid[i]<-as.character(df[(df$tissue==unlist(strsplit(temp()[i],'/'))[1])&(df$sample==unlist(strsplit(temp()[i],'/'))[2])&(df$sampletype==unlist(strsplit(temp()[i],'/'))[3]),2])
        ## "Ribosome_protected_fragment_data"
        state$inf[i]<-paste(datapath,unlist(strsplit(temp()[i],'/'))[1],'/',state$rpfid[i],'/',unlist(strsplit(temp()[i],'/'))[2],'/',unlist(strsplit(temp()[i],'/'))[3],sep='')
        ## "/data/riboseq/mouse/Muscle/Ribosome_protected_fragment_data/Sample_M-RP-C/Sample_M-RP-C.h5"
        state$transtemp[i]<-paste(unlist(strsplit(temp()[i],'/'))[1],state$rpfid[i],unlist(strsplit(temp()[i],'/'))[2],sep='_')
        ## "Muscle_Ribosome_protected_fragment_data_Sample_M-RP-C"
        state$attrdat<-c(state$attrdat,h5readAttributes(file=state$inf[i],name=paste(input$txt,'/',state$transtemp[i],'/reads',sep=''))$reads_by_len)
        ## name="ENSMUST00000000003.13/Muscle_Ribosome_protected_fragment_data_Sample_M-RP-C/reads"
        
        state$coord<-dim(h5read(state$inf[1],paste(input$txt,'/',state$transtemp[1],'/reads/data',sep='')))[2] #nucleotide coordinate e.g.1399 for ENSMUST00000000003.13
        state$d[[i]]<-h5read(state$inf[i],paste(input$txt,'/',state$transtemp[i],'/reads/data',sep='')) #e.g. for ENSMUST00000000003.13 is a 36x1399 matrix
        state$sumd<-c(state$sumd,lapply(data.frame(state$d[[i]]),sum)) #stores reads for all length, all positions
        
        if (grepl('lncrna',unlist(strsplit(state$temp[i],'/'))[3])==FALSE) #for lncrna samples, to avoid errors when the click to plot button is pressed(b/c once it's clicked these values will be stored)
        {state$cd[[i]]<-rbind(state$d[[i]][15:16,(250+df2[df2$transcript_id==input$txt,2]+1-15):(state$coord-247-df2[df2$transcript_id==input$txt,4]-15)],
                              state$d[[i]][17:19,(250+df2[df2$transcript_id==input$txt,2]+1-16):(state$coord-247-df2[df2$transcript_id==input$txt,4]-16)],
                              state$d[[i]][20:21,(250+df2[df2$transcript_id==input$txt,2]+1-17):(state$coord-247-df2[df2$transcript_id==input$txt,4]-17)]) #7X852 for ENSMUST00000056406.6, matrix to calculate the ribosome density
        state$num_of_codons<-(df2[df2$transcript_id==input$txt,3])/3 #number of codons= integer
        state$newcd<-c(state$newcd,lapply(data.frame(matrix(state$cd[[i]],nrow=nrow(state$cd[[i]])*3,ncol=state$num_of_codons)),sum))}
      }
      
      state$ddat<-list()
      output$plot1<-renderPlot({ #isolate in the renderPlot() so that it will only update when the input$clickplot button is fired
        input$clickplot
        isolate({
          if(input$select=='rbl'){
            matplot(15:50,matrix(state$attrdat,ncol=length(state$temp)),type='l',main='Read counts by Length',xlab='Read length',ylab='Read count',col =1:length(state$temp),lty=1,lwd=2)
            legend("topleft", legend = state$temp, col=1:length(state$temp),lty=1,lwd=2)
            myfun3(15:50,state$attrdat,'readLength','mouse_reads_based_on_location.csv')}
          
          else if(input$select=='rbp'){
            matplot(1:state$coord,matrix(state$sumd,ncol=length(state$temp)),type='l',main='Read counts by Position',xlab='Read position',ylab='Read count',col =1:length(state$temp),lty=1,lwd=2)
            legend("topleft", legend = state$temp, col=1:length(state$temp),lty=1,lwd=2)
            myfun3(1:state$coord,state$sumd,'readPosition','mouse_reads_based_on_position.csv')}
          
          else if(input$select=='rbpl'){
            for (i in 1:length(state$temp)){
              #print(i) always print two sets of i, strange!!!
              state$ddat[[i]]<-state$d[[i]][isolate((as.numeric(input$txt2)-14)),]}
            matplot(1:state$coord,matrix(unlist(state$ddat),ncol=length(state$temp)),type='l',main='Read counts by Position and length',xlab='Read position',ylab='Read count',lty=1,lwd=2,col =1:length(state$temp))
            legend("topleft", legend = state$temp, col=1:length(state$temp),lty=1,lwd=2)
            myfun3(1:state$coord,unlist(state$ddat),paste('readPosition_length',input$txt2,sep=''),'mouse_reads_based_on_position_location.csv')}
          
          else if(input$select=='rbc'){
            matplot(1:state$num_of_codons,matrix(unlist(state$newcd),ncol=length(state$temp)),type='l',main='Read counts by codon position',xlab='Codon position',ylab='Read count',lty=1,lwd=2,col =1:length(state$temp))
            legend("topleft", legend = state$temp, col=1:length(state$temp),lty=1,lwd=2)
            myfun3(1:state$num_of_codons,unlist(state$newcd),'codonPosition','mouse_reads_based_on_codons.csv')}
        })
      })
    })
  
  reactive({H5close()})
  
}
