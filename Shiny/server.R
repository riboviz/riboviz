library(shiny)
library(rhdf5)
library(ggplot2)
library(reshape2)
library(cowplot)
library(plotly)

shinyServer<-function(input,output){
  path<-'./'
  df<-read.table(paste0(path,'Supp/data_unq.tsv')) #stores the year/author/database information
  names(df)=c('year','author','journal','NA','database','NA','NA','NA','dbtype','condition','geoID')
  year_vector<-as.vector(unique(df$year)) #vector containing all years
  author_vector<-as.vector(unique(df$author)) #vector containing all authors from all years
  
  rpkmpath<-'./RPKM/'
  d1<-read.delim(paste0(rpkmpath,'F8_RPKMs_modified.tsv'))
  d2<-read.delim(paste0(rpkmpath,'bgdb.tsv'))
  
  max_db <- 3 #set the maximum databases to be selected
  color3<-c(rgb(67,162,202,maxColorValue=255),rgb(254,178,76,maxColorValue=255),rgb(136,86,167,maxColorValue = 255))
  
  #---------------------------------functions
  myfun1<-function(str_year){ #returns the author selections for UI year
    as.character(unique(df[df$year==str_year,2])) 
  }
  
  myfun_errormsg<-function(){ #when an error message displays, all the graphs and legends will disappear
    output$plot11<-renderPlot('')
    output$plot22<-renderPlot('')
    output$plots<-renderPlot('')
    output$l1<-renderText('')
    output$l2<-renderText('')
    output$l3<-renderText('')
  }
  #---------------------------------
  var1<-reactive({ #returns a list of author name choices corresponding to UI year
    i <- which(year_vector==input$selectyear)
    myfun1(year_vector[i])})
  
  var2<-reactive({ #returns a list of RPF database names corresponding toUI year and UI author
    i<-which(author_vector==input$var_author)
    as.character(df[(df$year==input$selectyear)&(df$author==author_vector[i]),5])
  })
  
  output$vauthor<-renderUI({ #create the dropdown menu for author selection
    selectInput('var_author','Select the author',choices=c(var1()))})
  
  output$vdb<-renderUI({ #create the dropdown menu for RPF databse selection
    selectInput('var_db','Select the databases',choices=c(var2()))})
  
  values<-reactiveValues() 
  values$selectdb<-vector()#create a reactive value to store the UI combination of year/author/database
  a<-eventReactive(input$clickadd,{if((input$var_db!='Choose')&(input$selectyear!='Choose')&(input$var_author!='Choose')){paste(input$selectyear,input$var_author,input$var_db,sep='/')}})
  
  observe({
    values$selectdb<-rbind(isolate(values$selectdb),a())
  })
  
  observe({ #when reset button(input$clickempty) is pressed, values$selectdb is reset to empty
    if (input$clickempty==0) return()
    values$selectdb <- NULL
  })
  
  temp<-reactive({tail(unique(values$selectdb),max_db)}) #stores the last 3 user selected databases
  #temp() e.g. 2016/Weinberg/unselected_total_RNA
  temp1<-reactive({gsub('/','\t',temp())})
  #temp1() e.g. 2016 Weinberg unselected_total_RNA
  output$sdb<-renderUI(if(length(values$selectdb)>0){HTML(paste(temp1(),collapse = '<br/>'))}) #print on the webpage what these <=3 databases are
  output$dbselect_text<-renderText(paste('You may select ',3-length(temp()),' more databases'))
  
  #-----------------------------------below is to plot
  observeEvent(
    input$clickplot,{ 
      if (length(temp())==0){ #check if database selection is empty
        output$errormsg<-renderText('Error: please select at least one database')
        myfun_errormsg() 
      }
      else{
        #store all the intermediate reactive values needed for making the plot
        state<-reactiveValues()
        state$temp<-temp() #use this reactive value because otherwise the temp() can't be isolated in the renderPlot function,leading to error
        state$temp1<-temp1()
        state$geoid<-vector()
        state$inf<-vector()
        state$transtemp<-vector()
        state$attrdat<-vector()
        
        state$coord<-vector()
        state$d<-list() 
        state$dd<-list()
        state$sumd<-vector()
        state$meanreadsrbp<-vector()
        state$meanreadsrbpl<-vector()
        state$meanreadsrbc<-vector()
        state$cd<-list()
        state$num_of_codons<-vector()
        state$newcd<-vector()
        
        state$vline<-vector()
        state$dtype<-vector()
        
        #legends are above the plots on top left
        if (length(state$temp1)==1){
          output$l1<-renderText(state$temp1[1])
          output$l2<-renderText('')
          output$l3<-renderText('')
        }
        if (length(state$temp1)==2){
          output$l1<-renderText(state$temp1[1])
          output$l2<-renderText(state$temp1[2])
          output$l3<-renderText('')
        }
        if (length(state$temp1)==3){
          output$l1<-renderText(state$temp1[1])
          output$l2<-renderText(state$temp1[2])
          output$l3<-renderText(state$temp1[3])
        }
        
        if(input$txt %in% names(d2)){ #check if gene name is correct
          output$errormsg<-renderText('')     
          #variables used to plot RPKM
          igene<-which(names(d2)==input$txt) #indexes for same gene in d1 and d2 are the same
          dat_mrna<-d2[d2$Type=='mRNA',c(5,igene)]
          dat_rpf<-d2[d2$Type=='RPF',c(5,igene)]
          #------------------
          myfun2<-function(x,y,lb,fname){ #download handler
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
          
          for (i in 1:length(temp())){ #temp() e.g. 2016/Weinberg/unselected_total_RNA
            state$geoid[i]<-as.character(df[(df$year==unlist(strsplit(temp()[i],'/'))[1])&(df$author==unlist(strsplit(temp()[i],'/'))[2])&(df$database==unlist(strsplit(temp()[i],'/'))[3]),11])
            state$inf[i]<-paste0(path,temp()[i],'/',state$geoid[i],'.h5')
            state$transtemp[i]<-paste(unlist(strsplit(temp()[i],'/'))[1],unlist(strsplit(temp()[i],'/'))[2],unlist(strsplit(temp()[i],'/'))[3],sep='_')
            #transtemp e.g. 2016_Weinberg_unselected_total_RNA
            state$attrdat<-c(state$attrdat,h5readAttributes(file=state$inf[i],name=paste(input$txt,'/',state$transtemp[i],'/reads',sep=''))$reads_by_len)
            #output$test <- renderText(state$attrdat)
            
            state$coord<-dim(h5read(state$inf[1],paste(input$txt,'/',state$transtemp[1],'/reads/data',sep='')))[2] #nucleotide coordinate including left(250) and right(247) buffer, e.g.3980 for YAL001C
            state$d[[i]]<-h5read(state$inf[i],paste(input$txt,'/',state$transtemp[i],'/reads/data',sep=''))
            state$dd[[i]]<-h5read(state$inf[i],paste(input$txt,'/',state$transtemp[i],'/reads/data',sep=''))[,251:(state$coord-247)] #e.g. for YAL001C is a 36x(3980-250-247) matrix
            state$sumd<-c(state$sumd,lapply(data.frame(state$dd[[i]]),sum)) #stores reads for all length, all positions
            state$meanreadsrbp<-c(state$meanreadsrbp,sum(state$dd[[i]])/(state$coord-250-247)) #mean reads=reads total/total positions
            
            state$cd[[i]]<-rbind(state$d[[i]][14:15,(251-15):(state$coord-247-15)],state$d[[i]][16,(251-16):(state$coord-247-16)]) #3X3483, matrix to calculate the ribosome density
            state$num_of_codons<-(state$coord-247-251+1)/3 #number of codons= integer
            state$newcd<-c(state$newcd,lapply(data.frame(matrix(state$cd[[i]],nrow=nrow(state$cd[[i]])*3,ncol=state$num_of_codons)),sum)) #9*(3483/3) matrix--->1161*1 (after summing up the 9 rows in a column)
            
            state$vline[i]<-d1[d1$Year==(unlist(strsplit(temp()[i],'/'))[1]) & d1$Author==(unlist(strsplit(temp()[i],'/'))[2]) & d1$Dataset==(unlist(strsplit(temp()[i],'/'))[3]),igene]
            state$dtype[i]<-as.character(d1[d1$Year==(unlist(strsplit(temp()[i],'/'))[1]) & d1$Author==(unlist(strsplit(temp()[i],'/'))[2]) & d1$Dataset==(unlist(strsplit(temp()[i],'/'))[3]),4]) #mRNA or RPF
          }
          
          #---------------------------------- plot functions so that no redundant codes for plotting the original and zooming plots
          myfun_ggplot<-function(data,title,labx,laby){
            ggplot(data,aes(x=x,y=value,color=variable)) +
              geom_line() +
              ggtitle(title) +
              labs(x=labx,y=laby) +
              theme_classic() +
              theme(
                axis.line.x = element_line(colour = 'black', size=0.6, linetype='solid'),
                axis.line.y = element_line(colour = 'black', size=0.6, linetype='solid'),
                legend.position = 'none', #remove legend
                axis.text = element_text(size=15), #axis label size
                axis.ticks.length=unit(0.3,"cm"), #axis ticks size
                axis.title=element_text(size=15), #axis title size
                plot.title = element_text(size=18)
              )
          }
          
          myfun_ignorelegend<-function(ix){
            if(1 %in% ix){output$l1<-renderText('')} #don't display the legends for the database that has meanreads=0 b/c it's not plotted
            if(2 %in% ix){output$l2<-renderText('')}
            if(3 %in% ix){output$l3<-renderText('')}
          }
          
          myfun_rbp<-function(){
            ########################## direct reads
            ptype$x<-1 #to ensure that plott22 only appears when corresponding plot11 is updated
            myfun2(0:(state$coord-250-247-1),state$sumd,'readPosition','yeast_reads_based_on_position.csv')
            x<-0:(state$coord-250-247-1)
            y<-data.frame(matrix(as.integer(state$sumd),ncol=length(state$temp)))
            dat<-melt(data.frame(x,y),id='x')
            myfun_ggplot(dat,"Read counts by position",'Read position','Read count') +
              scale_color_manual(values=color3)
          }
          ################################### normalized reads
          myfun_rbp_norm<-function(){
            myfun2(0:(state$coord-250-247-1),state$sumd,'readPosition','yeast_reads_based_on_position.csv')
            if (all(state$meanreadsrbp==0)){
              ptype$x<-0
              output$errormsg<-renderText(paste0('All the databases you selected for gene ',input$txt, ' at length=',input$txt2, ' have meanreads that equal to 0, therefore no normalized reads can be plotted'))
              myfun_errormsg()
            }
            
            else if (any(state$meanreadsrbp==0) & !all(state$meanreadsrbp==0)){
              ptype$x<-2
              ix<-which(state$meanreadsrbp==0)
              output$errormsg<-renderText(paste0('For database(s) ',paste(state$temp1[ix],collapse = ', '),', the meanreads are equal to 0 therefore no normalized reads can be plotted.'))
              myfun_ignorelegend(ix)
              
              ixx<-which(state$meanreadsrbp!=0)
              x<-0:(state$coord-250-247-1)
              normy<-data.frame(t(t(matrix(as.integer(state$sumd),ncol=length(state$temp))[,ixx])/state$meanreadsrbp[ixx]))
              #normalized reads, use t(t(a)/k) to divide each column of matrix a (a=actual read for db selection matrix) by each element in vector k(k=mean reads for db selection vector)
              dat<-melt(data.frame(x,normy),id='x')
              myfun_ggplot(dat,"Normalized read counts by position",'Read position','Normalized reads') +
                scale_color_manual(values=color3[ixx])
            }
            
            else {
              ptype$x<-2
              x<-0:(state$coord-250-247-1)
              normy<-data.frame(t(t(matrix(as.integer(state$sumd),ncol=length(state$temp)))/state$meanreadsrbp))
              dat<-melt(data.frame(x,normy),id='x')
              myfun_ggplot(dat,"Normalized read counts by position",'Read position','Normalized reads') +
                scale_color_manual(values=color3)
            }
          }
          
          myfun_rbpl<-function(){
            ptype$x<-3
            for (i in 1:length(state$temp)){
              state$ddat[[i]]<-state$dd[[i]][isolate((as.numeric(input$txt2)-14)),]}
            myfun2(0:(state$coord-250-247-1),unlist(state$ddat),paste('readPosition_length',input$txt2,sep=''),'yeast_reads_based_on_position_location.csv')
            x<-0:(state$coord-250-247-1)
            y<-data.frame(matrix(as.integer(unlist(state$ddat)),ncol=length(state$temp)))
            dat<-melt(data.frame(x,y),id='x')
            myfun_ggplot(dat,"Read counts by position and length",'Read position','Read count') +
              scale_color_manual(values=color3)
          }
          
          myfun_rbpl_norm<-function(){
            for (i in 1:length(state$temp)){
              state$ddat[[i]]<-state$dd[[i]][isolate((as.numeric(input$txt2)-14)),]
              state$meanreadsrbpl[[i]]<-sum(state$dd[[i]][isolate((as.numeric(input$txt2)-14)),])/(state$coord-250-247)}
            myfun2(0:(state$coord-250-247-1),unlist(state$ddat),paste('readPosition_length',input$txt2,sep=''),'yeast_reads_based_on_position_location.csv')
            
            if (all(state$meanreadsrbpl==0)){ # if all selected databases have meanreads=0
              ptype$x<-0
              output$errormsg<-renderText(paste0('All the databases you selected for gene ',input$txt, ' at length=',input$txt2, ' have meanreads that equal to 0, therefore no normalized reads can be plotted'))
              myfun_errormsg()
            }
            
            else if (any(state$meanreadsrbpl==0) & !all(state$meanreadsrbpl==0)){ # if some of the selected databases have meanreads=0
              ptype$x<-4
              ix<-which(state$meanreadsrbpl==0)
              output$errormsg<-renderText(paste0('For database(s) ',paste(state$temp1[ix],collapse = ', '),', the meanreads are equal to 0 therefore no normalized reads can be plotted.'))
              myfun_ignorelegend(ix)
              
              ixx<-which(state$meanreadsrbpl!=0)
              x<-0:(state$coord-250-247-1)
              normy<-data.frame(t(t(matrix(as.integer(unlist(state$ddat)),ncol=length(state$temp))[,ixx])/state$meanreadsrbpl[ixx])) # matrix of row=number of reads, col=number of selected db
              dat<-melt(data.frame(x,normy),id='x')
              myfun_ggplot(dat,"Normalized read counts by position and length",'Read position','Normalized reads') +
                scale_color_manual(values=color3[ixx]) # the color assignment should be same as when normalized reads for all databases can be plotted
            }
            
            else { #if none of the selected databases have meanreads =0
              ptype$x<-4
              x<-0:(state$coord-250-247-1)
              normy<-data.frame(t(t(matrix(as.integer(unlist(state$ddat)),ncol=length(state$temp)))/state$meanreadsrbpl)) # matrix of row=number of reads, col=number of selected db
              dat<-melt(data.frame(x,normy),id='x')
              myfun_ggplot(dat,"Normalized read counts by position and length",'Read position','Normalized reads') +
                scale_color_manual(values=color3)
            }
          }
          
          myfun_rbc<-function(){
            ptype$x<-5
            myfun2(1:state$num_of_codons,unlist(state$newcd),'codonPosition','yeast_reads_based_on_codons.csv')
            x<-1:state$num_of_codons
            y<-data.frame(matrix(as.numeric(unlist(state$newcd)),ncol=length(state$temp)))
            maxy<-max(y)
            names(y)<-state$temp1
            dat<-melt(data.frame(x,y),id='x')
            myfun_ggplot(dat,"Read counts by codon position",'Codon position','Read count') +
              scale_color_manual(values=color3)
          }
          
          myfun_rbc_norm<-function(){
            myfun2(1:state$num_of_codons,unlist(state$newcd),'codonPosition','yeast_reads_based_on_codons.csv')
            state$meanreadsrbc<-apply(matrix(as.numeric(unlist(state$newcd)),ncol=length(state$temp)),2,sum)/(state$num_of_codons)
            
            if (all(state$meanreadsrbc==0)){
              ptype$x<-0
              output$errormsg<-renderText(paste0('All the databases you selected for gene ',input$txt, ' at length=',input$txt2, ' have meanreads that equal to 0, therefore no normalized reads can be plotted'))
              myfun_errormsg()
            }
            
            else if (any(state$meanreadsrbc==0) & !all(state$meanreadsrbc==0)){
              ptype$x<-6
              ix<-which(state$meanreadsrbc==0)
              output$errormsg<-renderText(paste0('For database(s) ',paste(state$temp1[ix],collapse = ', '),', the meanreads are equal to 0 therefore no normalized reads can be plotted.'))
              myfun_ignorelegend(ix)
              
              ixx<-which(state$meanreadsrbc!=0)
              x<-1:state$num_of_codons
              normy<-data.frame(t(t(matrix(as.numeric(unlist(state$newcd)),ncol=length(state$temp))[,ixx])/state$meanreadsrbc[ixx]))
              dat<-melt(data.frame(x,normy),id='x')
              myfun_ggplot(dat,"Normalized read counts by codon position",'Codon position','Normalized reads') +
                scale_color_manual(values=color3[ixx])
            }
            
            else {
              ptype$x<-6
              x<-1:state$num_of_codons
              normy<-data.frame(t(t(matrix(as.numeric(unlist(state$newcd)),ncol=length(state$temp)))/state$meanreadsrbc))
              dat<-melt(data.frame(x,normy),id='x')
              myfun_ggplot(dat,"Normalized read counts by codon position",'Codon position','Normalized reads') +
                scale_color_manual(values=color3)
            }
          }
          
          myfun_rpkm<-function(databasetype){
            ii<-which(state$dtype==databasetype)
            vline<-state$vline[ii] #vertical line
            #lbl<-state$temp1[ii] #label
            
            if(databasetype=='mRNA'){dat<-dat_mrna}
            else if(databasetype=='RPF'){dat<-dat_rpf}
            
            p<-ggplot(dat, aes_string(x=input$txt)) +
              geom_histogram(bins = input$bins, position='dodge',alpha=0.5,color=rgb(189,189,189,maxColorValue=255),fill=rgb(189,189,189,maxColorValue=255)) +
              ggtitle(databasetype) +
              theme_classic() +
              theme(
                axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.text = element_text(size=15), #axis label size
                axis.ticks.length=unit(0.3,"cm"), #axis ticks size
                axis.title=element_text(size=15), #axis title size
                plot.title = element_text(size=18)) +
              geom_vline(xintercept = vline,color=color3[ii],size=1) 
            ggplotly(p)
            
          }
          #---------------------------------- plot11
          state$ddat  <-list()
          ranges <- reactiveValues(x = NULL)
          ptype <-reactiveValues(x=0) #this is used so that the zoomed plot(plot22) doesn' show up until plot11 of corresponding type shows first
          
          if(input$select=='rpkm'){
            output$plot11<-renderPlotly({#width = 800,height=800,{ #isolate in the renderPlot() so that it will only update when the input$clickplot button is fired
              input$clickplot
              isolate({
                p1<-myfun_rpkm('mRNA')
                p2<-myfun_rpkm('RPF')
                subplot(p1, p2)
                #plot_grid(p1, p2, ncol = 1, nrow = 2) #to display the two RPKM plots as one plot, function of library(cowplot)
              })
            }
            )
          }
          
          else {
            output$plot11<-renderPlotly({#width = 1000,height=500,{
              input$clickplot
              isolate({
                if(input$select=='rbl'){
                  myfun2(15:50,state$attrdat,'readLength','yeast_reads_based_on_location.csv')
                  x<-15:50
                  y<-data.frame(prop.table(matrix(as.integer(state$attrdat),ncol=length(state$temp)),2))
                  maxy<-max(y)
                  names(y)<-state$temp1
                  dat<- melt(data.frame(x, y),id='x')
                  p<-ggplot(dat, aes(x = x, y=value, fill = variable)) +
                    geom_bar(stat='identity',position = 'dodge') +
                    ggtitle("Reads by length") +
                    labs(x='Read length',y='Frequency of read lengths') +
                    scale_x_continuous(breaks=seq(15,50,by=1)) +
                    theme_classic() +
                    theme(
                      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                      legend.position = 'none', #remove legend
                      axis.text = element_text(size=15), #axis label size
                      axis.ticks.length=unit(0.3,"cm"), #axis ticks size
                      axis.title=element_text(size=15), #axis title size
                      plot.title = element_text(size=18)
                    ) +
                    scale_fill_manual(values=color3)
                  ggplotly(p)
                }
                
                else if(input$select=='rbp' & input$radioreads == '1'){
                  myfun_rbp()
                }
                
                else if(input$select=='rbp' & input$radioreads == '2'){
                  myfun_rbp_norm()
                }
                
                else if(input$select=='rbpl' & input$radioreads == '1'){
                  if (input$txt2 %in% c(15:50)){ #check if input$txt2 is valid
                    output$errormsg<-renderText('')
                    myfun_rbpl()
                  }
                  else{
                    output$errormsg<-renderText('Error: please enter a length between 15-50')
                    myfun_errormsg()
                  }
                }
                
                else if(input$select=='rbpl' & input$radioreads == '2'){
                  if (input$txt2 %in% c(15:50)){ #check if input$txt2 is valid
                    output$errormsg<-renderText('')
                    myfun_rbpl_norm()
                  }
                  else{
                    output$errormsg<-renderText('Error: please enter a length between 15-50')
                    myfun_errormsg()
                  }
                }
                
                else if(input$select=='rbc' & input$radioreads == '1'){
                  myfun_rbc()
                }
                
                else if(input$select=='rbc' & input$radioreads == '2'){
                  myfun_rbc_norm()
                }
              })
            }
            )
          }
          
          #---------------------------------------------------------------------plot22: zooming plots
          output$plot22<-renderPlot(width = 1000,height=500,{
            if(input$select=='rbp' & ptype$x==1 & input$radioreads=='1'){
              myfun_rbp() +
                coord_cartesian(xlim = ranges$x) +
                ggtitle('Select a region to zoom')
            }
            
            else if(input$select=='rbp' & ptype$x==2 & input$radioreads=='2'){
              myfun_rbp_norm() +
                coord_cartesian(xlim = ranges$x) +
                ggtitle('Select a region to zoom')
            }
            
            else if(input$select=='rbpl' & ptype$x==3 & input$radioreads=='1'){
              myfun_rbpl() +
                coord_cartesian(xlim = ranges$x) +
                ggtitle('Select a region to zoom')
            }
            
            else if(input$select=='rbpl' & ptype$x==4 & input$radioreads=='2'){
              myfun_rbpl_norm() +
                coord_cartesian(xlim = ranges$x) +
                ggtitle('Select a region to zoom')
            }
            
            else if(input$select=='rbc' & ptype$x==5 & input$radioreads=='1'){
              myfun_rbc() +
                coord_cartesian(xlim = ranges$x) +
                ggtitle('Select a region to zoom')
            }
            
            else if(input$select=='rbc' & ptype$x==6 & input$radioreads=='2'){
              myfun_rbc_norm() +
                coord_cartesian(xlim = ranges$x) +
                ggtitle('Select a region to zoom')
            }
          })
          #---------------------
          # Check if there's a brush on the plot.
          # If so, zoom to the brush bounds; if not, reset the zoom.
          observe({
            brush <- input$plot11_brush
            if (!is.null(brush)) {
              ranges$x <- c(brush$xmin, brush$xmax)
            } else {
              ranges$x <- NULL
            }
          })
        }
        #--------------------
        else{
          output$errormsg<-renderText('Error: the gene you entered does not exist in the database, please re-enter')
          myfun_errormsg()
        }
        #--------------------
      }
    }
  )
  reactive({H5close()})
}
