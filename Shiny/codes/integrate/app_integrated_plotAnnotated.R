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
  titlePanel('Explore the RPF and mRNA databases'),

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
          textInput(inputId = 'txt2',label='Write the read length of interest',placeholder='ranging from 15-50')),
        conditionalPanel(
          condition = "input.select == 'rpkm'",
          sliderInput(inputId = 'bins',label='Choose bins for the histogram',min=5,max=30, value = 10)),
        actionButton(inputId = 'clickplot',label='Click to plot'),
        downloadButton(outputId = 'down',label='Download data')
        )
      ),
    column(8,
      textOutput('errormsg'),
      tags$head(tags$style("#errormsg{color: red}")), 
      
      plotOutput('plot11',
                 dblclick = "plot11_dblclick",
                 brush = brushOpts(
                   id = "plot11_brush",
                   resetOnNew = TRUE),
                 width = 1000,height=500),
      
      plotOutput('plot22',width = 1000,height=500),
      
      uiOutput("plots",width = 1000,height=500)
    )
      )
  )
#------------------------------------------------------------------------shinyServer
shinyServer<-function(input,output){
  path<-'/data/riboseq/'
  df<-read.table(paste0(path,'Supp/data_unq.tsv')) #stores the year/author/database information
  names(df)=c('year','author','journal','NA','database','NA','NA','NA','dbtype','condition','geoID')
  year_vector<-as.vector(unique(df$year)) #vector containing all years
  author_vector<-as.vector(unique(df$author)) #vector containing all authors from all years
  
  rpkmpath<-'/home/txing/shinyProject/data/RPKM/'
  d1<-read.delim(paste0(rpkmpath,'F8_RPKMs_modified.tsv'))
  d2<-read.delim(paste0(rpkmpath,'bgdb.tsv'))

  max_db <- 3 #set the maximum databases to be selected
  color3<-c(rgb(67,162,202,maxColorValue=255),rgb(254,178,76,maxColorValue=255),rgb(136,86,167,maxColorValue = 255))
  
  #---------------------------------functions
  myfun1<-function(str_year){ #returns the author selections for UI year
    as.character(unique(df[df$year==str_year,2])) 
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
    values$selectdb<-NULL
  })
  
  temp<-reactive({tail(unique(values$selectdb),max_db)}) #stores the last 3 user selected databases
  temp1<-reactive({gsub('/','\t',temp())})
  output$sdb<-renderUI(if(length(values$selectdb)>0){HTML(paste(temp1(),collapse = '<br/>'))}) #print on the webpage what these <=3 databases are
  output$dbselect_text<-renderText(paste('You may select ',3-length(temp()),' more databases'))
  
  #-----------------------------------below is to plot
  observeEvent(
    input$clickplot,{ 
      if (length(temp())==0){ #check if database selection is empty
        output$errormsg<-renderText('Error: please select at least one database')
        output$plot11<-renderPlot('')
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
        state$cd<-list()
        state$num_of_codons<-vector()
        state$newcd<-vector()
        
        state$vline<-vector()
        state$dtype<-vector()
    
  
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
  
          for (i in 1:length(temp())){
            state$geoid[i]<-as.character(df[(df$year==unlist(strsplit(temp()[i],'/'))[1])&(df$author==unlist(strsplit(temp()[i],'/'))[2])&(df$database==unlist(strsplit(temp()[i],'/'))[3]),11])
            state$inf[i]<-paste0(path,temp()[i],'/',state$geoid[i],'.h5')
            state$transtemp[i]<-paste(unlist(strsplit(temp()[i],'/'))[1],unlist(strsplit(temp()[i],'/'))[2],unlist(strsplit(temp()[i],'/'))[3],sep='_')
            state$attrdat<-c(state$attrdat,h5readAttributes(file=state$inf[i],name=paste(input$txt,'/',state$transtemp[i],'/reads',sep=''))$reads_by_len)
  
            state$coord<-dim(h5read(state$inf[1],paste(input$txt,'/',state$transtemp[1],'/reads/data',sep='')))[2] #nucleotide coordinate including left(250) and right(247) buffer, e.g.3980 for YAL001C
            state$d[[i]]<-h5read(state$inf[i],paste(input$txt,'/',state$transtemp[i],'/reads/data',sep=''))
            state$dd[[i]]<-h5read(state$inf[i],paste(input$txt,'/',state$transtemp[i],'/reads/data',sep=''))[,251:(state$coord-247)] #e.g. for YAL001C is a 36x(3980-250-247) matrix
            state$sumd<-c(state$sumd,lapply(data.frame(state$dd[[i]]),sum)) #stores reads for all length, all positions
  
            state$cd[[i]]<-rbind(state$d[[i]][14:15,(251-15):(state$coord-247-15)],state$d[[i]][16,(251-16):(state$coord-247-16)]) #3X3483, matrix to calculate the ribosome density
            state$num_of_codons<-(state$coord-247-251+1)/3 #number of codons= integer
            state$newcd<-c(state$newcd,lapply(data.frame(matrix(state$cd[[i]],nrow=nrow(state$cd[[i]])*3,ncol=state$num_of_codons)),sum))
  
            state$vline[i]<-d1[d1$Year==(unlist(strsplit(temp()[i],'/'))[1]) & d1$Author==(unlist(strsplit(temp()[i],'/'))[2]) & d1$Dataset==(unlist(strsplit(temp()[i],'/'))[3]),igene]
            state$dtype[i]<-as.character(d1[d1$Year==(unlist(strsplit(temp()[i],'/'))[1]) & d1$Author==(unlist(strsplit(temp()[i],'/'))[2]) & d1$Dataset==(unlist(strsplit(temp()[i],'/'))[3]),4]) #mRNA or RPF
          }
  
          #---------------------------------- plot functions so that no redundant codes for plotting the original and zooming plots
          myfun_rbp<-function(){
            ptype$x<-1
            myfun2(1:state$coord,state$sumd,'readPosition','yeast_reads_based_on_position.csv')
            x<-0:(state$coord-250-247-1)
            y<-data.frame(matrix(as.integer(state$sumd),ncol=length(state$temp)))
            maxy<-max(y)
            names(y)<-state$temp1
            dat<-melt(data.frame(x,y),id='x')
            ggplot(dat,aes(x=x,y=value,color=variable)) +
              geom_line() +
              ggtitle("Read counts by position") +
              labs(x='Read position',y='Read count') +
              theme_classic() +
              theme(
                axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                legend.position = 'none' #remove legend
                ) +
              annotate("text", x = 3/5*(max(x)-min(x)), y = c(maxy-2*(maxy/13),maxy-maxy/13,maxy)[1:length(state$temp1)], label = state$temp, colour = color3[1:length(state$temp1)],hjust = 0, size=5) +
              scale_color_manual(values=color3)
          }
  
          myfun_rbpl<-function(){
            ptype$x<-2
            for (i in 1:length(state$temp)){
              state$ddat[[i]]<-state$dd[[i]][isolate((as.numeric(input$txt2)-14)),]}
            myfun2(1:state$coord,unlist(state$ddat),paste('readPosition_length',input$txt2,sep=''),'yeast_reads_based_on_position_location.csv')
            x<-0:(state$coord-250-247-1)
            y<-data.frame(matrix(as.integer(unlist(state$ddat)),ncol=length(state$temp)))
            maxy<-max(y)
            miny<-min(y)
            names(y)<-state$temp1
            dat<-melt(data.frame(x,y),id='x')
            if (maxy==0 & miny==0){
              ggplot(dat,aes(x=x,y=value,color=variable)) +
                geom_line() +
                ggtitle("Read counts by position and length") +
                labs(x='Read position',y='Read count') +
                theme_classic() +
                theme(
                  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                  legend.position = 'none' #remove legend
                ) +
                annotate("text", x = 3/5*(max(x)-min(x)), y = c(0.2,0.18,0.16)[1:length(state$temp1)], label = state$temp, colour = color3[1:length(state$temp1)],hjust = 0, size=5) +
                ylim(0,0.2) +
                scale_color_manual(values=color3)
            }
            else{
            ggplot(dat,aes(x=x,y=value,color=variable)) +
              geom_line() +
              ggtitle("Read counts by position and length") +
              labs(x='Read position',y='Read count') +
              theme_classic() +
              theme(
                axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                legend.position = 'none' #remove legend
                ) +
              annotate("text", x = 3/5*(max(x)-min(x)), y = c(maxy-miny-2*(maxy/13),maxy-miny-maxy/13,maxy-miny)[1:length(state$temp1)], label = state$temp, colour = color3[1:length(state$temp1)],hjust = 0, size=5) +
              scale_color_manual(values=color3)
            }
          }
  
          myfun_rbc<-function(){
            ptype$x<-3
            myfun2(1:state$num_of_codons,unlist(state$newcd),'codonPosition','yeast_reads_based_on_codons.csv')
            x<-1:state$num_of_codons
            y<-data.frame(matrix(as.numeric(unlist(state$newcd)),ncol=length(state$temp)))
            maxy<-max(y)
            names(y)<-state$temp1
            dat<-melt(data.frame(x,y),id='x')
            ggplot(dat,aes(x=x,y=value,color=variable)) +
              geom_line() +
              ggtitle("Read counts by codon position") +
              labs(x='Codon position',y='Read count') +
              theme_classic() +
              theme(
                axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                legend.position = 'none' #remove legend
                ) +
              annotate("text", x = 3/5*(max(x)-min(x)), y = c(maxy-2*(maxy/13),maxy-maxy/13,maxy)[1:length(state$temp1)], label = state$temp, colour = color3[1:length(state$temp1)],hjust = 0, size=5) +
              scale_color_manual(values=color3)
          }
  
          myfun_rpkm<-function(databasetype){
            ii<-which(state$dtype==databasetype)
            vline<-state$vline[ii]
            lbl<-state$temp1[ii]
  
            if(databasetype=='mRNA'){dat<-dat_mrna}
            else if(databasetype=='RPF'){dat<-dat_rpf}
  
            ggplot(dat, aes_string(x=input$txt)) +
              geom_histogram(bins = input$bins, position='dodge',alpha=0.5,color=rgb(189,189,189,maxColorValue=255),fill=rgb(189,189,189,maxColorValue=255)) +
              #theme(legend.position='top') +
              ggtitle(databasetype) +
              #labs(x='Read length',y='Read count') +
              theme_classic() +
              theme(
                axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
                geom_vline(xintercept = vline,color=color3[1:length(vline)],size=1) 
                annotate("text", x = 3, y = c(1,3,5)[1:length(vline)], label = lbl[1:length(lbl)], colour = color3[1:length(vline)], hjust = 0, size=5)

          }
          #---------------------------------- plot11
          state$ddat<-list()
          ranges <- reactiveValues(x = NULL)
          ptype<-reactiveValues(x=0) #this is used so that the zoomed plot(plot22) doesn' show up until plot11 of corresponding type shows first
          output$plot11<-renderPlot({ #isolate in the renderPlot() so that it will only update when the input$clickplot button is fired
            input$clickplot
            isolate({
              if(input$select=='rbl'){
                myfun2(15:50,state$attrdat,'readLength','yeast_reads_based_on_location.csv')
                x<-15:50
                y<-data.frame(matrix(as.integer(state$attrdat),ncol=length(state$temp)))
                maxy<-max(y)
                names(y)<-state$temp1
                dat<- melt(data.frame(x, y),id='x')
                ggplot(dat, aes(x = x, y=value, fill = variable)) +
                  geom_bar(stat='identity',position = 'dodge') +
                  ggtitle("Read counts by length") +
                  labs(x='Read length',y='Read count') +
                  scale_x_continuous(breaks=seq(15,50,by=1)) +
                  theme_classic() +
                  theme(
                    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                    legend.position = 'none' #remove legend
                    ) +
                  annotate("text", x = 3/4*max(x), y = c(maxy-2*(maxy/13),maxy-maxy/13,maxy)[1:length(state$temp1)], label = state$temp, colour = color3[1:length(state$temp1)],hjust = 0, size=5) +
                  scale_fill_manual(values=color3)
                  #theme(legend.text = element_text(colour = 'red', angle = 45, size = 10, hjust = 3, vjust = 3, face = 'bold'))
              }
  
              else if(input$select=='rbp'){
                myfun_rbp()
  
              }
              else if(input$select=='rbpl'){
                if (input$txt2 %in% c(15:50)){ #check if input$txt2 is valid
                  output$errormsg<-renderText('')
                  myfun_rbpl()
                }
                else{
                  output$errormsg<-renderText('Error: please enter a length between 15-50')
                }
                
              }
  
              else if(input$select=='rbc'){
                myfun_rbc()
              }
  
              else if(input$select=='rpkm'){
                p1<-myfun_rpkm('mRNA')
                p2<-myfun_rpkm('RPF')
  
                #plot_grid(p1, p2, labels=c("A", "B"), ncol = 1, nrow = 2)
                plot_grid(p1, p2, ncol = 1, nrow = 2) #to display the two RPKM plots as one plot
  
              }
            })
          }  )
    #---------------------------------------------------------------------plot22: zooming plots
          output$plot22<-renderPlot({
              if(input$select=='rbp' & ptype$x==1){
                myfun_rbp() +
                  coord_cartesian(xlim = ranges$x)
  
              }
              else if(input$select=='rbpl' & ptype$x==2){
                myfun_rbpl() +
                  coord_cartesian(xlim = ranges$x)
              }
  
              else if(input$select=='rbc' & ptype$x==3){
                myfun_rbc() +
                  coord_cartesian(xlim = ranges$x)
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
  
    #--------------------------------UIplot: multiple plots for 'reads by length' plot
          # Insert the right number of plot output objects into the web page
          output$plots <- renderUI({
            plot_output_list <- lapply(1:length(state$temp), function(i) {
              plotname <- paste("plot", i, sep="")
              plotOutput(plotname)
              #plotOutput(plotname, height = 320, width = 700)
            })
  
            do.call(tagList, plot_output_list)# Convert the list to a tagList - this is necessary for the list of items to display properly.
          })
  
          # Call renderPlot for each one. Plots are only actually generated when they are visible on the web page.
          for (i in 1:length(state$temp)) {
            # Need local so that each item gets its own number. Without it, the value of i in the renderPlot() will be the same across all instances, because of when the expression is evaluated.
            local({
              my_i <- i
              plotname <- paste("plot", my_i, sep="")
  
              output[[plotname]] <- renderPlot({
                input$clickplot
                isolate({
                  if(input$select=='rbl'){
                    #barplot(state$attrdat[(36*(my_i-1)+1):(36*my_i)],main='Read counts by Length',xlab='Read length',ylab='Read count',names.arg = c(15:50),col='blue')
                    x<-15:50
                    dat<-data.frame(x,y=state$attrdat[(36*(my_i-1)+1):(36*my_i)])
                    ggplot(dat, aes(x=x, y=y)) +
                      geom_bar(stat='identity') +
                      ggtitle("Read counts by length") +
                      labs(x='Read length',y='Read count') +
                      scale_x_continuous(breaks=seq(15,50,by=1)) +
                      theme_classic() +
                      theme(
                        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')
                        )
  
                  }
                })
              })
            })
          }
        }
        #--------------------
        else{
          output$errormsg<-renderText('Error: the gene you entered does not exist in the database, please re-enter')
          output$plot11<-renderPlot('')
        }
        #--------------------
      }
    }
  )
  reactive({H5close()})
}

shinyApp(ui = shinyUI, server = shinyServer)