library(shiny)
library(rhdf5)
library(ggplot2)
library(reshape2)
library(plotly)
library(RColorBrewer)
library(data.table) # To use rbindlist()
library(tidyr) # To use gather()

load("./riboViz.RData")     # Contains d1:F8_RPKMs_modified.tsv, d2:bgdb.tsv, df:data_unq.tsv

shinyServer <- function(input, output){
  path <- "./"
  yseq = './yeast_seq.h5'   # Stores nucleotide and codon sequences for all yeast genes
  codonseq = reactive({h5read(yseq,paste0(input$slct_gene, '/codon'))}) # Read in codons sequence for all genes
  ntseq = reactive({h5read(yseq,paste0(input$slct_gene, '/nt'))}) # Read in nucleotide sequence for all genes
  
  dbslash <- reactive({gsub(' ','/', input$slct_db)})          # e.g. "2016/Weinberg/unselected_total_RNA"
  myColors <- brewer.pal(9, "Set1")
  igene <- reactive({which(names(d2) == input$slct_gene)})      # Indexes for same gene in d1 and d2 are the same
  dat_mrna <- reactive({d2[d2$Type == 'mRNA', c(5,igene())]})    # RPKM values for selected gene for background mRNA databases 
  dat_rpf <- reactive({d2[d2$Type == 'RPF', c(5,igene())]})      # RPKM values for selected gene for background RPF databases 
  
  # To display legend texts with corresponding colors
  output$txt_legends <- renderText({ paste(unlist(lapply(1 : length(input$slct_db), function(i){
                        paste('<font color=\"', myColors[i],'\"><b>', input$slct_db[i],'<br>', '</b></font>')
      })
    ))
  })
  
  observe({
    if (length(input$slct_db)==0){ # Check if database selection is empty
      output$txt_errmsg <- renderText('Error: Please select at least one database')
    }
    else{ # If database selection isn't empty
      # Stores the intermediate reactive values needed for making the plot
      # state$geoid               e.g. "GSM1969533"
      # state$inf                 e.g. "./2016/Weinberg/unselected_total_RNA/GSM1969533.h5"
      # state$transdb             e.g. "2016_Weinberg_unselected_total_RNA"
      # state$attrdat             e.g. 0 0 0 0 0 0 0 1 2 1 16 35 41 33 47 58 126 221 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      # 
      # state$coord               Nucleotide coordinate including left(250) and right(247) buffer, e.g.3980 for YAL001C
      # state$d                   e.g. for YAL001C is a 36x3980 matrix
      # state$dd                  Excluding left and right buffer regions, e.g. for YAL001C is a 36x(3980-250-247)=36x125388 matrix
      # state$dbuffer             e.g. for YAL001C is a 36x(3980-250-247+25+25) matrix, with 25 nucleotides upstream start codon and 25 nt downstream stop codon
      # state$sumd                Stores reads for all length, all positions for state$d
      # state$sumdbuffer          Stores reads for all length, all positions for state$dbuffer
      # state$meanreadsrbp        Mean reads = reads total/total positions
      # state$meanreadsrbpl
      # state$meanreadsrbc
      # state$cd                  3 * 3483, matrix to calculate the ribosome density
      # state$num_of_codons       Number of codons = integer
      # state$newcd               9 * (3483/3) matrix--->1161 * 1 (after summing up the 9 rows in a column)
      # 
      # state$vline               RPKM values for selected databases
      # state$dtype               Database type - mRNA or RPF
      
      output$txt_errmsg <- renderText("")
      state <- reactiveValues()
      state$tempstr <- sapply(dbslash(), function(x){strsplit(x, '/')})  # splitted string of dbslash
      state$geoid <- as.character(sapply(state$tempstr, function(x){df[(df$year==x[1])&(df$author==x[2])&(df$database==x[3]), 11]}))
      state$inf <- as.character(mapply(function(x,y){paste0(path, x, '/', y, '.h5')}, dbslash(), state$geoid))
      state$transdb <- as.character(sapply(state$tempstr, function(x){paste(x[1], x[2], x[3], sep = "_")}))
      state$attrdat <- as.numeric(mapply(function(x,y){h5readAttributes(file = x, name = paste0(input$slct_gene, '/', y, '/reads'))$reads_by_len},
                                         state$inf, state$transdb))
      
      state$coord <- dim(h5read(state$inf[1], paste0(input$slct_gene, '/', state$transdb[1], '/reads/data')))[2]
      state$d <- mapply(function(x,y){return(data.frame(h5read(x, paste0(input$slct_gene, '/', y, '/reads/data'))))}, state$inf, state$transdb, SIMPLIFY = FALSE)
      state$dbuffer <- mapply(function(x,y){return(data.frame(h5read(x, paste0(input$slct_gene, '/', y, '/reads/data'))[,(251-25):(state$coord-247+25)]))},
                              state$inf, state$transdb, SIMPLIFY = FALSE)
      state$sumd <- lapply(state$d, function(x){lapply(x, sum)})
      state$sumdbuffer <- lapply(state$dbuffer, function(x){lapply(x, sum)})
      state$meanreadsrbp <- as.numeric(lapply(state$dbuffer, function(x){sum(x) / (state$coord-250-247+25+25)}))
      
      state$cd <- lapply(state$d, function(x){return(data.frame(rbindlist(list(x[14:15,(251-15) : (state$coord-247-15)], x[16,(251-16) : (state$coord-247-16)]))))})
      state$num_of_codons <- (state$coord-247-251+1) / 3
      state$newcd <- lapply(state$cd, function(x){lapply(data.frame(matrix(unlist(x), nrow = nrow(x)*3, ncol = state$num_of_codons)), sum)})
      
      state$vline <- as.numeric(sapply(state$tempstr, function(x){d1[d1$Year==x[1] & d1$Author==x[2] & d1$Dataset==x[3], igene()]}))
      state$dtype <- as.character(sapply(state$tempstr, function(x){d1[d1$Year==x[1] & d1$Author==x[2] & d1$Dataset==x[3], 4]}))
      
      # state$downdata_rbl <- matrix()     # Initiate the download content for each plot here b/c the data will be read in within functions
      # state$downdatar_bpl <- matrix()
      # state$downdata_rbc <- matrix()
      # state$downdata_rpkm <- matrix()
      # state$downdata_rpkmdens <- matrix()
      
      ###############################################################
      # Functions
      ###############################################################
      myfun_theme <- function(){
        theme_classic() +
          theme(
            axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.text = element_text(size=8), #axis label size
            axis.ticks.length = unit(0.2,"cm")) #axis ticks size

      }
      
      #--------------------------------------------------myfun_ggplot
      myfun_ggplot <- function(data,labx,laby){
        ggplot(data,aes(x = x,y = value,color = variable)) +
          geom_line() +
          labs(x = labx,y = laby) +
          myfun_theme() + 
          scale_color_manual(values = myColors) 
      }
      
      #--------------------------------------------------myfun_rbpl
      # Reads by position and length regular reads
      myfun_rbpl <- function(){
        if (input$slct_readlen == "All read lengths"){
          x <- 1:state$coord         # x-axis coordinates
          y <- data.frame(matrix(unlist(state$sumd),ncol = length(input$slct_db)))      # y-axis coordinates
        }
        
        else{ # If selected read length = 15:50
          state$ddat <- lapply(state$d, function(x){x[(as.numeric(input$slct_readlen) - 14),]})      # Subset data for selected read length
          x <- 1:state$coord
          y <- data.frame(matrix(as.integer(unlist(state$ddat)),ncol = length(dbslash())))
        }
        
        dat_rename <- data.frame(x, y)
        colnames(dat_rename) <- c("x",input$slct_db)
        dat<- melt(dat_rename, id='x')
        
        p <- myfun_ggplot(dat,'Read position','Read count') +
          geom_rect(aes(xmin = 251, xmax = (max(x)-247), ymin = 0, ymax = max(value), col = NA),   # Add grey shading to the coding region
                    fill = "gray",alpha = 0.2) 
        p <- ggplotly(p, dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 370, legend=list(orientation= "h"))
        
        
        #--------------- Create content to send to download handler
        temp <- matrix(c(x,unlist(y)),ncol = (length(input$slct_db))+1) # Intermediate matrix for adding the colnames
        colnames(temp) <- c("Read length",paste("Regular_",input$slct_db))
        state$downdata_rbpl <- temp    # Store the downloading content to state$downdata_rbpl
        #---------------
        
        p <- plotly_build(p)    # To override defaults provided by ggplotly
        
        for (i in 1:length(input$slct_db)){
          p$x$data[[i]]$name <- " "            # Change the legend text to empty
          
          p$x$data[[i]]$text <- paste("Read position:", dat$x[dat$variable==input$slct_db[i]], "<br>",               # Customize the content of tooltip
                                      "Read count:", round(dat$value[dat$variable==input$slct_db[i]],3), "<br>",
                                      "Database:", input$slct_db[i], "<br>",
                                      "Sequence:", ntseq())
        }
        
        p$x$data[[length(input$slct_db)+1]]$showlegend <- FALSE   # Turn off legend for the shading (b/c geom_rect is treated as a variable)
        
        p$x$layout$legend$y <- 1.2             # Set the legend above plot, range -2 to 3
        p
      }
      
      
      #--------------------------------------------------myfun_rbpl_norm
      # Reads by position and length normalized reads
      myfun_rbpl_norm <- function(){
        x <- 1:state$coord
        tempmatrix <- matrix(0, ncol = length(input$slct_db), nrow = state$coord) # Matrix filled with 0s
        
        if (input$slct_readlen == "All read lengths"){
          if(any(state$meanreadsrbp == 0)){
            ixx <- which(state$meanreadsrbp != 0)
            tempmatrix[,ixx] <- t(t(matrix(as.integer(unlist(state$sumd)),ncol = length(input$slct_db))[,ixx])/state$meanreadsrbp[ixx]) # Fill the matrix with normy for those meanreads !=0
            normy <- data.frame(tempmatrix)
          }
          
          else{       # If none of the selected database have meanreads that equals 0
            normy <- data.frame(t(t(matrix(as.integer(unlist(state$sumd)),ncol = length(input$slct_db)))/state$meanreadsrbp))
          }
          
          dat_rename <- data.frame(x, normy)
          colnames(dat_rename) <- c("x",input$slct_db)
          dat<- melt(dat_rename, id='x')
          
          p <- myfun_ggplot(dat,'Read position','Normalized reads') +
            geom_rect(aes(xmin = 251, xmax = (max(x)-247), ymin = 0, ymax = max(value), col = NA), fill = "gray",alpha = 0.2) # Add shading to the coding region
          
          p <- ggplotly(p, dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 370, legend=list(orientation= "h"))
        }
        
        
        else{       # If selected read length == 15:50
          state$ddat <- lapply(state$d, function(x){x[(as.numeric(input$slct_readlen)-14),]})  # e.g. 3980 x 1
          state$ddatbuffer <- lapply(state$dbuffer, function(x){x[(as.numeric(input$slct_readlen) - 14),]})    # e.g. 3533 x 1
          state$meanreadsrbpl <- as.numeric(lapply(state$dbuffer, 
                                                   function(x){sum(x[isolate((as.numeric(input$slct_readlen) - 14)),])/(state$coord-250-247+25+25)}))
          
          if(any(state$meanreadsrbpl) == 0){   # If mean reads for any selected databases is 0
            ixx <- which(state$meanreadsrbpl != 0)   # Index for databases that has mean reads not equal to 0
            tempmatrix[,ixx] <- t(t(matrix(as.integer(unlist(state$ddat)),     # Fill the matrix with normy for those meanreads !=0
                                           ncol=length(input$slct_db))[,ixx]) / state$meanreadsrbpl[ixx]) 
            normy <- data.frame(tempmatrix)
          }
          
          else{         # If none of the selected database have meanreads(at selected read length) that equals 0
            normy <- data.frame(t(t(matrix(as.integer(unlist(state$ddat)), ncol = length(input$slct_db))) / state$meanreadsrbpl))
          }
          
          dat_rename <- data.frame(x, normy)
          colnames(dat_rename) <- c("x",input$slct_db)
          dat<- melt(dat_rename, id='x')
          
          p <- myfun_ggplot(dat,'Read position','Normalized reads') +
            geom_rect(aes(xmin = 251, xmax = (max(x)-247), ymin = 0, ymax = max(value), col = NA), fill = "gray",alpha=0.2) # Add shading to the coding region
          
          p <- ggplotly(p, dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 370, legend = list(orientation = "h"))
        }
        
        #--------------- Create content to send to download handler
        temp <- matrix(c(x,unlist(normy)),ncol = (length(input$slct_db))+1) # Intermediate matrix for adding the colnames
        colnames(temp) <- c("Read length",paste("Normalized_",input$slct_db))
        state$downdata_rbpl <- temp
        #---------------
        
        p <- plotly_build(p)
        
        for (i in 1:length(input$slct_db)){
          p$x$data[[i]]$name <- " "            # Change the legend text to empty, only legend labels remaining
          
          p$x$data[[i]]$text <- paste("Read position:", dat$x[dat$variable==input$slct_db[i]], "<br>",               # Customize the content of tooltip
                                      "Normalized read:", round(dat$value[dat$variable==input$slct_db[i]],3), "<br>",
                                      "Database:", input$slct_db[i], "<br>",
                                      "Sequence:", ntseq())
        }
        
        p$x$data[[length(input$slct_db)+1]]$showlegend <- FALSE   # Turn off the legend for the shading (b/c geom_rect is treated as a variable)
        p$x$layout$legend$y <- 1.2             # Set the legend above plot, range -2 to 3
        p
      }
      
      #--------------------------------------------------myfun_rbc
      # Reads by codon regular reads
      myfun_rbc<-function(){
        x <- 1:state$num_of_codons
        y <- data.frame(matrix(as.numeric(unlist(state$newcd)),ncol = length(dbslash())))
        maxy <- max(y)
        
        dat_rename <- data.frame(x, y)
        colnames(dat_rename) <- c("x",input$slct_db)
        dat<- melt(dat_rename, id='x')
        
        p <- myfun_ggplot(dat,'Codon position','Read count') 
        p <- ggplotly(p, dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 370, legend=list(orientation= "h"))
        
        
        #--------------- Create content to send to download handler
        temp <- matrix(c(x,unlist(y)),ncol = (length(input$slct_db))+1) # Needs a intermediate matrix to add the colnames
        colnames(temp) <- c("Read length",paste("Regular_",input$slct_db))
        state$downdata_rbc <- temp
        #---------------
        
        p <- plotly_build(p)
        for (i in 1:length(input$slct_db)){
          p$x$data[[i]]$name <- " "            # Change the legend text to empty
          
          p$x$data[[i]]$text <- paste("Codon position:", dat$x[dat$variable==input$slct_db[i]], "<br>",               # Customize the content of tooltip
                                      "Read count:", round(dat$value[dat$variable==input$slct_db[i]],3), "<br>",
                                      "Database:", input$slct_db[i], "<br>",
                                      "Codon sequence:", codonseq())
        }
        
        p$x$layout$legend$y <- 1.2             # Set the legend above plot, range -2 to 3
        p
      }
      
      #--------------------------------------------------myfun_rbc_norm
      # Reads by codon normalized reads
      myfun_rbc_norm <- function(){
        state$meanreadsrbc <- apply(matrix(as.numeric(unlist(state$newcd)), 
                                           ncol = length(input$slct_db)), 2, sum)/(state$num_of_codons)
        x <- 1:state$num_of_codons
        tempmatrix <- matrix(0, ncol = length(input$slct_db), nrow = state$num_of_codons) # Matrix filled with 0s
        
        if(any(state$meanreadsrbc) == 0){
          ixx <- which(state$meanreadsrbc != 0)
          tempmatrix[,ixx] <- t(t(matrix(as.numeric(unlist(state$newcd)),ncol = length(input$slct_db))[,ixx]) / state$meanreadsrbc[ixx])
          normy <- data.frame(tempmatrix)
        }
        
        else{         # If none of the selected database have meanreads that equals 0
          normy <- data.frame(t(t(matrix(as.numeric(unlist(state$newcd)), ncol = length(input$slct_db))) / state$meanreadsrbc))
        }
        
        dat_rename <- data.frame(x, normy)
        colnames(dat_rename) <- c("x",input$slct_db)
        dat<- melt(dat_rename, id='x')
        
        p <- myfun_ggplot(dat,'Codon position','Normalized reads')
        p <- ggplotly(p, dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 370, legend = list(orientation = "h"))
        
        
        #--------------- Create content to send to download handler
        temp <- matrix(c(x,unlist(normy)),ncol=(length(input$slct_db))+1) # Intermediate matrix to add the colnames
        colnames(temp) <- c("Read length",paste("Normalized_",input$slct_db))
        state$downdata_rbc <- temp
        output$txt_test<-renderPrint(typeof(state$downdata_rbc))
        
        #---------------
        
        p <- plotly_build(p)
        
        for (i in 1:length(input$slct_db)){
          p$x$data[[i]]$name <- " "            # Change the legend text to empty
          
          p$x$data[[i]]$text <- paste("Codon position:", dat$x[dat$variable==input$slct_db[i]], "<br>",               # Customize the content of tooltip
                                      "Normalized read:", round(dat$value[dat$variable==input$slct_db[i]],3), "<br>",
                                      "Database:", input$slct_db[i], "<br>",
                                      "Codon sequence:", codonseq())
        }
        p$x$layout$legend$y <- 1.2             # Set the legend above plot, range -2 to 3
        p
      }
      
      #--------------------------------------------------myfun_rpkm
      # Reads Per Kilobase of transcript per Million mapped reads
      
      myfun_rpkm <- function(dbtype){              # Designated dbtype can be either mRNA or RPF
        ii <- which(state$dtype == dbtype)         # Indexes for selected database that is of designated dbtype
        vline <- state$vline[ii]                   # Vertical line = RPKM value for selected database of designated dbtype
        
        if(dbtype == 'mRNA'){dat <- dat_mrna()}
        else if(dbtype == 'RPF'){dat <- dat_rpf()}
        
        p <- ggplot(dat, aes_string(x = input$slct_gene)) +
          geom_histogram(bins = input$bin_rpkm, position = 'dodge',alpha = 0.5,color = rgb(189,189,189,maxColorValue = 255), 
                         fill = rgb(189,189,189, maxColorValue = 255)) +
          xlab(paste0("wt RPKM (",dbtype,")")) +
          ylab(if(dbtype=="mRNA"){"Count"}) +
          myfun_theme() 

        if (length(ii) == 0){    # If there's no database selected for desiganted dbtype, plot only the background
          p <- plotly_build(p)
          p$x$data[[1]]$text <- ""
        }   
        
        else {     # Otherwise plot both background and the vertical lines
          p<- p + geom_vline(xintercept = vline, color = myColors[ii], size = 1)
          p<-plotly_build(p)
          
          p$x$data[[1]]$text <- "" 
          
          for (i in 2:(length(ii) + 1)){          # This keeps messing with the $data index
            # To add database names, it's hard to read it from the plotly_build annotations, here is to read it from 
            # input$slct_db for which the index matches the index for vline in state$vline
            p$x$data[[i]]$text <- paste("RPKM:", unlist(strsplit(p$x$data[[i]]$text, ":"))[2], "<br>",
                                        "Database:", input$slct_db[which(state$vline == as.numeric(unlist(strsplit(p$x$data[[i]]$text, ":"))[2]))]

            )
          }
        }
        
        return(p)
      }
      ###############################################################
      # Making plots
      ###############################################################
      
      output$plotrbl<-renderPlotly({
        x <- 15:50
        y <- data.frame(prop.table(matrix(unlist(state$attrdat), ncol = length(input$slct_db)), 2))
        temp <- matrix(c(x, state$attrdat), ncol = (length(input$slct_db)) + 1) # Intermediate matrix for adding the colnames
        colnames(temp) <- c("Read length", input$slct_db)
        state$downdata_rbl <- temp
        
        maxy <- max(y) 
        dat_rename <- data.frame(x, y)
        colnames(dat_rename) <- c("x",input$slct_db)
        dat<- melt(dat_rename, id='x')
        
        p <- ggplot(dat, aes(x = x, y = value, fill = variable)) +
             geom_bar(stat = 'identity', position = 'dodge') +
             labs(x = 'Read length', y = 'Frequency of read lengths') +
             scale_x_continuous(breaks = seq(15,50, by = 1)) +
          
          theme_classic() +
          myfun_theme() +
          scale_fill_manual(values = myColors) 
        
        p <- ggplotly(p, dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 370, legend=list(orientation= "h"))
        p <- plotly_build(p)

        for (i in 1 : length(input$slct_db)){
          p$x$data[[i]]$name <- " "            # Change the legend text to empty
          p$x$data[[i]]$text <- paste("Read length:", dat$x[dat$variable==input$slct_db[i]], "<br>",               # Customize the content of tooltip
                                      "Frequency:", round(dat$value[dat$variable==input$slct_db[i]], 3), "<br>",
                                      "Database:", input$slct_db[i])
        }
        p$x$layout$legend$y <- 1.2             # Set the legend above plot, range -2 to 3
        p
      })
      
      output$plotrbpl <- renderPlotly({
        if(input$btn_rbpl == 1){
          myfun_rbpl()
        }
        else if(input$btn_rbpl == 2){
          myfun_rbpl_norm()
        }
      })
      
      output$plotrbc <- renderPlotly({
        if(input$btn_rbc == 1){
          myfun_rbc()
        }
        else if(input$btn_rbc == 2){
          myfun_rbc_norm()
        }
      })
      
      output$plotrpkm <- renderPlotly({
        p1 <- ggplotly(myfun_rpkm("mRNA"), dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 250)
        p2 <- ggplotly(myfun_rpkm("RPF"), dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 250)
        subplot(p1, p2, titleX = T, titleY = T)
      })
      
      output$plotrpkmdens <- renderPlotly({
        rpkm_all_gene <- data.frame(sapply(state$tempstr, function(x){d1[d1$Year==x[1] & d1$Author==x[2] & d1$Dataset==x[3],6:ncol(d1)]})) 
        # 5293(number of genes) X number of selected databases 
        
        state$downdata_rpkm_dens <- t(rpkm_all_gene)  # Create download content
        
        melted_rpkm <- rpkm_all_gene %>% gather(dbname, RPKM,1 : length(dbslash())) 
        # Convert the data frame to long format, col1=database names, col2=RPKM values
        
        melted_rpkm[2] <- lapply(melted_rpkm[2], as.numeric)   # Convert the 2nd column (RPKM values) to numeric
        p <- ggplot(melted_rpkm, aes(x = RPKM, fill = dbname)) + 
             geom_density(alpha = 0.2) + 
             scale_x_log10() +    # Covert to log data b/c original data is highly skewed
             theme_classic() +
             theme(legend.title = element_blank()) +
             myfun_theme() +
             scale_fill_manual(values = myColors) +
             geom_vline(xintercept = state$vline, color = myColors[1 : length(dbslash())], size = 1)
        
        p<-ggplotly(p, dynamicTicks = TRUE)%>%layout(autosize = F, width = 650, height = 370, legend=list(orientation= "h"))
        
        p<-plotly_build(p)
        
        p$x$layout$xaxis$title <-"log10 (RPKM)"     # Override the x-axis title
        for (i in 1 : length(input$slct_db)){    # e.g. if 3 database is selected, the first 3 data refers to the background density plot
          p$x$data[[i]]$name <- " "              # Change the legend text to empty
          p$x$data[[i]]$text <- ""               # Change the tooltip text to empty for background density plot
        }
        
        for (i in (length(input$slct_db) + 1) : (2 * length(input$slct_db))){  # e.g. if 3 database is selected, i=6, the last 3 data refers to the vertical lines
          
          # To add database names, it's hard to read it from the plotly_build annotations, here is to read it from
          # input$slct_db for which the index matches the index for vline in state$vline
          p$x$data[[i]]$text <- paste("log(RPKM):", unlist(strsplit(p$x$data[[i]]$text,":"))[2], "<br>",
                                      "Database:", input$slct_db[which(round(log10(state$vline), 6) == as.numeric(unlist(strsplit(p$x$data[[i]]$text,":"))[2]))]

          )
        }
        
        
        p$x$layout$legend$y <- 1.2               # Set the legend above plot, range -2 to 3
        p
      })
      
      
      ###############################################################
      # Download handler
      ###############################################################
      
      output$dwld_data <- downloadHandler(
        filename = function() {
          paste0(input$slct_gene, "_currentData.zip")
        },
        content = function(fname) {
          fs <-c("reads_by_length.csv","reads_by_position_length.csv","reads_by_codon.csv","rpkm_background.csv",
                 "rpkm_selected.csv", "rpkm_density_background.csv"
                 ) 
          # Shouldn't paste the gene name here as it'll create too many downloaded files on server
          tmpdir <- tempdir()         # A copy of downloaded files will be saved here, but kept being replaced
          setwd(tempdir())
          write.csv(state$downdata_rbl, fs[1])
          write.csv(state$downdata_rbpl, fs[2])
          write.csv(state$downdata_rbc, fs[3])
          write.csv(d2[,c(1,2,3,4,igene())], fs[4]) # RPKM for background databases
          write.csv(data.frame(Databases = input$slct_db, Type = state$dtype, RPKM = state$vline), fs[5]) # RPKM for selected databases
          write.csv(state$downdata_rpkm_dens, fs[6]) # RPKM density for background databases
          
          
          zip(zipfile=fname, files=fs)
        },
        contentType = "application/zip"
      )
      
      
    }
  })
  
  reactive({H5close()})
}
