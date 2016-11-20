library(shiny)
library(rhdf5)
inf='/Users/xtj/Documents/shahlab/shinyProject/data/test.h5'
shinyUI <- fluidPage(
  titlePanel('Explore the Ribosomal Profiling Databases'),
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = 'txt',label='Write the gene name of interest',value='YAL001C',placeholder = 'e.g. YAL001C'),
      selectInput(inputId = 'select',label='Select a type of plot',
                  choices = c('Reads by length'='rbl','Reads by position'='rbp','Reads by position and length'='rbpl','Reads by codon'='rbc')),
      conditionalPanel(
        condition = "input.select == 'rbpl'",
        textInput(inputId = 'txt2',label='Write the read length of interest',value='15',placeholder='ranging from 15-50')),
      submitButton('Update View')
    ),
    mainPanel(
      plotOutput('plot')
    )
  )
)

shinyServer <- function(input,output){
  
  attrloc<-reactive({h5readAttributes(file=inf,name=paste(input$txt,'/2016_Weinberg_RPF/reads',sep=''))})
  group<-reactive({paste(input$txt,'/2016_Weinberg_RPF/reads/data',sep='')})
  coord<-reactive({dim(h5read(inf,group()))[2]}) #nucleotide coordinate e.g.3980 for YAL001C
  endcoord<-reactive({coord()-247}) #end coordinate for the last nucleotide in the ORF
  d<-reactive({h5read(inf,group())}) #e.g. for YAL001C is a 36x3980 matrix
  sumd<-reactive({lapply(data.frame(d()),sum)}) #stores reads for all length, all positions
  
  d1<-reactive({d()[14:15,(251-15):(endcoord()-15)]}) #matrix 2 rows (ribosome protected length=28,29),columns=endcoord-251+1
  d2<-reactive({d()[16,(251-16):(endcoord()-16)]}) #matrix 1 row (ribosome protected length=30)
  cd<-reactive({rbind(d1(),d2())}) #3X7482, matrix to calculate the ribosome density
  num_of_codons<-reactive({(endcoord()-251+1)/3}) #number of codons= integer
  cd2<-reactive({matrix(cd(),nrow=nrow(cd())*3,ncol=num_of_codons())})
  newcd<-reactive({lapply(data.frame(cd2()),sum)}) 
  
  output$plot<-renderPlot({
    if(input$select=='rbl'){
      plot({attrloc()$lengths},{attrloc()$reads_by_len},type='l',main='Read counts by Length',xlab='Read length',ylab='Read count')}
    
    else if(input$select=='rbp'){
      plot(1:coord(),sumd(),type='l',main='Read counts by Position',xlab='Read position',ylab='Read count')}
    
    else if(input$select=='rbpl'){
      plot(1:coord(),d()[(as.numeric(input$txt2)-14),],type='l',main='Read counts by Position and length',xlab='Read position',ylab='Read count')}
    
    else if(input$select=='rbc'){
      plot(1:num_of_codons(),newcd(),type='l',main='Read counts by codon position',xlab='Codon position',ylab='Read count')}
  })
  reactive({H5close()})
}

shinyApp(server=shinyServer,ui=shinyUI)