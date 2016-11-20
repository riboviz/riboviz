library(shiny)

f<-read.table('dbnames.txt')

shinyUI<-fluidPage(
  titlePanel('Dynamic user interface -renderUI'),
  sidebarLayout(
    sidebarPanel(
      
      #select input with the list of datasets
      selectInput(inputId = 'data1',label='Select the dataset of your choice',
                  choices=c('iris','mtcars','trees')),
      
      uiOutput('vx'), #vx is coming from renderUI in server.r
      uiOutput('vy') #vy is coming from renderUI in server.r
    ),
    mainPanel(
      plotOutput('p')
    )
  )
)

shinyServer<-function(input,output){
  var<-reactive({
    switch(input$data1,
           'iris'=names(iris),
           'mtcars'=names(mtcars),
           'trees'=names(trees)
    )
  })
  #vx and vy are the output variable from renderUI containing the list of variable names in a dropdown/selectInput UI
  #renderUI is used in server side and is used along with uiOutput in the ui.r
  #uiOutput is used in ui.R to display the selectInput widget dynamically using the output variable vx&vy
  #using vx and vy in the ui.R we will dynamically create the drop down with the column names based on the dataset selected
  
  output$vx<-renderUI({
    selectInput('variablex','Select the First(X) variable',choices=var())
  })
  output$vy<-renderUI({
    selectInput('variabley','Select the second (Y) variable', choices=var())
  })
  
  #renderPlot is used to plot the ggplot and the plot output will be stored in the output variable p
  output$p <- renderPlot({
    attach(get(input$data1))
    plot(x=get(input$variablex),y=get(input$variabley),xlab=input$variablex,ylab=input$variabley)
  })
}

shinyApp(server=shinyServer,ui=shinyUI)
