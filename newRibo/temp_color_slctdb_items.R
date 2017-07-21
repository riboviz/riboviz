
library("shiny")

# appCSS <- 
#   "#color ~ .selectize-control.single .selectize-dropdown [data-value=blue] { color: blue }
#  #color ~ .selectize-control.single .selectize-dropdown [data-value=red] { color: red }"
# 
# runApp(shinyApp(
#   ui = fluidPage(
#     tags$head(tags$style(HTML(appCSS))),
#     selectInput("color", "Color", c("blue", "red"))
#   ),
#   server = function(input, output, session) {
#   }
# ))

####################################################

runApp(shinyApp(

 ui = fluidPage(

  # tags$head(
  #   tags$style(HTML("
  #                       .item {
  #                       background: #2196f3 !important;
  #                       }
  #                       .selectize-dropdown-content .active {
  #                       background: #2196f3 !important;
  #                       }
  #                       "))
  # ),

   tags$style(HTML(" .item {
                    background: #2196f3 !important;
                    }

                    #select2 + div> div>.item {
                    background:   #f3217a !important;
                    }
                    .selectize-dropdown-content .active {
                    background: #2196f3 !important;
                    }

                     #select2 + div> div>.selectize-dropdown-content .active {
                    background:   #f3217a !important;
                    }
                    ")),

  sidebarLayout(
    sidebarPanel(
      selectizeInput("select", label=NULL,
                     choices=c("a", "b", "c", "d"),
                     multiple=TRUE, options=list(placeholder="Wybierz")),
      selectizeInput("select2", label=NULL,
                     choices=c("a", "b", "c", "d"),
                     multiple=TRUE, options=list(placeholder="Wybierz"))),

    mainPanel())
),
  server = function(input, output, session) {
     }

))


######################################################

# shinyApp(
#   ui = fluidPage(
#     tags$head(
#       tags$style(HTML("
#           .item [data-value=\"Banana\"]{
#                             background: #2196f3 !important;
#                             color: yellow !important;
#                             }
#           .item [data-value=\"Tomato\"]{
#                             background: #2196f3 !important;
#                             color: red !important;
#                             }
#           .item [data-value=\"Kiwi\"]{
#                             background: #2196f3 !important;
#                             color: green !important;
#                             }
#                   "))
#     ),uiOutput("type")),
#   
#   server = function(input, output, session) {
#     output$type <- renderUI({
#       selectInput("color", "Color",as.list(fruits),multiple = T)
#     })
#   }
# )


  