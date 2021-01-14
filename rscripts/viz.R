#visualisation script
rmarkdown::render(
  input = "rmarkdown/new_visualization.Rmd",
  output_file = "Dashboard.html"
)