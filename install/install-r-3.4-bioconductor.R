r = getOption("repos")
r["CRAN"] = "https://cloud.r-project.org" 
options(repos = r)
source("https://bioconductor.org/biocLite.R")
biocLite("Rsamtools", lib=Sys.getenv("R_LIBS_USER"))
biocLite("rtracklayer", lib=Sys.getenv("R_LIBS_USER"))
biocLite("rhdf5", lib=Sys.getenv("R_LIBS_USER"))
