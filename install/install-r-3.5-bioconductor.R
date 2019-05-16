r = getOption("repos")
r["CRAN"] = "https://cloud.r-project.org" 
options(repos = r)
install.packages("BiocManager")
BiocManager::install("Rsamtools", lib=Sys.getenv("R_LIBS_USER"))
BiocManager::install("rtracklayer", lib=Sys.getenv("R_LIBS_USER"))
BiocManager::install("rhdf5", lib=Sys.getenv("R_LIBS_USER"))
