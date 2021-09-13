# Provenance-related utilities.

suppressMessages(library(getopt))
suppressMessages(library(git2r))

#' Get RiboViz version information.
#'
#' If the file is within the scope of a Git repository then a message
#' including the given file name, Git commit hash and date of HEAD is
#' returned.
#'
#' If the file is not within the scope of a Git repository then an
#' "unknown" message is returned.
#'
#' @param file_path R file path.
#' @return: message.
get_version <- function(file_path = getopt::get_Rscript_filename()) {
  location <- dirname(normalizePath(file_path))
  version <- tryCatch({
    repo <- git2r::repository(location)
    commit <- git2r::commits(repo, n = 1)[[1]]
    sha <- commit$sha
    time <- commit$author$when
    paste("commit", sha, "date", time)
  },
  error = function(cond) {
    version <- "unknown"
    return(version)
  })
  return(version)
}


#' Write a provenance header to a file with RiboViz date and
#' version information.
#'
#' @param file_path Path to R file invoking this function.
#' @param data_file Path to file into which header is to be written.
write_provenance_header <- function(file_path, data_file) {
  conx <- file(data_file, open = "w")
  writeLines("# Created by: RiboViz", conx)
  writeLines(paste("# Date:", toString(Sys.time())), conx)
  writeLines(paste("# File:", file_path), conx)
  writeLines(paste("# Version:", get_version(file_path)), conx)
  close(conx)
}

#' Print provenance information with RiboViz date and version
#' information.
#'
#' @param file_path Path to R file invoking this function.
print_provenance <- function(file_path) {
  print("Created by: RiboViz")
  print(paste("Date:", toString(Sys.time())))
  print(paste("File:", file_path))
  print(paste("Version:", get_version(file_path)))
}
