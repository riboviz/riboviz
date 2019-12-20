# Provenance-related utilities.

library(getopt)
library(git2r)

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
  file_name <- basename(file_path)
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
  return(paste(file_name, "version", version))
}
