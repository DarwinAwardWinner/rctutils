#' Find the first accessible file path from a vector of paths.
#'
#' This function is useful for finding a file that may be in one of
#' several different locations.
#'
#' @param paths
#'     "A character vector of paths to check for accessibility, in order."
#' @param mode "Has the same meaning as in \code{file.access()}."
#' @return The first element of \code{paths} for which
#'     \code{file.access()} returns TRUE. If none of \code{paths} is
#'     accessible, NA is returned.
#' @examples
#' exec_paths <- strsplit(Sys.getenv("PATH"), ":")[[1]]
#' first_accessible(exec_paths)
#' @export
first_accessible <- function(paths, mode=0) {
    for (path in paths) {
        if (file.access(path, mode) == 0) {
            return(path)
        }
    }
    return(NA_character_)
}
