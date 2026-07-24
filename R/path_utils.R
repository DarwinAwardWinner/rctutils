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
#' first_accessible_path(exec_paths)
#' @export
first_accessible_path <- function(paths, mode = 0) {
    for (path in paths) {
        if (file.access(path, mode) == 0) {
            return(path)
        }
    }
    return(NA_character_)
}

#' Remove all file extensions using `fs::path_ext_remove()` repeatedly.
#'
#' This function runs [fs::path_ext_remove()] as many times as needed to remove
#' all file extensions from `path`.
#'
#' @inheritParams fs::path_ext_remove
#'
#' @returns `path` with all extensions removed.
#' @export
#'
#' @examples
#' fs::path_ext_remove("a.csv.gz")
#' path_ext_remove_all("a.csv.gz")
path_ext_remove_all <- function(path) {
    req_ns("fs")
    new_path <- path_ext_remove(path)
    while (!all(na.omit(new_path == path))) {
        path <- new_path
        new_path <- path_ext_remove(path)
    }
    return(new_path)
}
