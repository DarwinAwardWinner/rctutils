#' Shortcut for the usual "requireNamespace" dance.
#'
#' This shortens the boilerplate required to use functios from a
#' suggested package.
#'
#' @importFrom stringr str_c
req_ns <- function(..., caller=as.character(sys.call(-1)[1])) {
    for (pkg in unlist(list(...))) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
        if (length(caller) != 1 || is.na(caller)) {
            caller <- "this function"
        } else {
            caller <- str_c(caller, "()")
        }
        stop(sprintf("Package '%s' must be installed to use %s", pkg, caller), call. = FALSE)
        }
    }
}
