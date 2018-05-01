#' Shortcut for the usual "requireNamespace" dance
#'
#' This shortens the boilerplate required to use functios from a
#' suggested package.
#'
#' @param ... Packages to require
#' @param caller The name of the function this was called from. Only
#'     used to create an appropriate error message.
req_ns <- function(..., caller = as.character(sys.call(-1)[1])) {
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

## Tell "R CMD check" not to worry about the magrittr pronoun
utils::globalVariables(".")

## Tell check not to worry about rex internal functions
globalVariables(c("one_or_more", "space", "zero_or_more", "capture", "maybe", "digit", "%if_prev_is%", "%if_next_isnt%", "or"))

#' Common imports
#'
#' This is just here to tell roxygen2 about all base package imports,
#' which were recommended by R CMD check, as well as some common
#' imports that are used in many functions. Adding these to every
#' individual function that uses these common functions is too
#' tedious, so I've just added them all here.
#'
#' @importFrom grDevices cairo_pdf dev.cur dev.list dev.off dev.set png
#' @importFrom graphics abline barplot lines par title
#' @importFrom methods as is new
#' @importFrom stats approx approxfun as.dist as.formula cmdscale end lowess model.matrix na.omit start
#' @importFrom utils read.csv read.table write.table
#' @import magrittr
#' @import dplyr
#' @import stringr
#' @import ggplot2
#' @importFrom assertthat assert_that
NULL
