#' Shortcut for the usual "requireNamespace" dance
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

#' Emulate default arguments from the body of a function
#'
#' This is useful if you need to do a requireNamespace before you can
#' evaluate the default.
#'
#' @importFrom assertthat assert_that
#' @importFrom rlang is_missing
arg_default <- function(...) {
    uneval_arglist <- as.list(substitute(list(...)))[-1L]
    caller_env <- parent.frame()
    current_vals <- mget(names(uneval_arglist), caller_env)
    missing_args <- sapply(current_vals, is_missing)
    for (i in names(uneval_arglist)[missing_args]) {
        assign(i, value = eval(uneval_arglist[[i]]), envir = caller_env)
    }
}

## Tell "R CMD check" not to worry about the magrittr pronoun
utils::globalVariables(".")

## Tell check not to worry about rex internal functions
globalVariables(c("one_or_more", "space", "zero_or_more", "capture", "maybe", "digit", "%if_prev_is%", "%if_next_isnt%", "or"))

#' This is just here to tell roxygen2 about all base package imports,
#' which were recommended by R CMD check. Adding these to every
#' individual function that uses these common functions is too
#' tedious, so I've just added them all here.
#'
#' @importFrom grDevices dev.cur dev.list dev.off dev.set png
#' @importFrom graphics abline barplot lines par title
#' @importFrom methods as is new
#' @importFrom stats approx approxfun as.dist as.formula cmdscale end lowess model.matrix na.omit start
#' @importFrom utils read.csv read.table write.table
NULL

#' This is here to tell roxygen2 about some common imports that are
#' used in many functions. Listing them in each of the individual
#' function's roxygen2 chunks would be overly tedious.
#'
#' @importFrom magrittr %<>% %>% %$%
#' @import dplyr
#' @import stringr
#' @import ggplot2
#' @importFrom assertthat assert_that
NULL
