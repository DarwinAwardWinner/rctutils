## Fixed versions of existing functions:

#' Deparse and then concatenate into a single string
#'
#' This is like the base R `deparse()`, except that it always returns
#' a single string, rather than a character vector of lines of text.
#'
#' @param ... See [deparse()].
#' @return A character vector of length 1.
#'
#' Note that depending on the different deparse options, the resulting
#' string may still not be a complete R expression that recreates the
#' original object. Also note that the returned string may contain
#' newlines.
#'
#' @examples
#'
#' # Make an object whose representation spans multiple lines
#' x <- seq(100)[-50]
#' length(deparse(x))
#' length(deparse_onestring(x))
#' identical(x, eval(parse(text=deparse_onestring(x))))
#'
#' @export
deparse_onestring <- function(...) {
    paste(deparse(...), collapse = "\n")
}

#' Improved match.arg with better errors and case folding.
#'
#' This function works identically to [match.arg()], except that the
#' error messages include the actual name of the argument. So for
#' example, if the argument name is "size", instead of a generic error
#' message "'arg' must be of length 1", you will get
#' "'size' must be of length 1". This makes it easier for users to
#' figure out which argument has been passed incorrectly.
#'
#' @param arg_name The name of the argument. Normally you do not need
#'     to specify this yourself, as it is detected automatically.
#'     However, if your variable name is not the name of the argument,
#'     you will need to specify this.
#' @param ignore.case If TRUE, case will be ignored when comparing the
#'     argument to the choices. Note that the matching choice will be
#'     returned in its original case, regardless of the case of arg.
#' @inheritParams base::match.arg
#'
#' @examples
#'
#' choices <- c("red", "green", "blue")
#' color <- "blue"
#' match_arg(color, choices)
#' color <- "pink"
#' try(match_arg(color, choices))
#' try(match_arg(color, choices, arg_name = "colour"))
#' color <- "BLUE"
#' match_arg(color, choices, ignore.case = TRUE)
#' try(match_arg(color, choices))
#'
#' @importFrom glue glue
#' @export
match_arg <- function (arg, choices, several.ok = FALSE, arg_name = substitute(arg), ignore.case = FALSE) {
    if (missing(choices)) {
        formal.args <- formals(sys.function(sys.parent()))
        choices <- eval(formal.args[[as.character(substitute(arg))]])
    }
    if (is.name(arg_name)) {
        arg_name <- as.character(arg_name)
    }
    arg_name_string <- deparse_onestring(arg_name)
    if (is.null(arg))
        return(choices[1L])
    else if (!is.character(arg))
        stop(glue("{arg_name_string} must be NULL or a character vector"))
    if (!several.ok) {
        if (identical(arg, choices))
            return(arg[1L])
        if (length(arg) > 1L)
            stop(glue("{arg_name_string} must be of length 1"))
    }
    else if (length(arg) == 0L)
        stop(glue("{arg_name_string} must be of length >= 1"))
    fold_case <- identity
    if (ignore.case) {
        fold_case <- tolower
    }
    i <- pmatch(fold_case(arg), fold_case(choices), nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        stop(gettextf("%s should be one of %s", arg_name_string, paste(dQuote(choices),
            collapse = ", ")), domain = NA)
    i <- i[i > 0L]
    ## TODO: This code seems unreachable? Is it also unreachable in
    ## match.arg? Possible R bug?
    if (!several.ok && length(i) > 1)
        stop(glue("there is more than one match for {arg_name_string} in \"match_arg\""))
    choices[i]
}
