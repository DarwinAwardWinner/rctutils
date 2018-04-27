#' Determine whether a list or vector has any non-missing values.
#'
#' A vector is "valueless" if it has zero length or contains only NA
#' values. A list is empty if has zero length or contains only empty
#' vectors.
#'
#' @param x The object to check for emptiness.
#' @param empty_values A vector of values that should be considered
#'     "empty", in addition to NA. Any vector containing only NA and
#'     these values will be considered "valueless". For example, if
#'     you want to know whether a numeric vector contains any non-NA,
#'     non-zero values, you could use `empty_values = 0`. Another
#'     common value to use is `""`, the empty string.
#'
#' @seealso [rlang::is_empty()], which only checks for zero length.
#' @export
is_valueless <- function(x, empty_values = NA) {
    x %>% unlist %>% na.omit %>% length %>% equals(0)
}

#' Add numbered colnames with a common prefix.
#'
#' @param x The variable to add column names to.
#' @param prefix The prefix to use for each column name. Each column's
#'     number will be appended to this prefix to generate the column
#'     name.
#'
#' @return `x`, with its column names set to the column number
#'     appended to `prefix`.
#'
#' @examples
#'
#' # TODO: PCA example
#'
#' @importFrom glue glue
#' @export
add_numbered_colnames <- function(x, prefix="C") {
    x %>% set_colnames(glue("{prefix}{num}", num=seq(from=1, length.out=ncol(x))))
}

#' Intelligently convert character columns to factors in a data frame
#'
#' For each column of a data frame, if it is a character vector with
#' at least one non-unique value, convert that column to a factor.
#' Character columns with no repeated values (e.g. sample IDs) will not
#' be modified.
#'
#' @param df The data frame to operate on.
#'
#' @return `df`, possibly with some character columns replaced by
#'     factors
#'
#' @examples
#'
#' # Initially, both columns are character vectors
#' x <- data.frame(sample=letters, group=LETTERS[1:2], stringsAsFactors=FALSE)
#' sapply(x, class)
#'
#' # This converts group, but not sample, into a factor
#' x2 <- auto_factorize_columns(x)
#' sapply(x2, class)
#'
#' @export
auto_factorize_columns <- function(df) {
    for (i in colnames(df)) {
        if (is.character(df[[i]]) && anyDuplicated(df[[i]])) {
            df[[i]] <- factor(df[[i]])
        }
    }
    df
}

#' Remove valueless mcols from an object
#'
#' Any mcols for which [is_valueless()] returns TRUE will be removed.
#'
#' @param object The object whose mcols should be cleaned.
#' @param mcols_df Instead of passing `object`, you can pass a data
#'     frame (or [S4Vectors::DataFrame]) using this named argument.
#'
#' @return If `object` was passed, it is returned with any valueless
#'     mcols removed. If `mcols_df` was passed and `object` was not,
#'     the data frame is returned directly.
#'
#' Note that `mcols_df` can only be passed by name. An unnamed
#' argument is always interpreted as `object`.
#'
#' @examples
#' #TODO
#'
#' @seealso [S4Vectors::mcols()], [is_valueless()]
#' @export
cleanup_mcols <- function(object, mcols_df = S4Vectors::mcols(object)) {
    req_ns("S4Vectors")
    nonempty <- !sapply(mcols_df, is_valueless)
    mcols_df <- mcols_df[nonempty]
    if (!missing(object)) {
        S4Vectors::mcols(object) <- mcols_df
        return(object)
    } else {
        return(mcols_df)
    }
}

#' Variant of code_control that generates more verbose column names
#'
#' @param ... All arguments are passed to
#'     `codingMatrices::code_control()`.
#'
#' @return The same as `codingMatrices::code_control()`, but with more
#'     descriptive column names.
#'
#' @seealso [codingMatrices::code_control()]
#' @examples
#' library(codingMatrices)
#'
#' # Compare the resulting column names
#' code_control(3)
#' code_control_named(3)
#'
#' @export
code_control_named <- function (...) {
    req_ns("codingMatrices")
    x <- codingMatrices::code_control(...)
    colnames(x) %<>% str_replace("(\\d+)-(\\d+)", "\\2.vs.\\1")
    x
}

#' Convert a list to an atomic vector.
#'
#' If `x` is already an atomc vector, it is not modified. If it is a
#' list, non-scalar elements are converted to scalars, and these are
#' concatenated to produce an atomic vector the same length as `x`.
#' Specifically, list elements with 0 length are replaced by NA, and
#' elements with length greater than one are pasted together into
#' strings using `sep` as the separator.
#'
#' @param x The list to convert.
#' @param sep Separator to use for collapsing vectors to scalars.
#'
#' @return An atomic vector that same length and names as `x`.
#'
#'
#' @export
collapse_to_atomic <- function(x, sep=",") {
    if (is.atomic(x)) {
        return(x)
    } else {
        y <- lapply(x, str_c, collapse=sep)
        y[lengths(y) == 0] <- NA
        assert_that(all(lengths(y) == 1))
        ## Ensure that unlisting does not mess up the names
        y <- set_names(unlist(y), names(y))
        assert_that(length(y) == length(x))
        return(y)
    }
}

#' Ensure that all columns of a data frame are atomic vectors.
#'
#' Columns that are lists have each of their elements collapsed into
#' strings using the sepcified separator.
#'
#' @param df The data frame whose list columns should be collapsed
#' @param sep The separator to use when collapsing vectors to scalars
#'
#' @return `df` with all columns converted to atomic vectors.
#'
#' @examples
#'
#' library(S4Vectors)
#' x <- DataFrame(a=1:5, b=List(split(1:25, 1:5)))
#' ensure_atomic_columns(x)
#'
#' @seealso [collapse_to_atomic()]
#'
#' @export
ensure_atomic_columns <- function(df, sep=",") {
    df[] %<>% lapply(collapse_to_atomic, sep=sep)
    df
}

#' Convert all factors in a data frame to character vectors
#'
#' This is useful when you forgot to specify `stringsAsFactors =
#' FALSE` when creating your data frame.
#'
#' @param df The data frame in which fo convert factors.
#' @return `df`, with all factor columns converted to character
#'     vectors.
#'
#' @examples
#'
#' x <- data.frame(letters, 1:26, LETTERS)
#' sapply(x, class)
#' x2 <- fac2char(x)
#' sapply(x2, class)
#'
#' @seealso [as.character()], [data.frame()]
#' @export
fac2char <- function(df) {
    df[sapply(df, class) == "factor"] %<>% lapply(as.character)
    df
}

#' Parse a number of base pairs with optional units.
#'
#' @examples
#' parse_bp(c("100", "100bp", "100 kbp", "1Mbp"))
#'
#' @export
parse_bp <- function(size) {
    suppressWarnings({
        result <- si2f(size, "bp")
        ## Fall back to just parsing a number without the "bp" suffix
        ## (which also includes bare numbers with no unit at all)
        result[is.na(result)] <- si2f(size[is.na(result)])
    })
    ## Require everyting to parse successfully or throw an error
    assert_that(!any(is.na(result)))
    result
}

#' Format a number of base pairs using the most appropriate unit.
#'
#' @param x A numeric vector representing numbers of base pairs.
#'
#' @return A character vector containing string representations of the
#'     number of base pairs, using appropriate SI prefixes.
#'
#' @examples
#'
#' format_bp(c(100, 1000, 1e6))
#'
#' @export
format_bp <- function(x) {
    req_ns("sitools", "rex")
    x %>% round %>% sitools::f2si("bp") %>% str_replace_all(rex::rex(one_or_more(space)), "")
}

#' Perform mutations only if specific column names are present
#'
#' Note that columns required to be present need not be involved in the mutations.
#'
#' @param .data,... These have the same meaning as in
#'     [dplyr::mutate()]
#' @param names If `.data` does not contain all of these names as
#'     columns, it will be returned unmodified.
#'
#' @seealso [dplyr::mutate()]. Not to be confused with
#'     [dplyr::mutate_if()], which mutates columns that match a
#'     specific filter criteria.
#'
#' @export
mutate_if_present <- function(.data, names, ...) {
    if (all(names %in% base::names(.data))) {
        mutate(.data, !!!quos(...))
    } else {
        .data
    }
}

#' Escape a string for use in a Perl-compatible regular expression
#'
#' This is an implementation of Perl's `quotemeta` operator. It is
#' meant for escaping any characters that might have special meanings
#' in a PCRE.
#'
#' @param string The string to escape
#' @return `string`, with all non-word characters escaped with backslashes.
#'
#' Note that this strategy does not work for R's own regular
#' expression syntax, so you should specify `perl = TRUE`, or an
#' equivalent option, when using this to construct a regular
#' expression.
#'
#' @examples
#'
#' s <- "A (string) with [many] regex *characters* in it."
#' grepl(s, s, perl=TRUE)
#' quotemeta(s)
#' # After escaping, the string matches itself.
#' grepl(quotemeta(s), s, perl=TRUE)
#'
#' @export
quotemeta <- function (string) {
  str_replace_all(string, "(\\W)", "\\\\\\1")
}

#' Relevel many factor columns in a data frame at once
#'
#' This is a shortcut for calling [forcats::fct_relevel()] on multiple
#' columns of a data frame.
#'
#' @param df The data frame in which to re-level some factor columns.
#' @param ... Additional arguments control the re-leveling of factors.
#'     The name of each argument indicates the column of `df` to
#'     re-level (this column should be a factor). The value should be
#'     the new levels in the order that they should appear (or any
#'     valid list of arguments to [forcats::fct_relevel()]).
#' @return `df`, with the specified factor columns re-leveled
#'     according to the provided specifications.
#' @examples
#'
#' x <- data.frame(a=letters[1:2], b=LETTERS[1:6])
#' sapply(x, levels)
#' x2 <- relevel_columns(x, a=c("b", "a"), b=list("B", "A", after=2))
#' sapply(x2, levels)
#'
#' @importFrom rlang is_named
#' @export
relevel_columns <- function(df, ...) {
    req_ns("forcats", "rlang")
    relevel_specs <- list(...)
    assert_that(is_named(relevel_specs))
    for (i in names(relevel_specs)) {
        spec <- relevel_specs[[i]]
        after <- 0L
        if (!is.list(spec)) {
            spec <- list(spec)
        }
        arglist <- c(list(.f=df[[i]]), spec)
        df[[i]] <- do.call(forcats::fct_relevel, arglist)
    }
    df
}

#' Like `sprintf()` but using the same value in each specification
#'
#' This is equivalent to calling [sprintf()] with `value` passed after
#' `fmt` as many times as there are conversion specifications in
#' `fmt`.
#'
#' @param fmt This has the same meaning as in `sprintf()`.
#' @param value The single value use for every conversion
#'     specification in `fmt`.
#' @return See [sprintf()].
#'
#' Note that vector arguments are still recycled normally.
#'
#' @examples
#'
#' fmt <- "Hello %s. Your name is %s. Isn't it good to be named %s?"
#' # With normal sprintf, you need to specify the right number of
#' # arguments:
#' name <- c("John Doe", "Jane Doe")
#' sprintf(fmt, name, name, name)
#' # Instead, this figures out the right number of times to pass it:
#' sprintf_single_value(fmt, name)
#'
#' @export
sprintf_single_value <- function(fmt, value) {
    ## Max function arguments is 100
    arglist = c(list(fmt=fmt), rep(list(value), 99))
    do.call(sprintf, arglist)
}

#' Simplify the column names of a simple design matrix
#'
#' The standard [model.matrix()] function generates column names for
#' factor columns by concatenating the factor's name with the level.
#' This helps guarantee non-duplicated column names in the general
#' case, but for simple cases it results in overly long and complex
#' column names. This function removes the factor names, leaving only
#' the level names.
#'
#' @param design The design matrix whose column names should be
#'     cleaned up (usually created by `model.matrix()`).
#' @param prefixes A character vector of prefixes to be deleted from
#'     column names. Normally, this argument is set automatically
#'     based on the attributes of `design`, but you can set it
#'     manually.
#' @return The same matrix as `design`, possibly with some column
#'     names shortened by deleting one of `prefixes` from the start of
#'     the column names.
#'
#' Column names for columns not corresponding to factors are left
#' unchanged.
#'
#' Note that this is only guaranteed to work correctly if none of your
#' variable names are prefixes of other variable names.
#'
#' A warning will be issued if, after stripping factor names, the
#' column names of `design` are no longer unique.
#'
#' @examples
#'
#' targets <- data.frame(group=letters[1:2], batch=LETTERS[1:5], x=rnorm(10))
#' design <- model.matrix(~0 + group + x + batch, targets)
#' colnames(design)
#' design2 <- strip_design_factor_names(design)
#' colnames(design2)
#'
#' @export
strip_design_factor_names <- function(design, prefixes=names(attr(design, "contrasts"))) {
    req_ns("rex")
    if (!is.null(prefixes)) {
        regex.to.delete <- rex::rex(or(prefixes) %if_prev_is% or(start, ":") %if_next_isnt% end)
        colnames(design) %<>% str_replace_all(regex.to.delete, "")
        if (anyDuplicated(colnames(design))) {
            warning("Stripping design names resulted in non-unique names")
        }
    }
    design
}
