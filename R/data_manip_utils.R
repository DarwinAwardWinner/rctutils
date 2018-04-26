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

#' Add numbered colnames with a commo prefix.
#'
#' @examples
#' # TODO: PCA example
#' @importFrom glue glue
#' @export
add_numbered_colnames <- function(x, prefix="C") {
    x %>% set_colnames(glue("{prefix}{num}", num=seq(from=1, length.out=ncol(x))))
}


#' Intelligently convert character columns to factors in a data frame
#'
#' For each column of a data frame, if it is a character vector with
#' at least one repeated value, convert that column to a factor.
#' Character columns with all unique values (e.g. sample IDs) will not
#' be modified.
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
#' Emptiness is defined as in [is_valueless()].
#'
#' @export
cleanup_mcols <- function(object, mcols_df) {
    req_ns("S4Vectors")
    if (missing(mcols_df)) {
        mcols_df <- S4Vectors::mcols(object)
    }
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
#' @export
code_control_named <- function (n, contrasts = TRUE, sparse = FALSE) {
    req_ns("codingMatrices")
    if (is.numeric(n) && length(n) == 1L) {
        if (n > 1L)
            levels <- codingMatrices:::.zf(seq_len(n))
        else stop("not enough degrees of freedom to define contrasts")
    }
    else {
        levels <- as.character(n)
        n <- length(n)
    }
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if (!contrasts) {
        if (sparse)
            B <- codingMatrices:::.asSparse(B)
        return(B)
    }
    B <- B - 1/n
    B <- B[, -1, drop=FALSE]
    colnames(B) <- paste(levels[1], levels[-1], sep = ".vs.")
    if (sparse) {
        codingMatrices:::.asSparse(B)
    }
    else {
        B
    }
}

#' Convert a list to an atomic vector.
#'
#' TODO: If sep is NA, throw an error for lists.
#'
#' @export
collapseToAtomic <- function(x, sep=",") {
    if (is.atomic(x)) {
        return(x)
    } else {
        y <- lapply(x, str_c, collapse=sep)
        y[lengths(y) == 0] <- NA
        assert_that(all(lengths(y) == 1))
        y <- unlist(y)
        assert_that(length(y) == length(x))
        return(y)
    }
}

#' Ensure that all columns of a data frame are atomic vectors.
#'
#' Columns that are lists have each of their elements collapsed into
#' strings using the sepcified separator. TODO: If sep is NA, throw an
#' error for lists.
#'
#' @export
ensureAtomicColumns <- function(df, sep=",") {

    df[] %<>% lapply(collapseToAtomic, sep=sep)
    df
}

#' Convert all factors to character in a data frame
#'
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
#' @export
format_bp <- function(x) {
    req_ns("sitools", "rex")
    x %>% round %>% sitools::f2si("bp") %>% str_replace_all(rex::rex(one_or_more(space)), "")
}

#' Print a readable summary of a list of values.
#'
#' This is useful for printing out a list of the parsed command-line
#' arguments for a script.
#'
#' TODO: Find a better home for this?
#'
#' @export
print_var_vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep="")
    }
    invisible(v)
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

#' @export
quotemeta <- function (string) {
  str_replace_all(string, "(\\W)", "\\\\\\1")
}

#' @importFrom rlang is_named
#' @export
relevel_columns <- function(df, ...) {
    req_ns("forcats")
    relevel_specs <- list(...)
    assert_that(is_named(relevel_specs))
    for (i in names(relevel_specs)) {
        df[[i]] %<>% forcats::fct_relevel(relevel_specs[[i]])
    }
    df
}

#' @export
sprintf.single.value <- function(fmt, value) {
    ## Max function arguments is 100
    arglist = c(list(fmt=fmt), rep(list(value), 99))
    do.call(sprintf, arglist)
}

#' @export
strip_design_factor_names <- function(design, prefixes=names(attr(design, "contrasts"))) {
    req_ns("rex")
    if (!is.null(prefixes)) {
        regex.to.delete <- res::rex(or(prefixes) %if_prev_is% or(start, ":") %if_next_isnt% end)
        colnames(design) %<>% str_replace_all(regex.to.delete, "")
        if (anyDuplicated(colnames(design))) {
            warning("Stripping design names resulted in non-unique names")
        }
    }
    design
}
