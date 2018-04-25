add_numbered_colnames <- function(x, prefix="C") {
    x %>% set_colnames(glue("{prefix}{num}", num=seq(from=1, length.out=ncol(x))))
}

# For each column of a data frame, if it is a character vector with at least one
# repeated value, convert that column to a factor
auto_factorize <- function(df) {
    for (i in colnames(df)) {
        if (is.character(df[[i]]) && anyDuplicated(df[[i]])) {
            df[[i]] %<>% factor
        }
    }
    df
}

cleanup_mcols <- function(object, mcols_df=mcols(object)) {
    nonempty <- !sapply(mcols_df, is.empty)
    mcols_df %<>% .[nonempty]
    if (!missing(object)) {
        mcols(object) <- mcols_df
        return(object)
    } else {
        return(mcols_df)
    }
}

## TODO: Figure out how to correctly reference internal functions from
## another pkg
library(codingMatrices)
# Variant of code_control that generates more verbose column names
code_control_named <- function (n, contrasts = TRUE, sparse = FALSE) {
    if (is.numeric(n) && length(n) == 1L) {
        if (n > 1L)
            levels <- .zf(seq_len(n))
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
            B <- .asSparse(B)
        return(B)
    }
    B <- B - 1/n
    B <- B[, -1, drop=FALSE]
    colnames(B) <- paste(levels[1], levels[-1], sep = ".vs.")
    if (sparse) {
        .asSparse(B)
    }
    else {
        B
    }
}
environment(code_control_named) <- new.env(parent = environment(code_control))

# Helper function for ensureAtomicColumns
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

# Function to ensure that all columns of a data frame are atomic
# vectors. Columns that are lists have each of their elements
# collapsed into strings using the sepcified separator. TODO: If sep
# is NA, throw an error for lists.
ensureAtomicColumns <- function(df, sep=",") {

    df[] %<>% lapply(collapseToAtomic, sep=sep)
    df
}

fac2char <- function(df) {
    df[sapply(df, class) == "factor"] %<>% lapply(as.character)
    df
}


#' @include si2f.R
parse_bp <- function(size) {
    suppressWarnings({
        result <- si2f(size, "bp")
        ## Fall back to just parsing a number without the "bp" suffix
        result[is.na(result)] <- si2f(size[is.na(result)])
    })
    assert_that(!any(is.na(result)))
    result
}

#' @importFrom sitools f2si
format_bp <- function(x) {
    x %>% round %>% f2si("bp") %>% str_replace_all(rex(one_or_more(space)), "")
}

## TODO: Find a better home for this?
print_var_vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep="")
    }
    invisible(v)
}


## Given a GRangesList whose underlying ranges have mcols, find mcols
## of the ranges that are constant within each gene and promote them
## to mcols of the GRangesList. For example, if exons are annotated with
promote.common.mcols <- function(grl, delete.from.source=FALSE, ...) {
    colnames.to.promote <- get.gene.common.colnames(mcols(unlist(grl)), rep(names(grl), lengths(grl)), ...)
    promoted.df <- mcols(unlist(grl))[cumsum(lengths(grl)),colnames.to.promote, drop=FALSE]
    if (delete.from.source) {
        mcols(grl@unlistData) %<>% .[setdiff(names(.), colnames.to.promote)]
    }
    mcols(grl) %<>% cbind(promoted.df)
    grl
}

## Given a GRangesList whose underlying ranges have mcols, find mcols
## of the ranges that are constant within each gene and promote them
## to mcols of the GRangesList. For example, if exons are annotated with
promote.common.mcols <- function(grl, delete.from.source=FALSE, ...) {
    colnames.to.promote <- get.gene.common.colnames(mcols(unlist(grl)), rep(names(grl), lengths(grl)), ...)
    promoted.df <- mcols(unlist(grl))[cumsum(lengths(grl)),colnames.to.promote, drop=FALSE]
    if (delete.from.source) {
        mcols(grl@unlistData) %<>% .[setdiff(names(.), colnames.to.promote)]
    }
    mcols(grl) %<>% cbind(promoted.df)
    grl
}

#' @importFrom dplyr mutate
mutate_if_present <- function(.data, names, ...) {
    if (all(names %in% base::names(.data))) {
        mutate(.data, !!!quos(...))
    } else {
        .data
    }
}

#' @importFrom stringr str_replace_all
quotemeta <- function (string) {
  str_replace_all(string, "(\\W)", "\\\\\\1")
}

relevel_columns <- function(df, ...) {
    relevel_specs <- list(...)
    assert_that(is_named(relevel_specs))
    for (i in names(relevel_specs)) {
        df[[i]] %<>% fct_relevel(relevel_specs[[i]])
    }
    df
}

sprintf.single.value <- function(fmt, value) {
    ## Max function arguments is 100
    arglist = c(list(fmt=fmt), rep(list(value), 99))
    do.call(sprintf, arglist)
}

is_empty <- function(x) {
    x %>% unlist %>% na.omit %>% length %>% equals(0)
}

is_fully_named <- function(obj, namefun=names) {
    the.names <- namefun(obj)
    if (is.null(the.names)) {
        return(FALSE)
    } else if (!is.character(the.names)) {
        return(FALSE)
    } else if (any(is.na(the.names))) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

strip_design_factor_names <- function(design, prefixes=names(attr(design, "contrasts"))) {
    if (!is.null(prefixes)) {
        regex.to.delete <- rex(or(prefixes) %if_prev_is% or(start, ":") %if_next_isnt% end)
        colnames(design) %<>% str_replace_all(regex.to.delete, "")
        if (anyDuplicated(colnames(design))) {
            warning("Stripping design names resulted in non-unique names")
        }
    }
    design
}
