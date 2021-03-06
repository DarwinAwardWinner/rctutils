#' Determine the column name of the p-value column in a table
#'
#' limma uses "P.Value" while edgeR uses "PValue", so we need an
#' abstraction.
#'
#' @param ttab A table of results produced by [limma::topTable()],
#'     [edgeR::topTags()], [DESeq2::results()], or any table with a
#'     p-value column in it.
#'
#' @return The name of the column containing p-values. If no such
#'     column is found, an error is thrown.
#'
#' @export
get_pval_colname <- function(ttab) {
    if (is.character(ttab)) {
        cnames <- ttab
    } else {
        cnames <- colnames(ttab)
    }
    ## Match case-insensitively, but return the name in its original
    ## case.
    pcol <- match(c("p.value", "pvalue", "pval", "p", "p value"),
        tolower(cnames)) %>%
        na.omit %>% .[1]
    pcolname <- cnames[pcol]
    if (length(pcolname) != 1)
        stop("Could not determine p-value column name")
    return(pcolname)
}

#' Add a q-value column to any table with a p-value column
#'
#' This function takes a results table and augments it with columns
#' containing q-values and local FDR values computed using
#' [qvalue::qvalue()]
#'
#' If a p-value column cannot be found or if the q-value computation
#' fails, `ttab` will be returned with no new columns added, with a
#' warning. This might happen if the p-value distribution is highly
#' atypical, such that it is impossible to meaningfully estimate
#' `pi0`.
#' @param ttab A table of results produced by [limma::topTable()],
#'     [edgeR::topTags()], [DESeq2::results()], or any table with a
#'     p-value column in it.
#' @param ... Additional arguments are passed to [qvalue::qvalue()].
#' @return `ttab` possibly with two new columns added, named "QValue"
#'     and "LocFDR".
#'
#' @examples
#'
#' #TODO steal from topTable
#'
#' @export
add_qvalue <- function(ttab, ...) {
    req_ns("qvalue")
    tryCatch({
        P <- ttab[[get_pval_colname(ttab)]]
        qobj <- qvalue::qvalue(P, ...)
        ttab$QValue <- qobj$qvalues
        ttab$LocFDR <- qobj$lfdr
        attr(ttab, "qvalue") <- qobj
    }, error = function(e) {
        warning(str_c("Failed to compute q-values: ", e$message))
    })
    ttab
}

#' Compute Bayesian FDR values from limma's B-statistics
#'
#' @param B A vector of B-statistics
#' @return A data frame containing 3 columns: "B", the input; "PP",
#'     the posterior probability of each gene being differentially
#'     expressed; and "BFDR", the estimated FDR.
#'
#' @seealso [limma::topTable()].
#'
#' @examples
#'
#' # TODO Merge this documentation with that of add_bfdr().
#'
#' @export
bfdr <- function(B) {
    o <- order(B, decreasing = TRUE)
    ro <- order(o)
    B <- B[o]
    PP <- exp(B)/(1+exp(B))
    # Computing FDR from from 1-PP gives better numerical precision
    # for the most significant genes (large B-values)
    oneMinusPP <- 1/(1+exp(B))
    BayesFDR <- cummean(oneMinusPP)
    data.frame(B, PP, BayesFDR)[ro,]
}

#' Add Bayesian FDR values to a limma top table
#'
#' This function uses [bfdr()] to compute the posterior probability of
#' differential expression and the FDR based on limma's B-statistics
#' and add them as columns named "PP" and "BFDR" respectively.
#'
#' @param ttab A data frame returned by [limma::topTable()] with
#'     `n=Inf` containing B-statistics.
#' @return The same table with additional columns named "PP" and
#'     "BFDR".
#'
#' Note that limma only produces B-statistics for single-contrast
#' tests (i.e. moderated t-tests). As a convenience, if this fuction
#' is called on a data frame with no column named "B", it will issue a
#' warning and return the data frame unmodified.
#'
#' @examples
#'
#' # TODO copy examples from limma::topTable()
#'
#' @seealso [limma::topTable()].
#'
#' @export
add_bfdr <- function(ttab) {
    B <- ttab[["B"]]
    if (is.null(B)) {
        warning("Cannot add BFDR to table with no B-statistics")
        return(ttab)
    }
    btab <- bfdr(B)[c("PP", "BayesFDR")]
    for (i in names(btab)) {
        ttab[[i]] <- btab[[i]]
    }
    ttab
}
