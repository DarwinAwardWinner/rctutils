#' Determine the column name of the p-value column in a table
#'
#' limma uses "P.Value", edgeR uses "PValue", so we need an
#' abstraction.
#'
#' @export
get_pval_colname <- function(ttab) {
    if (is.character(ttab)) {
        cnames <- ttab
    } else {
        cnames <- colnames(ttab)
    }
    pcol <- match(c("p.value", "pvalue", "pval", "p"),
        tolower(cnames)) %>%
        na.omit %>% .[1]
    pcolname <- cnames[pcol]
    if (length(pcolname) != 1)
        stop("Could not determine p-value column name")
    return(pcolname)
}

#' Add a q-value column to any table with a p-value column
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
    }, error=function(e) {
        warning(str_c("Failed to compute q-values: ", e$message))
    })
    ttab
}

#' Compute Bayesian FDR values from limma's B statistics
#'
#' @importFrom dplyr cummean
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
#' @export
add.bfdr <- function(ttab) {
    B <- ttab[["B"]]
    if (is.null(B)) {
        warning("Cannot add BFDR to table with no B statistics")
        return(ttab)
    }
    btab <- bfdr(B)[c("PP", "BayesFDR")]
    for (i in names(btab)) {
        ttab[[i]] <- btab[[i]]
    }
    ttab
}

#' Variant of propTrueNullByLocalFDR that does not restrict to the range [0,1]
#'
#' Useful for identifying potential atypical p-value distributions
#' (e.g. too many large p-values).
#'
#' @export
propTrueNullByLocalFDR_unrestricted <- function (p) {
    n <- length(p)
    i <- n:1L
    p <- sort(p, decreasing = TRUE)
    q <- n/i * p
    n1 <- n + 1L
    sum(i * q)/n/n1 * 2
}
