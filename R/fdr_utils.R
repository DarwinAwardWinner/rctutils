# limma uses "P.Value", edgeR uses "PValue", so we need an abstraction
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

# Add a q-value column to any table with a p-value column
add.qvalue <- function(ttab, ...) {
    tryCatch({
        P <- ttab[[get_pval_colname(ttab)]]
        qobj <- qvalue(P, ...)
        ttab$QValue <- qobj$qvalues
        ttab$LocFDR <- qobj$lfdr
        attr(ttab, "qvalue") <- qobj
    }, error=function(e) {
        warning(str_c("Failed to compute q-values: ", e$message))
    })
    ttab
}

# Compute posterior probabilities and Bayesian FDR values from limma's B
# statistics
bfdr <- function(B) {
    o <- order(B, decreasing = TRUE)
    ro <- order(o)
    B <- B[o]
    ## TODO: Finish this function?
    positive <- which(B > 0)
    PP <- exp(B)/(1+exp(B))
    # Computing from 1-PP gives better numerical precision for the
    # most significant genes (large B-values)
    oneMinusPP <- 1/(1+exp(B))
    BayesFDR <- cummean(oneMinusPP)
    data.frame(B, PP, BayesFDR)[ro,]
}

# Add Bayesian FDR to a limma top table
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

# Variant that does not restrict to the range [0,1]. Useful for identifying
# potential atypical p-value distributions (too many large p-values).
propTrueNullByLocalFDR_unrestricted <- function (p) {
    n <- length(p)
    i <- n:1L
    p <- sort(p, decreasing = TRUE)
    q <- n/i * p
    n1 <- n + 1L
    sum(i * q)/n/n1 * 2
}

# Add corrected p-value, q-value, and locfdr columns to any table with a p-value
# column. Works by the possibly questionable method of converting p-values to
# equivalent z-scores and running fdrtool on the z-scores.
add.fdrtool <- function(ttab, verbose=FALSE, plot=TRUE, convert.to.zscores, cutoff.method,
 ...) {
    P <- ttab[[get_pval_colname(ttab)]]
    assert_that(!is.null(P), length(P) > 0, !any(is.na(P)))
    if (missing(convert.to.zscores)) {
        # If pval dist is high-biased, use normal modelling instead.
        ptn <- propTrueNullByLocalFDR_unrestricted(P)
        convert.to.zscores <- ptn > 1
    }
    if (convert.to.zscores) {
        # Try normal modelling instead
        Zscore <- qnorm(1-(P/2))
        if (missing(cutoff.method)) {
            co <- fndr.cutoff(Zscore, statistic="normal")
            if (co >= 0.3) {
                cutoff.method <- "fndr"
            } else {
                cutoff.method <- "locfdr"
            }
        }
        fdrmod <- fdrtool(Zscore, statistic="normal", verbose=verbose,
            plot=plot, cutoff.method=cutoff.method, ...)

    } else {
        if (missing(cutoff.method)) {
            cutoff.method <- "fndr"
        }
        fdrmod <- fdrtool(P, statistic="pvalue", verbose=verbose,
            plot=plot, cutoff.method=cutoff.method,...)
    }
    fdrdf <- do.call(data.frame, fdrmod[c("pval", "qval", "lfdr")])
    for (i in names(fdrdf)) {
        ttab[[str_c("fdrtool.", i)]] <- fdrdf[[i]]
    }
    attr(ttab, "fdrtool") <- fdrmod
    ttab
}

# Variant of eBayes that uses propTrueNull (or a custom function) to
# set the proportion argument automatically
eBayes_autoprop <- function(..., prop.method="lfdr") {
    eb <- eBayes(...)
    if (is.function(prop.method)) {
        ptn <- prop.method(eb$p.value)
    } else {
        ptn <- propTrueNull(eb$p.value, method=prop.method)
    }
    eBayes(..., proportion=1-ptn)
}
