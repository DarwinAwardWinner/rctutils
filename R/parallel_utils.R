BPselectModel <- function (y, designlist, criterion = "aic", df.prior = 0, s2.prior = NULL,
                           s2.true = NULL, ..., BPPARAM=bpparam())
{
    ym <- as.matrix(dge)
    if (any(is.na(ym)))
        stop("NAs not allowed")
    narrays <- ncol(ym)
    rm(ym)
    nmodels <- length(designlist)
    models <- names(designlist)
    if (is.null(models))
        models <- as.character(1:nmodels)
    if (df.prior > 0 & is.null(s2.prior))
        stop("s2.prior must be set")
    if (df.prior == 0)
        s2.prior <- 0
    criterion <- match.arg(criterion, c("aic", "bic", "mallowscp"))
    criterion <- match.arg(criterion, c("aic", "bic", "mallowscp"))
    if (criterion == "mallowscp") {
        if (is.null(s2.true))
            stop("Need s2.true values")
        fits <- bplapply(designlist, lmFit, object=y, BPPARAM=BPPARAM)
        for (i in 1:nmodels) {
            fit <- fits[[i]]
            npar <- narrays - fit$df.residual[1]
            if (i == 1) {
                IC <- matrix(nrow = nrow(fit), ncol = nmodels,
                    dimnames = list(Probes = rownames(fit), Models = models))
                if (length(s2.true) != nrow(fit) && length(s2.true) != 1)
                    stop("s2.true wrong length")
            }
            IC[, i] <- fit$df.residual * fit$sigma^2/s2.true +
                npar * 2 - narrays
        }
    }
    else {
        ntotal <- df.prior + narrays
        penalty <- switch(criterion, bic = log(narrays), aic = 2)
        fits <- bplapply(designlist, lmFit, object=y, BPPARAM=BPPARAM)
        for (i in 1:nmodels) {
            fit <- fits[[i]]
            npar <- narrays - fit$df.residual[1] + 1
            s2.post <- (df.prior * s2.prior + fit$df.residual *
                            fit$sigma^2)/ntotal
            if (i == 1)
                IC <- matrix(nrow = nrow(fit), ncol = nmodels,
                    dimnames = list(Probes = rownames(fit), Models = models))
            IC[, i] <- ntotal * log(s2.post) + npar * penalty
        }
    }
    pref <- factor(apply(IC, 1, which.min), levels = 1:nmodels,
        labels = models)
    list(IC = IC, pref = pref, criterion = criterion)
}

combineFCResults <- function(fcreslist) {
    combfuncs <- list(
        counts=cbind,
        counts_junction=cbind,
        annotation=function(x, ...) x,
        targets=c,
        stat=function(...) {
            statlist <- list(...)
            firstcol <- statlist[[1]][,1, drop=FALSE]
            restcols <- statlist %>% lapply(. %>% .[,-1, drop=FALSE])
            cbind(firstcol, do.call(cbind, restcols))
        })
    res <- list()
    for (i in names(combfuncs)) {
        if (i %in% names(fcreslist[[1]])) {
            res[[i]] <- fcreslist %>%
                lapply(`[[`, i) %>%
                do.call(what=combfuncs[[i]])
        }
    }
    res
}

#' @importFrom withr with_output_sink
featureCountsQuiet <- function(...) {
    with_output_sink("/dev/null", featureCounts(...))
}

featureCountsParallel <- function(files, ...) {
    # Let featureCounts handle the degenerate case itself
    if (length(files) == 0) {
        return(featureCounts(files, ...))
    }
    bplapply(files, featureCountsQuiet, ..., nthreads=1) %>%
        combineFCResults
}

windowCountsParallel <- function(bam.files, ..., filter=10, BPPARAM=bpparam()) {
    reslist <- bplapply(X=bam.files, FUN=windowCounts, ..., filter=0, BPPARAM=BPPARAM)
    assert_that(all(sapply(reslist, is, "SummarizedExperiment")))
    res <- do.call(cbind, reslist)
    rm(reslist)
    keep <- rowSums(assay(res)) >= filter
    res[keep,]
}

# Like lapply, but returns future objects. The results can be fetched all at
# once with values().
future.lapply <- function(X, FUN, ...) {
    FUTUREFUN <- function(x) future(FUN(x, ...))
    lapply(X, FUTUREFUN, ...)
}
