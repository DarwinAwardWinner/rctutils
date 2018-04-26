#' Parallelized version of [limma::selectModel()]
#'
#' @include internal.R
#' @export
BPselectModel <- function (y, designlist, criterion = "aic", df.prior = 0, s2.prior = NULL,
                           s2.true = NULL, ..., BPPARAM)
{
    req_ns("BiocParallel", "limma")
    if (missing(BPPARAM)) {
        BPPARAM <- BiocParallel::bpparam()
    }
    ym <- as.matrix(y)
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
    criterion <- match_arg(criterion, c("aic", "bic", "mallowscp"))
    if (criterion == "mallowscp") {
        if (is.null(s2.true))
            stop("Need s2.true values")
        fits <- BiocParallel::bplapply(designlist, limma::lmFit, object=y, BPPARAM=BPPARAM)
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
        fits <- BiocParallel::bplapply(designlist, limma::lmFit, object=y, BPPARAM=BPPARAM)
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

#' Combine multiple featureCounts results
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

#' @include internal.R
featureCountsQuiet <- function(...) {
    req_ns("withr", "Rsubread")
    withr::with_output_sink("/dev/null", Rsubread::featureCounts(...))
}

#' Parallel version of [Rsubread::featureCounts()]
#'
#' @include internal.R
#' @export
featureCountsParallel <- function(files, ...) {
    req_ns("Rsubread", "BiocParallel")
    # Let featureCounts handle the degenerate case itself
    if (length(files) == 0) {
        return(Rsubread::featureCounts(files, ...))
    }
    ## We pass nthreads=1 to tell featureCounts not to parallelize
    ## itself, since we are already handling parallelization.
    BiocParallel::bplapply(files, featureCountsQuiet, ..., nthreads=1) %>%
        combineFCResults
}

#' Parallel version of [csaw::windowCounts()]
#'
#' @include internal.R
#' @export
windowCountsParallel <- function(bam.files, ..., filter=10, BPPARAM) {
    req_ns("BiocParallel", "csaw", "BiocGenerics")
    if (missing(BPPARAM)) {
        BPPARAM <- BiocParallel::bpparam()
    }
    reslist <- BiocParallel::bplapply(X=bam.files, FUN=csaw::windowCounts, ..., filter=0, BPPARAM=BPPARAM)
    res <- do.call(cbind, reslist)
    rm(reslist)
    keep <- rowSums(assay(res)) >= filter
    res[keep,]
}

#' Like lapply, but returns future objects.
#'
#' The results can be fetched all at once with [future::values()].
#'
#' @include internal.R
#' @export
future.lapply <- function(X, FUN, ...) {
    req_ns("future")
    FUTUREFUN <- function(x) future::future(FUN(x, ...))
    lapply(X, FUTUREFUN, ...)
}