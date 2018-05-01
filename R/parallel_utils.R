## https://github.com/HenrikBengtsson/future/issues/162
## https://github.com/DarwinAwardWinner/future/blob/6a000af1e9ea41674c85a5476cf5e8c6c9e75d80/R/use_futures.R

#' Set up foreach and BiocParallel to use futures
#'
#' This function ensures that all parallel functions of the
#' BiocParallel and foreach packages will use the selected future.
#' This includes functions like [foreach::foreach()] with `%dopar%`,
#' [BiocParallel::bplapply()], and [[plyr::llply()]].
#'
#' @param strategy,... If provided, these are passed directly to
#'     [future::plan()]. If `strategy` is not provided,
#'     `future::plan()` is not called at all, leaving the current
#'     execution strategy unchanged and ignoring any additional
#'     arguments. This is useful if you have already called
#'     `future::plan()` and you just want to set up other parallel
#'     packages to use futures.
#' @param quiet If FALSE (the default), indicate what is being done.
#'     If TRUE, do not issue any messages.
#'
#' Note that this forces loading of the BiocParallel and foreach
#' packages. Ideally, this function would only set up a hook to run
#' the appropriate setup code after these packages are loaded, thus
#' saving the setup time if they are never loaded. However, this does
#' not appear to be possible. The available "on load" hook mechanism
#' seems to only trigger when a package is attached, not when it is
#' loaded. Hence, if another pacakge uses BiocParallel or foreach
#' internally, (e.g. [SummarizedExperiment::summarizeOverlaps()]),
#' this would not trigger the hook.
#'
#' @examples
#'
#' # Set up BiocParallel and foreach to use futures, keeping the
#' # existing future execution strategy.
#' use_futures()
#'
#' # Same, but switch to the multicore strategy.
#' use_futures("multisession")
#'
#' # Same, but pass additional options to the multicore strategy.
#' use_futures("multisession", workers=2)
#'
#' @export
use_futures <- function(strategy, ..., quiet = FALSE) {
    req_ns("future")
    if (quiet) {
        message <- identity
    }
    if (!missing(strategy)) {
        message("Setting up new future execution strategy.")
        future::plan(strategy, ...)
    } else {
        message("Using existing future execution strategy.")
    }
    use_futures_for_foreach(quiet = quiet)
    use_futures_for_BiocParallel(quiet = quiet)
}

use_futures_for_foreach <- function(quiet = FALSE) {
    if (quiet) {
        message <- identity
    }
    if (requireNamespace("foreach", quietly = TRUE)) {
        if (requireNamespace("doFuture", quietly=TRUE)) {
            doFuture::registerDoFuture()
            message("Foreach will now use the doFuture backend.")
        } else {
            warning("Install the doFuture package to allow foreach to use futures for parallel operation.")
        }
    } else {
        message("Not setting up foreach to use futures because it is not installed.")
    }
    NULL
}

use_futures_for_BiocParallel <- function(quiet = FALSE, via_foreach=FALSE) {
    if (!via_foreach && !requireNamespace("BiocParallel.FutureParam", quietly = TRUE)) {

    }
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
        if (via_foreach) {
            suppressMessages(use_futures_for_foreach(quiet = TRUE))
            BiocParallel::register(BiocParallel::DoparParam())
            message("BiocParallel will now use the DoparParam (i.e. foreach) backend, which should in turn use the doFuture backend.")
        } else if (requireNamespace("BiocParallel.FutureParam", quietly = TRUE)) {
            BiocParallel::register(BiocParallel.FutureParam::FutureParam())
            message("BiocParallel will now use the FutureParam backend")
        } else {
            warning("Install the BiocParallel.FutureParam package to allow BiocParallel to use futures for parallel operation.")
        }
    } else {
        message("Not setting up BiocParallel to use futures because it is not installed.")
    }
}

#' Parallelized version of [limma::selectModel()]
#'
#' @param y,designlist,criterion,df.prior,s2.prior,s2.true,... These
#'     arguments all have the same meaning as in
#'     [limma::selectModel()].
#' @param BPPARAM A BiocParallelParam instance that determines how to
#'     parallelize the operation
#'
#' @return See [limma::selectModel()].
#'
#' @seealso [limma::selectModel()], [BiocParallel::bplapply()]
#'
#' @export
selectModelParallel <- function (y, designlist, criterion = "aic", df.prior = 0, s2.prior = NULL,
                                 s2.true = NULL, ..., BPPARAM = BiocParallel::bpparam())
{
    req_ns("BiocParallel", "limma")
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
#'
#' This is a helper function for `featureCountsParallel()`.
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

#' `Rsubread::featureCounts()` with output suppressed
#'
#' This is a helper function for `featureCountsParallel()`.
#'
#' This redirects the output to `/dev/null`, so it assumes a UNIX-like
#' system.
featureCountsQuiet <- function(...) {
    req_ns("withr", "Rsubread")
    withr::with_output_sink("/dev/null", Rsubread::featureCounts(...))
}

#' Parallel version of `Rsubread::featureCounts()`
#'
#' @param files This has the same meaning as in
#'     [Rsubread::featureCounts()]
#' @param ... Other arguments are passed to
#'     [Rsubread::featureCounts()]
#' @param BPPARAM A BiocParallelParam instance that determines how to
#'     parallelize the operation
#'
#' @return See [Rsubread::featureCounts()].
#'
#' @seealso [Rsubread::featureCounts()], [BiocParallel::bplapply()]
#'
#' @export
featureCountsParallel <- function(files, ...,
                                  BPPARAM = BiocParallel::bpparam()) {
    req_ns("Rsubread", "BiocParallel")
    # Let featureCounts handle the degenerate case itself
    if (length(files) == 0) {
        return(Rsubread::featureCounts(files, ...))
    }
    ## We pass nthreads=1 to tell featureCounts not to parallelize
    ## itself, since we are already handling parallelization.
    BiocParallel::bplapply(files, featureCountsQuiet, ..., nthreads=1, BPPARAM = BPPARAM) %>%
        combineFCResults
}

#' Parallel version of `csaw::windowCounts()`
#'
#' @param bam.files,filter These arguments have the same meaning
#'     as in [csaw::windowCounts()]
#' @param ... Other arguments are passed to [csaw::windowCounts()]
#' @param BPPARAM A BiocParallelParam instance that determines how to
#'     parallelize the operation
#'
#' @return See [csaw::windowCounts()].
#'
#' @seealso [csaw::windowCounts()], [BiocParallel::bplapply()]
#' @export
windowCountsParallel <- function(bam.files, ..., filter=10,
                                 BPPARAM = BiocParallel::bpparam()) {
    req_ns("BiocParallel", "csaw", "SummarizedExperiment")
    reslist <- BiocParallel::bplapply(X=bam.files, FUN=csaw::windowCounts, ..., filter=0, BPPARAM=BPPARAM)
    res <- do.call(cbind, reslist)
    rm(reslist)
    keep <- rowSums(SummarizedExperiment::assay(res)) >= filter
    res[keep,]
}

#' Like lapply, but returns a list of future objects.
#'
#' @param X,FUN,... See [lapply()].
#'
#' @return A list of futures. You can fetch all the results in a list
#'     using [future::values()].
#'
#' @seealso [lapply()], [future::future()]
#'
#' @export
future.lapply <- function(X, FUN, ...) {
    req_ns("future")
    FUTUREFUN <- function(x) future::future(FUN(x, ...))
    lapply(X, FUTUREFUN, ...)
}
