#' Alternative interface to `limma::removeBatchEffect()`
#'
#' Instead of accepting the batch effects and covariates separately,
#' this function accepts a single design matrix, and an index of which
#' columns of that design matrix represent batch effects.
#'
#' @param x This has the same meaning as in
#'     [limma::removeBatchEffect()].
#' @param design The full design matrix, including both covariates
#'     (whose effects are to be kept) and batch effects (whose effects
#'     are to be subtracted out).
#' @param coefsToSubtract An index vector indicating which columns of
#'     `design` should be subtracted out. This can be any vector such
#'     that `design[,coefsToSubtract]` selects the appropriate
#'     columns.
#' @param ... Additional arguments are passed to
#'     [limma::removeBatchEffect()].
#'
#' @return See [limma::removeBatchEffect()].
#'
#' @examples
#'
#' #TODO: Copy from removeBatchEffect
#'
#' @export
subtractCoefs <- function(x, design, coefsToSubtract, ...) {
    req_ns("limma")
    assert_that(!anyDuplicated(colnames(design)))
    subtract.design <- design[,coefsToSubtract]
    keep.design <- design[,setdiff(colnames(design), colnames(subtract.design))]
    limma::removeBatchEffect(x, design=keep.design, covariates=subtract.design, ...)
}

#' Implementation of `limma::voom()` that uses an offset matrix
#'
#' Unlike [limma::voom()], this function can use an offset matrix
#' instead of normalization factors. This is useful for using limma
#' with gene-level estimated counts from Salmon or Kallisto, for which
#' an offset matrix can be produced from the effective gene lengths.
#'
#' @param counts,design,normalize.method,span,plot,save.plot,...
#'     These arguments have the same meaning as in `limma::voom()`.
#' @param offset An offset matrix to be used in computing log2 CPM
#'     values. This is optional only if `counts` is a
#'     [edgeR::DGEList()].
#'
#' @return An EList object like that returned by [limma::voom()], but
#'     with log2 CPM values computed using the offset matrix.
#'
#' If `counts` is a DGEList object that does not contain an offset
#' matrix, one will be generated from the normalized library sizes,
#' yielding identical behavior to `limma::voom()`. Hence, this is
#' generally usable as a drop-in replacement for `limma::voom()`.
#'
#' @examples
#'
#' # TODO: Copy from voom
#'
#' @seealso [edgeR::getOffset()]
#'
#' @export
voomWithOffset <-
    function (counts, design = NULL, offset, normalize.method = "none",
              span = 0.5, plot = FALSE, save.plot = FALSE, ...)
{
    ## TODO: Port over new save.plot option


    ## TODO: Rewrite not to require a DGEList
    req_ns("limma", "edgeR", "BiocGenerics")
    out <- list()
    if (is(counts, "DGEList")) {
        out$genes <- counts$genes
        out$targets <- counts$samples
        if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >
            0)
            design <- model.matrix(~group, data = counts$samples)
        if (missing(offset)) {
            offset <- edgeR::expandAsMatrix(edgeR::getOffset(counts), dim(counts))
        }
        counts <- counts$counts
    }
    else {
        isExpressionSet <- suppressPackageStartupMessages(is(counts,
            "ExpressionSet"))
        if (isExpressionSet) {
            if (length(Biobase::fData(counts)))
                out$genes <- Biobase::fData(counts)
            if (length(Biobase::pData(counts)))
                out$targets <- Biobase::pData(counts)
            counts <- Biobase::exprs(counts)
        }
        else {
            counts <- as.matrix(counts)
        }
    }
    n <- nrow(counts)
    if (n < 2L)
        stop("Need at least two genes to fit a mean-variance trend")
    if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
    }

    effective.lib.size <- exp(offset)

    y <- log2((counts + 0.5)/(effective.lib.size + 1) * 1e+06)
    y <- limma::normalizeBetweenArrays(y, method = normalize.method)
    fit <- limma::lmFit(y, design, ...)
    if (is.null(fit$Amean))
        fit$Amean <- BiocGenerics::rowMeans(y, na.rm = TRUE)
    sx <- fit$Amean + mean(log2(effective.lib.size + 1)) - log2(1e+06)
    sy <- sqrt(fit$sigma)
    allzero <- rowSums(counts) == 0
    if (any(allzero)) {
        sx <- sx[!allzero]
        sy <- sy[!allzero]
    }
    l <- lowess(sx, sy, f = span)
    if (plot) {
        plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
            pch = 16, cex = 0.25)
        title("voom: Mean-variance trend")
        lines(l, col = "red")
    }
    f <- approxfun(l, rule = 2)
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coef[, j, drop = FALSE] %*%
            t(fit$design[, j, drop = FALSE])
    } else {
        fitted.values <- fit$coef %*% t(fit$design)
    }
    fitted.cpm <- 2^fitted.values
    # fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.count <- 1e-06 * fitted.cpm * (effective.lib.size + 1)
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    out$E <- y
    out$weights <- w
    out$design <- design
    out$effective.lib.size <- effective.lib.size
    if (is.null(out$targets)) {
        out$targets <- data.frame(lib.size = exp(BiocGenerics::colMeans(offset)))
    } else {
        out$targets$lib.size <- exp(BiocGenerics::colMeans(offset))
    }
    if (save.plot) {
        out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )",
            ylab = "Sqrt( standard deviation )")
        out$voom.line <- l
    }
    new("EList", out)
}

# TODO: Combine help text for all voom-related functions into one page

#' Implementation of `limma::voomWithQualityWeights()` that uses an offset matrix
#'
#' This is a combination of [limma::voomWithQualityWeights()] and
#' [voomWithOffset()].
#'
#' @param
#'     counts,design,normalize.method,span,plot,var.design,method,maxiter,tol,trace,replace.weights,col,...
#'     These arguments have the same meaning as in
#'     `limma::voomWithQualityWeights()`.
#' @param offset An offset matrix to be used in computing log2 CPM
#'     values. This is optional only if `counts` is a
#'     [edgeR::DGEList()].
#'
#' @return An EList object with weights, like that returned by
#'     [limma::voomWithQualityWeights()], but with log2 CPM values
#'     computed using the offset matrix.
#'
#' If `counts` is a DGEList object that does not contain an offset
#' matrix, one will be generated from the normalized library sizes,
#' yielding identical behavior to `limma::voom()`. Hence, this is
#' generally usable as a drop-in replacement for `limma::voom()`.
#'
#' @examples
#'
#' # TODO: Copy from voomWithQualityWeights
#'
#' @seealso [edgeR::getOffset()]
#'
#' @export
voomWithQualityWeightsAndOffset <-
    function (counts, design = NULL, offset, normalize.method = "none",
              plot = FALSE, span = 0.5, var.design = NULL, method = "genebygene",
              maxiter = 50, tol = 1e-10, trace = FALSE, replace.weights = TRUE,
              col = NULL, ...)
{
    req_ns("limma", "edgeR", "BiocGenerics")
    if (missing(offset) && is(counts, "DGEList")) {
        offset <- edgeR::expandAsMatrix(edgeR::getOffset(counts), dim(counts))
    }
    if (plot) {
        oldpar <- par(mfrow = c(1, 2))
        on.exit(par(oldpar))
    }
    v <- voomWithOffset(counts, design = design, offset = offset, normalize.method = normalize.method,
        plot = FALSE, span = span, ...)
    aw <- limma::arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, var.design = var.design)
    v <- voomWithOffset(counts, design = design, weights = aw, offset = offset,
        normalize.method = normalize.method, plot = plot, span = span,
        ...)
    aw <- limma::arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, trace = trace, var.design = var.design)
    wts <- limma::asMatrixWeights(aw, dim(v)) * v$weights
    attr(wts, "arrayweights") <- NULL
    if (plot) {
        barplot(aw, names = 1:length(aw), main = "Sample-specific weights",
            ylab = "Weight", xlab = "Sample", col = col)
        abline(h = 1, col = 2, lty = 2)
    }
    if (replace.weights) {
        v$weights <- wts
        v$sample.weights <- aw
        return(v)
    } else {
        return(wts)
    }
}

#' Alternate `limma::duplicateCorrelation()` and `limma::voom()` until convergence
#'
#' This function runs [limma::voom()] and
#' [limma::duplicateCorrelation()] repeatedly until two consecutive
#' runs yield the same correlation value (within `tol` of each other),
#' indicating that the iteration has converged.
#'
#' @param counts,design,plot,... These are passed to [limma::voom()].
#' @param block,trim These are passed to
#'     [limma::duplicateCorrelation()].
#' @param voom.fun Function to use in place of `limma::voom()`. For
#'     example, you could use `limma::voomWithQualityWeights`.
#'     Obvoiusly, this function should conform to the same interface
#'     as `limma::voom()`, accepting the same arguments and returning
#'     an EList.
#' @param dupCor.fun Function to use in place of
#'     `limma::duplicateCorrelation()`. I don't currently know of any
#'     suitable values for this other than
#'     `limma::duplicateCorrelation()`.
#' @param initial.correlation Initial correlation value to use for the
#'     first run of `limma::voom()`. There is no need to set this
#'     value, but setting it to a reasonable guess might allow the
#'     iteration to converge sooner.
#' @param maxiter Maximum number of iterations. If convergence has not
#'     been achieved after this many iterations, the final result will
#'     be returned with a warning.
#' @param tol Convergence tolerance. If a run of
#'     `limma::duplicateCorrelation()` yeilds a correlation value
#'     within `tol` of the previous correlation value, the result is
#'     considered converged, and iteration is stopped.
#' @param verbose If TRUE, give messages indicating the progress of
#'     the iteration.
#' @return The EList object returned by the final iteration of
#'     `limma::voom()`.
#'
#' It has previously been recommended in the past to run voom and
#' duplicateCorrelation twice each, with the implication that further
#' iterations would see diminishing returns (see
#' https://support.bioconductor.org/p/59700/#67620). However, rather
#' than assume two iterations is always sufficient, this function
#' instead keeps iterating until convergence is actually observed.
#'
#' Note that I am not aware of any proof that iterating to convergence
#' is guaranteed to produce a unique solution. The optimization may
#' well be non-convex, in which case this naive algorithm could
#' converge to a local optimum rather than a global one. Empirically,
#' differences between successive runs get very small after the second
#' run, so it seems stable in practice.
#'
#' @examples
#'
#' # TODO: Steal from voom/dupcor examples
#'
#' @export
voomWithDuplicateCorrelation <- function(counts, design = NULL, plot = FALSE, block = NULL, trim = 0.15,
                                         voom.fun=voomWithOffset,
                                         dupCor.fun=limma::duplicateCorrelation,
                                         initial.correlation=0,
                                         maxiter=5, tol=1e-6, verbose=TRUE, ...) {
    req_ns("limma", "glue")
    assert_that(maxiter >= 1)
    assert_that(is.finite(maxiter))

    if (maxiter < 2) {
        warning("Using less than 2 iterations of voom and duplicateCorrelation is not recommended.")
    }
    # Do first iteration
    prev.cor <- 0
    iter.num <- 0
    if (initial.correlation != 0 && !is.null(block)) {
        elist <- voom.fun(counts, design=design, plot=FALSE, block=block, correlation=initial.correlation, ...)
    } else {
        elist <- voom.fun(counts, design=design, plot=FALSE, ...)
    }
    if (is.null(block)) {
        warning("Not running duplicateCorrelation because block is NULL.")
    }
    if (verbose) {
        message(glue::glue("Initial guess for duplicate correlation before 1st iteration: {initial.correlation}"))
    }
    dupcor <- dupCor.fun(elist, design, block=block, trim=trim)
    iter.num <- iter.num + 1
    if (verbose) {
        message(glue::glue("Duplicate correlation after {toOrdinal(iter.num)} iteration: {dupcor$consensus.correlation}"))
    }
    while (iter.num < maxiter) {
        if (!is.null(tol) && is.finite(tol)) {
            if (abs(dupcor$consensus.correlation - prev.cor) <= tol) {
                if (verbose) {
                    message(glue::glue("Stopping after {toOrdinal(iter.num)} iteration because tolerance threshold was reached."))
                }
                break
            }
        }
        prev.cor <- dupcor$consensus.correlation
        elist <- voom.fun(counts, design=design, plot=FALSE, block=block, correlation=prev.cor, ...)
        dupcor <- dupCor.fun(elist, design, block=block, trim=trim)
        iter.num <- iter.num + 1
        if (verbose) {
            message(glue::glue("Duplicate correlation after toOrdinal(iter.num) iteration: dupcor$consensus.correlation"))
        }
    }
    elist <- voom.fun(counts, design=design, plot=plot, block=block, correlation=dupcor$consensus.correlation, ...)
    for (i in names(dupcor)) {
        elist[[i]] <- dupcor[[i]]
    }
    return(elist)
}

#' Variant of `limma::eBayes()` that sets the `proportion` argument automatically.
#'
#' Limma's [limma::eBayes()] function can calculate log-odds scores
#' (a.k.a. B statistics), which depend on an estimated proportion of
#' differentially expressed genes that the user specifies *a priori*.
#' Limma also provides methods for calculating this proportion via
#' [limma::propTrueNull()], and other packages also provide similar
#' functions. This function combines the two in order to set the
#' proportion automatically.
#'
#' The B statistics produced when running `limma::topTable()` on the
#' value returned by this function should be directly interpretable as
#' estimated log odds, as long as the distribituion of p-values is
#' well-behaved.
#'
#' Note that if the p-value distribution is highly atypical, a
#' proportion that is not between 0 and 1 could be estimated. If this
#' happens, and error will be thrown.
#'
#' @param prop.method The method by which to calculate the proportion
#'     of true null hypotheses (i.e. non-differentially expressed
#'     genes). This can be a string to be used as the `method`
#'     argument to [limma::propTrueNull()], or a function that accepts
#'     a vector of p-values and returns a single value between 0 and
#'     1. This argument can only be passed by name.
#' @param ... All other arguments are passed to [limma::eBayes()].
#'
#' @return See [limma::eBayes()].
#'
#' @examples
#'
#' #TODO Steal from eBayes
#'
#' @export
eBayes_auto_proportion <- function(..., prop.method="lfdr") {
    req_ns("limma")
    eb <- limma::eBayes(...)
    if (is.function(prop.method)) {
        ptn <- prop.method(eb$p.value)
    } else {
        ptn <- limma::propTrueNull(eb$p.value, method=prop.method)
    }
    assert_that(ptn > 0,ptn < 1)
    limma::eBayes(..., proportion=1-ptn)
}

#' Get a table of MDS values, with proper column names.
#'
#' This runs [limma::plotMDS()], but suppresses the generation of the
#' plot and instead returns the MDS dimensions in a matrix.
#'
#' @param x The object to run [limma::plotMDS()] on.
#' @param k The number of MDS dimensions to return. If not specified,
#'     the maximum possible number will be returned.
#' @param ... Additional arguments to [limma::plotMDS()].
#'
#' @return A matrix with `k` columns, and `ncol(x)` rows containing
#'     the MDS dimensions, with each column named "DimN", where N is
#'     the number of that dimension.
#'
#' @examples
#'
#' # TODO Steal from plotMDS
#'
#' @export
get_mds <- function(x, k, ...) {
    req_ns("limma")
    dmat <- limma::plotMDS(x, ..., plot = FALSE)$distance.matrix %>% as.dist
    max_k <- attr(dmat, "Size") - 1
    if (missing(k)) {
        k <- attr(dmat, "Size") - 1
    } else if (k > max_k) {
        warning(glue("Number of requested dimensions ({k}) is greater than the number available ({max_k}). Returning all dimensions."))
        k <- max_k
    }
    mds <- cmdscale(dmat, k=k, eig=TRUE)
    mds$points %>% add_numbered_colnames("Dim")
}
