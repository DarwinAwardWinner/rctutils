## limma:

#' @export
subtractCoefs <- function(x, design, coefsToSubtract, ...) {
    req_ns("limma")
    assert_that(!anyDuplicated(colnames(design)))
    subtract.design <- design[,coefsToSubtract]
    keep.design <- design[,setdiff(colnames(design), colnames(subtract.design))]
    limma::removeBatchEffect(x, design=keep.design, covariates=subtract.design, ...)
}

#' @importFrom graphics abline barplot lines par title
#' @export
voomWithOffset <- function (dge, design = NULL, offset=expandAsMatrix(getOffset(dge), dim(dge)),
                            normalize.method = "none", plot = FALSE, span = 0.5, ...) {
    req_ns("limma")
    req_ns("BiocGenerics")
    out <- list()
    out$genes <- dge$genes
    out$targets <- dge$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 0)
        design <- model.matrix(~group, data = counts$samples)
    counts <- dge$counts
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
    new("EList", out)
}

#' Version of voom that uses an offset matrix instead of lib sizes
#'
#' @export
voomWithQualityWeightsAndOffset <-
    function (dge, design = NULL,
              offset=expandAsMatrix(getOffset(dge), dim(dge)),
              normalize.method = "none",
              plot = FALSE, span = 0.5, var.design = NULL, method = "genebygene",
              maxiter = 50, tol = 1e-10, trace = FALSE, replace.weights = TRUE,
              col = NULL, ...)
{
    counts <- dge$counts
    if (plot) {
        oldpar <- par(mfrow = c(1, 2))
        on.exit(par(oldpar))
    }
    v <- voomWithOffset(dge, design = design, offset = offset, normalize.method = normalize.method,
        plot = FALSE, span = span, ...)
    aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, var.design = var.design)
    v <- voomWithOffset(dge, design = design, weights = aw, offset = offset,
        normalize.method = normalize.method, plot = plot, span = span,
        ...)
    aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, trace = trace, var.design = var.design)
    wts <- asMatrixWeights(aw, dim(v)) * v$weights
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

#' Convenience function for alternating dupCor and voom until convergence
#'
#' @export
voomWithDuplicateCorrelation <-
    function(counts, design = NULL, plot = FALSE, block = NULL, trim = 0.15,
             voom.fun=voom, dupCor.fun=duplicateCorrelation, initial.correlation=0,
             niter=5, tol=1e-6, verbose=TRUE, ...)
{
    assert_that(niter >= 1)
    assert_that(is.finite(niter))

    if (niter < 2) {
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
        message(glue("Initial guess for duplicate correlation before 1st iteration: {initial.correlation}"))
    }
    dupcor <- dupCor.fun(elist, design, block=block, trim=trim)
    iter.num <- iter.num + 1
    if (verbose) {
        message(glue("Duplicate correlation after {toOrdinal(iter.num)} iteration: {dupcor$consensus.correlation}"))
    }
    while (iter.num < niter) {
        if (!is.null(tol) && is.finite(tol)) {
            if (abs(dupcor$consensus.correlation - prev.cor) <= tol) {
                if (verbose) {
                    message(glue("Stopping after {toOrdinal(iter.num)} iteration because tolerance threshold was reached."))
                }
                break
            }
        }
        prev.cor <- dupcor$consensus.correlation
        elist <- voom.fun(counts, design=design, plot=FALSE, block=block, correlation=prev.cor, ...)
        dupcor <- dupCor.fun(elist, design, block=block, trim=trim)
        iter.num <- iter.num + 1
        if (verbose) {
            message(glue("Duplicate correlation after toOrdinal(iter.num) iteration: dupcor$consensus.correlation"))
        }
    }
    elist <- voom.fun(counts, design=design, plot=plot, block=block, correlation=dupcor$consensus.correlation, ...)
    for (i in names(dupcor)) {
        elist[[i]] <- dupcor[[i]]
    }
    return(elist)
}
