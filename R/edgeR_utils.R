#' Version of `edgeR::aveLogCPM()` that uses offsets by default
#'
#' This is the same as [edgeR::aveLogCPM()], except that when called
#' on a [edgeR::DGEList()] object that contains an offset matrix, that
#' matrix will be used for normalization instead of the library sizes
#' and normalization factors. If called on a matrix of counts, the
#' offset must be supplied manually.
#'
#' @param y,offset,prior.count,dispersion,weights,... All arguments
#'     have the same meaning as in [edgeR::aveLogCPM()], except that
#'     if `y` is a `edgeR::DGEList()`, the offset matrix will be used
#'     in favor of the sample normalized library sizes.
#' @return See [edgeR::aveLogCPM()].
#'
#' Note that `edgeR::aveLogCPM()` also accepts an offset argument. The
#' only difference is that it doesn't use the offset matrix from a
#' DGEList by default, whereas this function does.
#'
#' @examples
#'
#' # TODO Steal from aveLogCPM and maybe scaleOffset
#'
#' @seealso [edgeR::getOffset()]
#'
#' @export
aveLogCPMWithOffset <- function(y, ...) {
    UseMethod("aveLogCPMWithOffset")
}

#' @export
aveLogCPMWithOffset.default <- function(y, offset = NULL, prior.count = 2,
                                        dispersion = NULL, weights = NULL, ...) {
    req_ns("edgeR")
    edgeR::aveLogCPM(y, lib.size = NULL, offset = offset, prior.count = prior.count,
                     dispersion = dispersion, weights = weights, ...)
}

#' @export
aveLogCPMWithOffset.DGEList <- function(y, offset = edgeR::expandAsMatrix(edgeR::getOffset(y), dim(y)),
                                        prior.count = 2, dispersion = NULL, ...) {
    req_ns("edgeR")
    if (is.null(dispersion)) {
        dispersion <- y$common.dispersion
    }
    offsetMat <- offset
    aveLogCPMWithOffset(
        y$counts, offset = offsetMat, prior.count = prior.count,
        dispersion = dispersion, weights = y$weights, ...)
}

#' Version of `edgeR::cpm()` that uses offsets
#'
#' This is the same as [edgeR::cpm()], except that when called on a
#' [edgeR::DGEList()] object that contains an offset matrix, that
#' matrix will be used by default for normalization instead of the
#' library sizes and normalization factors. If called on a matrix of
#' counts, the offset must be supplied manually.
#'
#' @param y,log,prior.count,... These arguments have the same meaning
#'     as in [edgeR::cpm()], except that if `y` is a
#'     `edgeR::DGEList()`, the offset matrix will be used in favor of
#'     the sample normalized library sizes. Note that `lib.size` is
#'     not an accepted argument.
#' @param offset The offset matrix to use. If `y` is a DGEList, this
#'     is determined automatically.
#' @param preserve.mean If TRUE (the default), scale each row of
#'     `offset` such that the average CPM is the same as it would be
#'     without the offset.
#' @param ... Additional arguments are passed to [edgeR::cpm()].
#' @return A numeric matrix of CPM values, normalized using the
#'     provided offset matrix.
#'
#' If `y` is a DGEList with no offset matrix, `offset` is generated
#' from the normalized library sizes.
#'
#' @examples
#'
#' # TODO Steal from cpm and maybe scaleOffset
#'
#' @seealso [edgeR::getOffset()], [edgeR::scaleOffset()], [edgeR::cpm()]
#'
#' @export
cpmWithOffset <- function(y, offset, ...) {
    UseMethod("cpmWithOffset")
}

#' @export
cpmWithOffset.default <- function(y, offset, log = FALSE, prior.count = 0.25, ..., preserve.mean = TRUE) {
    req_ns("edgeR")
    if (preserve.mean) {
        offset <- edgeR::scaleOffset(y, offset)
    }
    edgeR::cpm(y$counts, lib.size = exp(offset), log = log,
               prior.count = prior.count, ...)
}

#' @export
cpmWithOffset.DGEList <- function(y, offset = edgeR::expandAsMatrix(edgeR::getOffset(y), dim(y)),
                                  log = FALSE, prior.count = 0.25, ..., preserve.mean = TRUE) {
    req_ns("edgeR")
    cpmWithOffset(t$counts, offset, log = log,
                  prior.count = prior.count, preserve.mean = preserve.mean, ...)
}

#' Estimate edgeR dispersions separately for each group
#'
#' This takes a DGEList and a grouping factor and estimates the dispersion
#'
#' @param dge A DGEList object.
#' @param group A factor to use in spliting the samples in `dge` into
#'     groups.
#' @param batch A vector or list of multiple vectors to be treated as
#'     batch effects when estimating dispersions.
#' @param ... Additional arguments are passed to
#'     [edgeR::estimateDisp()].
#' @param BPPARAM A BiocParallelParam instance to be used for
#'     parallelization.
#' @return A named list with one DGEList for each group. Each DGEList
#'     contains the samples and estimated dispersions for that group.
#'
#' @seealso [edgeR::estimateDisp()]
#'
#' @export
estimateDispByGroup <- function(dge, group = as.factor(dge$samples$group), batch, ...,
                                BPPARAM = BiocParallel::bpparam()) {
    req_ns("edgeR", "BiocParallel", "magrittr")
    assert_that(nlevels(group) > 1)
    assert_that(length(group) == ncol(dge))
    if (!is.list(batch)) {
        batch <- list(batch = batch)
    }
    batch <- as.data.frame(batch)
    assert_that(nrow(batch) == ncol(dge))
    colnames(batch) %>% make.names(unique = TRUE)
    igroup <- seq_len(ncol(dge)) %>% split(group)
    BiocParallel::bplapply(igroup, function(i) {
        group.dge <- dge[,i]
        group.batch <- droplevels(batch[i,, drop = FALSE])
        group.batch <- group.batch[sapply(group.batch, . %>% unique %>% length %>%
                                                       is_greater_than(1))]
        group.vars <- names(group.batch)
        if (length(group.vars) == 0)
            group.vars <- "1"
        group.batch.formula <- as.formula(str_c("~", str_c(group.vars, collapse = "+")))
        des <- model.matrix(group.batch.formula, group.batch)
        edgeR::estimateDisp(group.dge, des, ...)
    }, BPPARAM=BPPARAM)
}

#' Get a table of genewise BCV values from a DGEList
#'
#' Given a DGEList object, this will extract the average logCPM as
#' well as the common, trended, and genewise dispersion values. In
#' addition, it will re-estimate the dispersions with no empirical
#' Bayes shrinkage (prior.df = 0) to obtain the raw genewise
#' dispersion values.
#'
#' @param y A DGEList. If it does not already contain dispersion
#'     values, they will be estimated.
#' @param design A design matrix to use when estimating dispersions.
#'     This is optional if `y` already contains a design.
#' @param ... Additional arguments are passed to
#'     [edgeR::estimateDisp()], if dispersion estimation is required
#'     (see below).
#' @param rawdisp Another DGEList, which should be identical to `y`
#'     except that it had dispersions estimated with `prior.count` =
#'     0. This argument is completely optional, and is only used for
#'     speeding up the function if the raw dispersions have already
#'     been estimated. To suppress estimation of raw dispersions
#'     entirely, set this to NA, NULL, or FALSE.
#'
#' (TODO: Take a long hard look at the logic in this function and make
#' sure it matches this documentation.)
#'
#' If `y` already contains estimated dispersions, then no other
#' arguments are required, and the existing dispersion values will be
#' used. If `y` does not contain dispersions, then they will be
#' estimated using `design` and `...`. Similarly, `design` is only
#' required if `y` does not already contain a design (check if
#' `y$design` is non-NULL).
#'
#' @export
getBCVTable <- function(y, design, ..., rawdisp) {
    req_ns("edgeR")
    assert_that(is(y, "DGEList"))
    design.passed <- !missing(design)
    if (!design.passed) {
        # Can be NULL
        design <- y$design
    }
    # Estimate dispersions now if they are not already present
    if (! all(c("common.dispersion", "trended.dispersion", "tagwise.dispersion") %in%
                  names(y))) {
        if (is.null(design) && !design.passed) {
            warning("Estimating dispersions with no design matrix")
        }
        y <- edgeR::estimateDisp(y, design = design, ...)
    }
    assert_that(!is.null(y$prior.df))
    # If y has prior.df == 0, then y is actually rawdisp
    if (all(y$prior.df == 0)) {
        if (missing(rawdisp) || is.null(rawdisp)) {
            rawdisp <- y
        }
        y <- edgeR::estimateDisp(y, design = design, ...)
    }
    # Get raw (unsqueezed dispersions)
    if (missing(rawdisp) || is.na(rawdisp)) {
        # Estimate raw disperions using given design
        y.raw <- edgeR::estimateDisp(y, design = design, prior.df = 0)
    } else if (is.null(rawdisp)) {
        # Explicitly passing NULL means no raw disperions are desired.
        y.raw <- NULL
    } else if (is(rawdisp, "DGEList")) {
        # Assume DGEList already contains raw dispersions
        assert_that(all(dim(rawdisp) == dim(y)))
        y.raw <- rawdisp
    } else {
        # Assume anything else is a numeric vector of raw dispersions
        assert_that(is.numeric(rawdisp),
            length(rawdisp) == nrow(y))
        y.raw <- y
        y.raw$tagwise.dispersion <- rawdisp
        y.raw$prior.df <- y.raw$prior.n <- rep(0, length(rawdisp))
    }
    assert_that(all(y.raw$prior.df == 0))
    y %<>% as.list
    disptable <- data.frame(
        logCPM = y$AveLogCPM,
        CommonBCV = y$common.dispersion %>% sqrt,
        TrendBCV = y$trended.dispersion %>% sqrt,
        PriorDF = y$prior.df,
        eBayesBCV = y$tagwise.dispersion %>% sqrt)
    if (!is.null(y.raw)) {
        disptable$RawBCV <- y.raw$tagwise.dispersion %>% sqrt
    }
    return(disptable)
}
