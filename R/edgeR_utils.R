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
        dispersion = dispersion, weights = y$weights)
}


#' @export
cpmWithOffset <- function(dge, offset = edgeR::expandAsMatrix(edgeR::getOffset(dge), dim(dge)),
                          log = FALSE, prior.count = 0.25, preserve.mean=TRUE, ...) {
    req_ns("edgeR")
    if (preserve.mean) {
        dge <- edgeR::scaleOffset(dge, offset)
    } else {
        dge$offset <- offset
    }
    edgeR::cpm(dge$counts, lib.size=exp(edgeR::getOffset(dge)), log=log, prior.count=prior.count, ...)
}

#' @importFrom magrittr is_greater_than
#' @export
estimateDispByGroup <- function(dge, group=as.factor(dge$samples$group), batch, ...) {
    req_ns("edgeR", "BiocParallel", "magrittr")
    assert_that(nlevels(group) > 1)
    assert_that(length(group) == ncol(dge))
    if (!is.list(batch)) {
        batch <- list(batch=batch)
    }
    batch <- as.data.frame(batch)
    assert_that(nrow(batch) == ncol(dge))
    colnames(batch) %>% make.names(unique=TRUE)
    igroup <- seq_len(ncol(dge)) %>% split(group)
    BiocParallel::bplapply(igroup, function(i) {
        group.dge <- dge[,i]
        group.batch <- droplevels(batch[i,, drop=FALSE])
        group.batch <- group.batch[sapply(group.batch, . %>% unique %>% length %>%
                                                       is_greater_than(1))]
        group.vars <- names(group.batch)
        if (length(group.vars) == 0)
            group.vars <- "1"
        group.batch.formula <- as.formula(str_c("~", str_c(group.vars, collapse="+")))
        des <- model.matrix(group.batch.formula, group.batch)
        edgeR::estimateDisp(group.dge, des, ...)
    })
}

#' @importFrom magrittr %>% %$% %<>%
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
        y <- edgeR::estimateDisp(y, design=design, ...)
    }
    assert_that(!is.null(y$prior.df))
    if (all(y$prior.df == 0)) {
        if (missing(rawdisp) || is.null(rawdisp)) {
            rawdisp <- y
        }
        y <- edgeR::estimateDisp(y, design=design, ...)
    }
    # Get raw (unsqueezed dispersions)
    if (missing(rawdisp) || is.na(rawdisp)) {
        # Estimate raw disperions using given design
        y.raw <- edgeR::estimateDisp(y, design=design, prior.df=0)
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
        logCPM=y$AveLogCPM,
        CommonBCV=y$common.dispersion %>% sqrt,
        TrendBCV=y$trended.dispersion %>% sqrt,
        PriorDF=y$prior.df,
        eBayesBCV=y$tagwise.dispersion %>% sqrt)
    if (!is.null(y.raw)) {
        disptable$RawBCV <- y.raw$tagwise.dispersion %>% sqrt
    }
    return(disptable)
}

#' Returns mds values, with proper colnames.
#'
#' @importFrom magrittr %>%
#' @export
get_mds <- function(y, k) {
    req_ns("limma")
    dmat <- suppressPlot(limma::plotMDS(y, top=5000)$distance.matrix) %>% as.dist
    if (missing(k)) {
        k <- attr(dmat, "Size") - 1
    }
    mds <- cmdscale(dmat, k=k, eig=TRUE)
    mds$points %>% add_numbered_colnames("Dim")
}
