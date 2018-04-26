#' @export
aveLogCPMWithOffset <- function(y, ...) {
    UseMethod("aveLogCPMWithOffset")
}

#' @export
aveLogCPMWithOffset.default <- function(y, offset = NULL, prior.count = 2,
                                        dispersion = NULL, weights = NULL, ...) {
    aveLogCPM(y, lib.size = NULL, offset = offset, prior.count = prior.count,
              dispersion = dispersion, weights = weights, ...)
}

#' @importFrom edgeR getOffset expandAsMatrix
#' @export
aveLogCPMWithOffset.DGEList <- function(y, offset = expandAsMatrix(getOffset(y), dim(y)),
                                        prior.count = 2, dispersion = NULL, ...) {
    if (is.null(dispersion)) {
        dispersion <- y$common.dispersion
    }
    offsetMat <- offset
    aveLogCPMWithOffset(
        y$counts, offset = offsetMat, prior.count = prior.count,
        dispersion = dispersion, weights = y$weights)
}


#' @export
cpmWithOffset <- function(dge, offset=getOffset(dge),
                         log = FALSE, prior.count = 0.25, preserve.mean=TRUE, ...) {
    if (preserve.mean) {
        dge <- scaleOffset(dge, offset)
    } else {
        dge$offset <- offset
    }
    cpm(dge$counts, lib.size=exp(getOffset(dge)), log=log, prior.count=prior.count, ...)
}

#' @export
estimateDispByGroup <- function(dge, group=as.factor(dge$samples$group), batch, ...) {
    assert_that(nlevels(group) > 1)
    assert_that(length(group) == ncol(dge))
    if (!is.list(batch)) {
        batch <- list(batch=batch)
    }
    batch <- as.data.frame(batch)
    assert_that(nrow(batch) == ncol(dge))
    colnames(batch) %>% make.names(unique=TRUE)
    igroup <- seq_len(ncol(dge)) %>% split(group)
    bplapply(igroup, function(i) {
        group.dge <- dge[,i]
        group.batch <- droplevels(batch[i,, drop=FALSE])
        group.batch <- group.batch[sapply(group.batch, . %>% unique %>% length %>% is_greater_than(1))]
        group.vars <- names(group.batch)
        if (length(group.vars) == 0)
            group.vars <- "1"
        group.batch.formula <- as.formula(str_c("~", str_c(group.vars, collapse="+")))
        des <- model.matrix(group.batch.formula, group.batch)
        estimateDisp(group.dge, des, ...)
    })
}

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
        y %<-% edgeR::estimateDisp(y, design=design, ...)
    }
    assert_that(!is.null(y$prior.df))
    if (all(y$prior.df == 0)) {
        if (missing(rawdisp) || is.null(rawdisp)) {
            rawdisp <- y
        }
        y %<-% edgeR::estimateDisp(y, design=design, ...)
    }
    # Get raw (unsqueezed dispersions)
    if (missing(rawdisp) || is.na(rawdisp)) {
        # Estimate raw disperions using given design
        y.raw %<-% estimateDisp(y, design=design, prior.df=0)
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
    disptable <- y %>% as.list %$% data.frame(
        logCPM=AveLogCPM,
        CommonBCV=common.dispersion %>% sqrt,
        TrendBCV=trended.dispersion %>% sqrt,
        PriorDF=prior.df,
        eBayesBCV=tagwise.dispersion %>% sqrt)
    if (!is.null(y.raw)) {
        disptable$RawBCV <- y.raw$tagwise.dispersion %>% sqrt
    }
    return(disptable)
}

#' Returns mds values, with proper colnames.
#'
#' @export
get_mds <- function(x, k) {
    dmat <- suppressPlot(plotMDS(x, top=5000)$distance.matrix) %>% as.dist
    if (missing(k)) {
        k <- attr(dmat, "Size") - 1
    }
    mds <- cmdscale(dmat, k=k, eig=TRUE)
    mds$points %>% add.numbered.colnames("Dim") %>% cbind(sample.table, .)
}
