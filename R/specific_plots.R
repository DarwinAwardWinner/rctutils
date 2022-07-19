#' Create an annotatied p-value histogram plot
#'
#' This function takes a vector of p-values and plots a histogram of
#' them using ggplot. It will add horizontal lines at 1 and the
#' estimated proportion of true null hypotheses for reference.
#'
#' @param pvals A vector of p-values
#' @param nbins Number of bins in the histogram
#' @param ptn Estimated proportion of true null hypotheses. If not
#'     provided, it will be estimated using [limma::propTrueNull()].
#'     This can also be a function that takes one argument (the vector
#'     of p-values) and returns the estimated proportion of true null
#'     hypotheses, which should be a single number between 0 and 1.
#' @return a ggplot object
#'
#' @importFrom rlang is_scalar_double
#' @importFrom glue glue
#' @export
plot_pval_hist <- function(pvals, nbins = 100, ptn = limma::propTrueNull) {
    req_ns("limma")
    if (is.function(ptn)) {
        ptn <- ptn(pvals)
    }
    assert_that(is_scalar_double(ptn))
    assert_that(ptn >= 0, ptn <= 1)
    df <- data.frame(p = pvals)
    linedf <- data.frame(y = c(1, ptn), Line = c("Uniform", "Est. Null") %>% factor(levels = unique(.)))
    ggplot(df) + aes_(x = ~p) +
        geom_histogram(aes_(y = ~..density..), binwidth = 1/nbins, boundary = 0) +
        geom_hline(aes_(yintercept = ~y, color = ~Line),
                   data = linedf, alpha = 0.5, show.legend = TRUE) +
        scale_color_manual(name = "Ref. Line", values = c("blue", "red")) +
        xlim(0,1) +
        ggtitle(
            "P-value distribution",
            subtitle = glue("(Est. {format(100 * (1-ptn), digits = 3)}% non-null.)")
        ) +
        expand_limits(y = c(0, 1.25)) +
        xlab("p-value") + ylab("Relative frequency") +
        theme(legend.position = c(0.95, 0.95),
              legend.justification = c(1,1))
}

#' Variant of `GGally::ggduo()` with separate arguments for  `dataX` and `dataY`
#'
#' @export
ggduo_dataXY <- function(dataX, dataY, extraData = NULL, ...) {
    req_ns("GGally")
    assert_that(ncol(dataX) > 0)
    assert_that(ncol(dataY) > 0)
    ## TODO: Ensure no repeated colnames
    alldata <- cbind(dataX, dataY)
    if (!is.null(extraData)) {
        alldata <- cbind(alldata, extraData)
    }
    GGally::ggduo(alldata, columnsX = colnames(dataX), columnsY = colnames(dataY), ...)
}

#' ggplot version of `edgeR::plotBCV()`.
#'
#' @param y A DGEList, or a data frame of the kind returned by
#'     [getBCVTable()], which should contain columns named "logCPM",
#'     "CommonBCV", "TrendBCV", "eBayesBCV", and optionally "RawBCV".
#' @param xlab,ylab Axis labels
#' @param rawdisp,... Additional arguments are passed to
#'     [getBCVTable()]. Note that this function is only called if `y`
#'     is a DGEList. Otherwise, these additional arguments are
#'     ignored.
#' @return A ggplot object.
#'
#' @export
ggplotBCV <- function(y, xlab = "Average log CPM", ylab = "Biological coefficient of variation", rawdisp = NULL, ...) {
    req_ns("reshape2")
    if (is(y, "DGEList")) {
        disptable <- getBCVTable(y, ..., rawdisp = rawdisp)
    } else {
        disptable <- as.data.frame(y)
        assert_that(all(c("logCPM", "CommonBCV", "TrendBCV", "eBayesBCV") %in% names(disptable)))
    }

    # Reduce the number of points to plot for each line for performance
    # reasons
    npoints <- c(Common = 2, Trend = 500)
    disp.line.table <-
        disptable %>%
        select_(~logCPM, ~TrendBCV, ~CommonBCV) %>%
        ## TODO: Replace with tidyr::gather
        reshape2::melt(id.vars = "logCPM", variable.name = "DispType", value.name = "BCV") %>%
        mutate_(DispType = ~str_replace(DispType, "BCV$", "")) %>%
        group_by_(~DispType) %>%
        do_(~{
            approx(x = .$logCPM, y = .$BCV, n = npoints[.$DispType[1]]) %>%
                data.frame(logCPM = .$x, BCV = .$y)
        })

    p <- ggplot(disptable) +
        aes_(x = ~logCPM)
    if ("RawBCV" %in% names(disptable)) {
        p <- p +
            geom_point(aes_(y = ~RawBCV), size = 0.4, color = "black") +
            geom_density2d(aes_(y = ~RawBCV), color = "gray30", n = 512) +
            labs(subtitle = "Raw BCV (black) and eBayes-squeezed (blue)")
    }
    p <- p +
        geom_point(aes_(y = ~eBayesBCV), size = 0.1, color = "darkblue") +
        geom_density2d(aes_(y = ~eBayesBCV), color = "blue", n = 512) +
        geom_line(data = disp.line.table, aes_(x = ~logCPM, y = ~BCV, group = ~DispType), color = "white", size = 1.5, alpha = 0.5) +
        geom_line(data = disp.line.table, aes_(x = ~logCPM, y = ~BCV, linetype = ~DispType), color = "darkred", size = 0.5) +
        scale_linetype_manual(name = "Dispersion Type", values = c(Trend = "solid", Common = "dashed")) +
        labs(title = "BCV plot", x = xlab, y = ylab)
    p
}

# Get log2-fold-change y-axis position of normalization lines between two
# samples for a list of DGEList objects.
#' @export
getNormLineData <- function(dgelists, s1, s2) {
    req_ns("rlang")
    assert_that(length(dgelists) >= 1)
    assert_that(rlang::is_dictionaryish(dgelists))
    dgelists %>%
        sapply(. %>% {.$samples$norm.factors} %>% log2 %>% {.[s2] - .[s1]}) %>%
        data.frame(NormFactor=., NormType=names(dgelists))
}

# Get a curve representing the loess-based normaliation implied by the offsets
# of two samples in a DGEList object. N is the number of points to interpolate the curve at.
#' @export
getOffsetNormCurveData <- function(dge, s1, s2, n=1000) {
    req_ns("edgeR")
    assert_that(is.numeric(dge$offset))
    a <- edgeR::aveLogCPM(dge, dispersion=0.05, prior.count=0.5)
    # Need to subtract the library size difference out of the offset
    raw.offset <- dge$offset %>% {.[,s2] - .[,s1]} %>% divide_by(log(2))
    lib.size.offset <- dge$samples$lib.size %>% {.[s2] / .[s1]} %>% log2
    x <- data.frame(A=a, Offset=raw.offset - lib.size.offset)
    f <- approxfun(x$A, x$Offset)
    data.frame(A=seq(from=min(x$A), to=max(x$A), length.out = n)) %>%
        mutate(M=f(A))
}
