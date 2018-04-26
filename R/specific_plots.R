# TODO: Make ggplot, GGally optional

#' Create an annotatied p-value histogram plot
#'
#' @export
plotpvals <- function(pvals, ptn=propTrueNull(pvals)) {
    df <- data.frame(p=pvals)
    linedf <- data.frame(y=c(1, ptn), Line=c("Uniform", "Est. Null") %>% factor(levels=unique(.)))
    ggplot(df) + aes(x=p) +
        geom_histogram(aes(y = ..density..), binwidth=0.01, boundary=0) +
        geom_hline(aes(yintercept=y, color=Line),
            data=linedf, alpha=0.5, show.legend=TRUE) +
        scale_color_manual(name="Ref. Line", values=c("blue", "red")) +
        xlim(0,1) + ggtitle(glue("P-value distribution (Est. {format(100 * (1-ptn), digits=3)}% signif.)")) +
        expand_limits(y=c(0, 1.25)) +
        xlab("p-value") + ylab("Relative frequency") +
        theme(legend.position=c(0.95, 0.95),
            legend.justification=c(1,1))
}

#' Variant of `GGally::ggduo()` with separate arguments for  `dataX` and `dataY`
#'
#' @export
ggduo.dataXY <- function(dataX, dataY, extraData=NULL, ...) {
    req_ns("GGally")
    assert_that(ncol(dataX) > 0)
    assert_that(ncol(dataY) > 0)
    alldata <- cbind(dataX, dataY)
    if (!is.null(extraData)) {
        alldata <- cbind(alldata, extraData)
    }
    GGally::ggduo(alldata, columnsX=colnames(dataX), columnsY=colnames(dataY), ...)
}

#' ggplot version of `edgeR::plotBCV()`.
#'
#' Additional arguments passed to [getBCVTable()].
#'
#' @include edgeR_utils.R
#' @export
ggplotBCV <- function(y, xlab="Average log CPM", ylab="Biological coefficient of variation", rawdisp=NULL, ...) {
    req_ns("ggplot2")
    if (is(y, "DGEList")) {
        disptable <- getBCVTable(y, ..., rawdisp=rawdisp)
    } else {
        disptable <- as.data.frame(y)
        assert_that(all(c("logCPM", "CommonBCV", "TrendBCV", "eBayesBCV") %in% names(disptable)))
    }

    # Reduce the number of points to plot for each line for performance
    # reasons
    npoints <- c(Common=2, Trend=500)
    disp.line.table <-
        disptable %>%
        select(logCPM, TrendBCV, CommonBCV) %>%
        melt(id.vars="logCPM", variable.name="DispType", value.name = "BCV") %>%
        mutate(DispType=str_replace(DispType, "BCV$", "")) %>%
        group_by(DispType) %>%
        do({
            approx(x=.$logCPM, y=.$BCV, n=npoints[.$DispType[1]]) %>%
                data.frame(logCPM=.$x, BCV=.$y)
        })

    p <- ggplot(disptable) +
        aes(x=logCPM)
    if ("RawBCV" %in% names(disptable)) {
        p <- p +
            geom_point(aes(y=RawBCV), size=0.4, color="black") +
            geom_density2d(aes(y=RawBCV), color="gray30", n=512) +
            labs(subtitle="Raw BCV (black) and eBayes-squeezed (blue)")
    }
    p <- p +
        geom_point(aes(y=eBayesBCV), size=0.1, color="darkblue") +
        geom_density2d(aes(y=eBayesBCV), color="blue", n=512) +
        geom_line(data=disp.line.table, aes(x=logCPM, y=BCV, group=DispType), color="white", size=1.5, alpha=0.5) +
        geom_line(data=disp.line.table, aes(x=logCPM, y=BCV, linetype=DispType), color="darkred", size=0.5) +
        scale_linetype_manual(name="Dispersion Type", values=c(Trend="solid", Common="dashed")) +
        labs(title="BCV plot", x=xlab, y=ylab)
    p
}
