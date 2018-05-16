## TODO: Make ggplot an optional dep

# Utilities for ggplot2 corrdinate transformation

#' @importFrom glue glue
#' @export
power_trans <- function(pow) {
    req_ns("scales")
    name <- glue("^{pow}")
    scales::trans_new(
        name,
        transform = function(x) x ^ pow,
        inverse = function(x) x ^ (1/pow),
        domain = c(0,Inf))
}

#' @importFrom glue glue
#' @export
clamp_trans <- function(lower_threshold = 0, upper_threshold = 1) {
    req_ns("scales")
    name <- glue("Clamp values outside of [{lower_threshold}, {upper_threshold}]")
    scales::trans_new(
        name,
        transform = function(x) pmin(upper_threshold, pmax(lower_threshold, x)),
        # transform is only invertible for part of the range
        inverse = identity)
}

#' @export
neglog_trans <- function(base = exp(1)) {
    req_ns("scales")
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(
        paste0("negativelog-", format(base)), trans, inv,
        scales::log_breaks(base = base),
        domain = c(1e-100, Inf))
}

#' @export
discrete_gradient <- function(n) {
    req_ns("scales")
    scales::seq_gradient_pal(low = "#132B43", high = "#56B1F7")(seq(0,1, length.out = n))
}

# Always returns a list of ggplot objects. Flattens nested lists,
# encapsulates single plots into a 1-element list, ensures that all
# elements are ggplots.

#' @export
get.ggplots <- function(plots) {
    UseMethod("get.ggplots")
}

#' @export
get.ggplots.default <- function(plots) {
    stop(glue("Don't know how to get ggplots from an object of class {deparse(class(plots)[1])}"))
}

#' @export
get.ggplots.gg <- function(plots) {
    list(plots)
}

#' @export
get.ggplots.list <- function(plots) {
    plotlists <- lapply(plots, get.ggplots)
    do.call(c, plotlists)
}

#' @export
ggprint <- function(plots, device = dev.cur(), closedev, printfun = print) {
    req_ns("withr")
    p <- get.ggplots(plots)
    with_dev(device, lapply(p, printfun), closedev)
    invisible(p)
}

# Printer function for ggplotly, to be passed as the prinfun for
# ggprint.

#' @export
ggplotly.printer <- function(...) {
    req_ns("plotly")
    dots <- list(...)
    function(p) {
        args <- c(list(p = p), dots)
        print(do.call(plotly::ggplotly, args))
    }
}
