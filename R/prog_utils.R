
#' @export
tsmsg <- function(...) {
    message(date(), ": ", ...)
}

# Use to assign to complex sub-expressions in the middle of a dplyr pipeline.
# For example, you can't easily do the following in the middle of a pipeline:
# "assays(x[[1]])$counts[3,5] <- 45". But now you can do it like: "x %>%
# assign_into(assays(.[[1]])$counts[3,5], 45) %>% another_fun() %>% ..."

#' @export
assign_into <- function(x, expr, value) {
    req_ns("lazyeval")
    expr <- lazyeval::lazy(expr)$expr
    f_eval(f_interp(~ x %>% { lazyeval::uq(expr) <- lazyeval::uq(value); . }))
}

#' @export
withGC <- function(expr) {
    on.exit(gc())
    return(expr)
}
