#' Like `message` but with a timestamp
#'
#' This is identical to the [message()] function, except that it
#' prepends a timestamp to every message.
#'
#' @param ... All arguments are passed directly to `message`.
#'
#' The time stamp is generated using [date()]. If you need anything
#' more complicated than this, it's probably easier to abondon this
#' function and just use `date()` and `message()` to manually add the
#' timestamp in the format you want.
#'
#' @examples
#'
#' tsmsg("Hello world")
#'
#' @seealso [message()], [date()]
#'
#' @export
tsmsg <- function(...) {
    message(date(), ": ", ...)
}

#' Assign into complex sub-expressions and return the whole object
#'
#' @noMd
#'
#' This function exists to facilitate assignment to parts of an object
#' in the middle of a magrittr pipeline. Normally this is disruptive,
#' since assignment returns the value that was assigned, rather than
#' the whole object.
#'
#' @param x The object to assign into. (Typically this argument is
#'     delivered via \code{\%\>\%()}.)
#' @param expr The left hand side of the assignment operation to be
#'     performed on \code{x}.
#' @param value The right hand side of the assignment operation to be
#'     performed.
#' @return The value of \code{x} after performing the assignment.
#'
#' As usual, \code{x}, the object being passed in, should be
#' referenced in both \code{expr} and \code{value} as \code{.}. In
#' fact, \code{expr} *must* contain at least one \code{.}, or else the
#' object will not be modified. (This is currently not checked.)
#'
#' Note that this function uses the lazyeval package rather than its
#' apparent successor, rlang, because rlang doesn't support
#' interpolating expressions on the left-hand-side of an assignment
#' operation: https://github.com/r-lib/rlang/issues/212.
#'
#' @examples
#'
#' library(magrittr)
#'
#' # Returns the entire list, not just the value of x
#' list(a = 1, b = 2, c = 3) %>% assign_into(.$x, 5)
#'
#' \dontrun{
#' # A more complex example of how this might be used in the middle of
#' # a pipeline. Imagine that x is a list of SummarizedExperiment
#' # objects, and for some reason one of the values in one of the
#' # assays in one of the objects is wrong and needs to be modified.
#' x %>% assign_into(assays(.[[1]])$counts[3,5], 45) %>% ...
#' }
#'
#' @seealso \code{\%\>\%} from the magrittr package.
#'
#' @export
assign_into <- function(x, expr, value) {
    req_ns("lazyeval")
    ## Required because `f_interp` doesn't seem to like `::`
    uq <- lazyeval::uq
    expr <- lazyeval::lazy(expr)$expr
    lazyeval::f_eval(lazyeval::f_interp(~( x %>% { uq(expr) <- uq(value); . })))
}

#' Evaluate an arbitrary mutating expression in a magrittr pipe
#'
#' @noMd
#'
#' This function exists to facilitate assignment to parts of an object
#' in the middle of a magrittr pipeline. Normally this is disruptive,
#' since assignment returns the value that was assigned, rather than
#' the whole object.
#'
#' @param x The object to assign into. (Typically this argument is
#'     delivered via \code{\%\>\%()}.)
#' @param expr The expression to be evaluated with \code{.} set to
#'     \code{x}. Typically this expression modifies \code{.} in some
#'     way, e.g. by performing an assignment with \code{.} on the
#'     left-hand side.
#' @return The value of \code{x} after evaluating \code{expr}.
#'
#' As usual, \code{x}, the object being passed in, should be
#' referenced in \code{expr} as \code{.}.
#'
#' Note that this function uses the lazyeval package rather than its
#' apparent successor, rlang, because rlang doesn't support
#' interpolating expressions on the left-hand-side of an assignment
#' operation: https://github.com/r-lib/rlang/issues/212.
#'
#' @examples
#'
#' library(magrittr)
#'
#' # Returns the entire list, not just a single element
#' list(a = 1, b = 2, c = 3) %>%
#'   dot_mutate(.$d <- 4) %>%
#'   dot_mutate({ .$e <- 5; .$f <- 6})
#'
#' # Compare:
#' list(a = 1, b = 2, c = 3) %>%
#'   {.$d <- 4} %>%
#'   print
#'
#' @export
dot_mutate <- function(x, expr) {
    req_ns("lazyeval")
    ## Required because `f_interp` doesn't seem to like `::`
    uq <- lazyeval::uq
    expr <- lazyeval::lazy(expr)$expr
    lazyeval::f_eval(lazyeval::f_interp(~( x %>% { uq(expr); . })))
}

#' Evaluate an expression and then collect garbage.
#'
#' Logically, this function is equivalent to [identity()], simply
#' returning the value of the expression passed to it. However, it
#' also ensures that after evaluating this expression, a garbage
#' collection run is performed, even if the expression throws an
#' error.
#'
#' I have found this function occasionally useful when dealing with
#' very large objects that bump up against the memory capacity of the
#' computer I am using. One common use case is loading a very large R
#' data file and keeping only one object from it while discarding the
#' rest. However, overuse of this function when it is not needed will
#' simply slow down your code.
#'
#' @param expr The expression to evaluate.
#' @param ... Further arguments are passed to [gc()].
#'
#' @examples
#'
#' with_gc({
#'   # Create a "large" object (this one is not actually large because
#'   # examples need to be kept small)
#'   large_object <- 1:5000
#'   # Return only a small piece of the object.
#'   large_object[5]
#' })
#' # large_object has now been garbage-collected
#'
#' @seealso [gc()]
#'
#' @export
with_gc <- function(expr, ...) {
    on.exit(gc(...))
    return(expr)
}

#' Print a readable summary of a list of values.
#'
#' This is useful for printing out a list of the parsed command-line
#' arguments for a script. The output is generally more compact than
#' [print()] or [str()].
#'
#' @param v A named list or vector
#'
#' @return Returns `v` invisibly.
#'
#' @examples
#'
#' x <- list(verbose = TRUE, infile = "a.txt", outfile = "b.txt", ncores = 8)
#' print_var_vector(x)
#'
#' @importFrom rlang is_named
#' @export
print_var_vector <- function(v) {
    assert_that(is_named(v))
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep = "")
    }
    invisible(v)
}

#' Construct futures for a list of expressions
#'
#' Instead of evaluating the expressions in the list, this function
#' constructs a future for each expression and returns the list of
#' these futures.
#'
#' This is the standard-evaluation helper for [make_futures()], which
#' uses non-standard evaluation.
#'
#' @param expressions The list of unevaluated expressions for which to
#'     construct futures. Each element can be anything that
#'     [rlang::as_quosure()] will accept.
#' @param .future.args Additional arguments to pass to
#'     [future::future()].
#'
#' @return A list of futures, whose values will be the result of
#'     evaluating each of `expressions`.
#'
#' @seealso [rlang::eval_tidy()], [rlang::quo()], [rlang::as_quosure()], [future::value()]
#'
#' @examples
#'
#' library(future)
#' library(rlang)
#' expressions <- list(a = quo(1+1), b = quo(2+2))
#'
#' flist <- make_futures_(expressions)
#' flist
#' value(flist)
#'
#' # Same result without futures
#' lapply(expressions, eval_tidy)
#'
#' @importFrom rlang eval_tidy quo as_quosure as_list
#' @export
make_futures_ <- function(expressions, .future.args = list()) {
    .future.args <- as_list(.future.args)
    if (! "lazy" %in% names(.future.args)) {
        .future.args$lazy <- TRUE
    }
    req_ns("future")
    lapply(expressions, . %>% as_quosure %>% {quo(future::future(!!., !!!.future.args))} %>% eval_tidy)
}

#' Construct futures for several expressions
#'
#' Instead of evaluating its arguments, this function converts them
#' into futures using [future::future()]. Calling [future::value()] on
#' each one will yield the value of that expression.
#'
#' @param ... Expressions for which to construct futures.
#' @param .future.args Additional arguments to pass to
#'     [future::future()].
#' @return A list of futures, one for each expression in `...`.
#'
#' @examples
#'
#' library(future)
#' flist <- make_futures(a = 1+1, b = 2+2)
#' flist
#' value(flist)
#'
#' @importFrom rlang enquos
#' @export
make_futures <- function(..., .future.args = list()) {
    make_futures_(enquos(...), .future.args)
}
