#' Variant of `save.image()` that allows excluding specific names
#'
#' If you want to save your entire workspace except for a specific set
#' of variable names that you wish to exclude, this function makes it
#' easy to do so.
#'
#' @param file,version,ascii,compress,safe These arguments have the
#'     same meaning as in `save.image()`.
#' @param exclude Character vector of variable names not to save.
#' @examples
#' x <- 5
#' y <- 3
#' datafile <- tempfile(fileext=".rda")
#' # Only saves y, not x
#' save.image.filtered(datafile, exclude="x")
#' @seealso [save.image()]
#' @export
save.image.filtered <-
    function (file = ".RData", version = NULL, ascii = FALSE,
              compress = !ascii, safe = TRUE, exclude = NULL)
{
    if (!is.character(file) || file == "")
        stop("'file' must be non-empty string")
    opts <- getOption("save.image.defaults")
    if (is.null(opts))
        opts <- getOption("save.defaults")
    if (missing(safe) && !is.null(opts$safe))
        safe <- opts$safe
    if (missing(ascii) && !is.null(opts$ascii))
        ascii <- opts$ascii
    if (missing(compress) && !is.null(opts$compress))
        compress <- opts$compress
    if (missing(version))
        version <- opts$version
    if (safe) {
        outfile <- paste0(file, "Tmp")
        i <- 0
        while (file.exists(outfile)) {
            i <- i + 1
            outfile <- paste0(file, "Tmp", i)
        }
    }
    else outfile <- file
    on.exit(file.remove(outfile))
    vars.to.save <- setdiff(names(.GlobalEnv), exclude)
    save(list = vars.to.save, file = outfile, version = version,
         ascii = ascii, compress = compress, envir = .GlobalEnv,
         precheck = FALSE)
    if (safe)
        if (!file.rename(outfile, file)) {
            on.exit()
            stop(gettextf("image could not be renamed and is left in %s",
                          outfile), domain = NA)
        }
    on.exit()
}

#' Load an R data file into a new environment
#'
#' Unlike `load()`, this function avoids polluting the global
#' namespace by loading objects into separate environment.
#'
#' This function is useful for loading data files that potentially
#' contain many objects with common names that might otherwise
#' overwrite existing objects with the same name.
#'
#' @param file This argument has the same meaning as in `load()`
#' @param envir If you pass an existing environment, objects will be
#'     loaded into that environment instead of a newly-created one,
#'     and this same environment will be returned.
#' @param ... Further arguments are passed to `load()`.
#' @return An environment containing all the objects loaded from
#'     `file`.
#' @examples
#' x <- 5
#' y <- 3
#' datafile <- tempfile(fileext=".rda")
#' save(file=datafile, list=c("x", "y"))
#' rm(x,y)
#' loaded_vars <- load.in.new.env(datafile)
#' as.list(loaded_vars)
#' @seealso [save.image()]
#' @export
load.in.new.env <- function(file, envir=new.env(), ...) {
    load(file, envir, ...)
    return(envir)
}

#' Variant of `load()` that allows excluding specific names
#'
#' If you want to load every object in an R data file except for a
#' specific set of variable names that you wish to exclude, this
#' function makes it easy to do so.
#'
#' @param file,envir,... These arguments have the same meaning as in
#'     `load()`.
#' @param exclude Character vector of variable names not to load.
#' @examples
#' x <- 5
#' y <- 3
#' datafile <- tempfile(fileext=".rda")
#' # Only saves y, not x
#' save.image(datafile)
#' # Set new values which can be overridden by loading the old ones.
#' x <- 10
#' y <- 10
#' load.filtered(datafile, exclude="x")
#' # x is still 10, but y is now 3 again
#' mget(c("x", "y"))
#' @seealso [load()]
#' @export
load.filtered <- function(file, envir = parent.frame(), ..., exclude=NULL) {
    if (!length(exclude)) {
        return(load(file, envir, ...))
    }
    tempenv <- load.in.new.env(file=file, ...)
    for (i in setdiff(names(tempenv), exclude)) {
        envir[[i]] <- tempenv[[i]]
    }
}

#' Read a single R object from an RDA file.
#'
#' If run on an RDA file containing more than one object, throws an
#' error.
#'
#' @export
read.single.object.from.rda <- function(filename) {
    objects <- within(list(), suppressWarnings(load(filename)))
    if (length(objects) != 1) {
        stop("RDA file should contain exactly one object")
    }
    return(objects[[1]])
}

#' Read a single object from RDS or RDA file
#'
#' @export
read.RDS.or.RDA <- function(filename, expected.class="ANY") {
    object <- suppressWarnings(tryCatch({
        readRDS(filename)
    }, error=function(...) {
        read.single.object.from.rda(filename)
    }))
    if (!any(sapply(expected.class, is, object=object))) {
        object <- as(object, expected.class)
    }
    return(object)
}
