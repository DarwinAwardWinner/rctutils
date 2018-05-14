#' Variant of `save.image()` that allows excluding specific names
#'
#' If you want to save your entire workspace except for a specific set
#' of variable names that you wish to exclude, this function makes it
#' easy to do so.
#'
#' @param file,version,ascii,compress,safe These arguments have the
#'     same meaning as in `save.image()`.
#' @param exclude Character vector of variable names not to save.
#'
#' @examples
#'
#' x <- 5
#' y <- 3
#' datafile <- tempfile(fileext = ".rda")
#' # Only saves y, not x
#' save_image_filtered(datafile, exclude = "x")
#'
#' @seealso [save.image()]
#'
#' @export
save_image_filtered <-
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
    if (safe) {
        if (!file.rename(outfile, file)) {
            on.exit()
            stop(gettextf("image could not be renamed and is left in %s",
                          outfile), domain = NA)
        }
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
#'
#' @examples
#'
#' x <- 5
#' y <- 3
#' datafile <- tempfile(fileext = ".rda")
#' save(file = datafile, list = c("x", "y"))
#' rm(x,y)
#' loaded_vars <- load_in_new_env(datafile)
#' as.list(loaded_vars)
#'
#' @seealso [save.image()]
#'
#' @export
load_in_new_env <- function(file, envir = new.env(), ...) {
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
#'
#' x <- 5
#' y <- 3
#' datafile <- tempfile(fileext = ".rda")
#' # Only saves y, not x
#' save.image(datafile)
#' # Set new values which can be overridden by loading the old ones.
#' x <- 10
#' y <- 10
#' load_filtered(datafile, exclude = "x")
#' # x is still 10, but y is now 3 again
#' mget(c("x", "y"))
#'
#' @seealso [load()]
#' @export
load_filtered <- function(file, envir = parent.frame(), ..., exclude = NULL) {
    if (!length(exclude)) {
        return(load(file, envir, ...))
    }
    tempenv <- load_in_new_env(file = file, ...)
    for (i in setdiff(names(tempenv), exclude)) {
        envir[[i]] <- tempenv[[i]]
    }
}

#' Read a single R object from an RDA file.
#'
#' If run on an RDA file containing more than one object, throws an
#' error.
#'
#' @param filename The file to read an R object from.
#' @return The R object read from the file.
#'
#' This function allows you to read an RData file as if it were an RDS
#' file, as long as it satisfies the requirement of containing only
#' one object.
#'
#' @seealso [readRDS()], [load()]
#'
#' @export
read_single_object_from_rda <- function(filename) {
    objects <- within(list(), suppressWarnings(load(filename)))
    if (length(objects) != 1) {
        stop("RDA file should contain exactly one object")
    }
    return(objects[[1]])
}

#' Read a single object from an RDS or RDA file
#'
#' This is like [readRDS()] except that it can also read from an RData
#' file containing a single object (i.e. the kind of file that is read
#' using [load()]). Use this in place of `readRDS()` if you want to be
#' slightly more forgiving about what kind of R data file you accept.
#'
#' @param filename The file to read an R object from.
#' @param expected.class If specified, the object will be coerced into
#'     this class using `as()`, which will throw an error as normal if
#'     the coercion is not possible. This allows you to restrict what
#'     kind of objects you will accept. This can also be a function
#'     that accepts a single argument and performs the proper coercion
#'     itself.
#' @return The object read from the file, possibly after coercing it
#'     into another class.
#'
#' @examples
#'
#' tmpf <- tempfile()
#' saveRDS(1:10, tmpf)
#' read_RDS_or_RDA(tmpf)
#' read_RDS_or_RDA(tmpf, "character")
#' # Using a function instead of a class name.
#' read_RDS_or_RDA(tmpf, as.character)
#' read_RDS_or_RDA(tmpf, "factor")
#' read_RDS_or_RDA(tmpf, "data.frame")
#'
#' \dontrun{
#' # This will throw an error because the coercion to "lm" is not
#' # possible.
#' read_RDS_or_RDA(tmpf, "lm")
#' }
#'
#' @seealso [readRDS()], [read_single_object_from_rda()], [as()]
#'
#' @importFrom glue glue
#' @export
read_RDS_or_RDA <- function(filename, expected.class = "ANY") {
    object <- suppressWarnings(tryCatch({
        readRDS(filename)
    }, error = function(...) {
        read_single_object_from_rda(filename)
    }))
    if (is.function(expected.class)) {
        object <- do.call(expected.class, list(object))
    } else if (!is(object, expected.class)) {
        ## Try to use as.[CLASS] if it exists. If not use as(.,
        ## "[CLASS]").
        coerce_func <- tryCatch(
            get(glue("as.{expected.class}")),
            error = function(...) . %>% as(expected.class))
        object <- coerce_func(object)
    }
    return(object)
}

#' Save a single object to am RDS or RDA file
#'
#' If the output file name ends in ".rda" or ".rdata" (with any
#' capitalization), the file will be saved using [save()]. Otherwise,
#' it will be saved using [saveRDS()].
#'
#' @param object R object to serialize.
#' @param file A file name or connection object to save to.
#' @param ascii,version,compress See `saveRDS()`.
#' @param savetype Either "rda" or "rds". If this argument is
#'     supplied, it overrides the automatic detection of whether to
#'     use `save()` or `saveRDS()` based on the file extension.
#'
#' @examples
#'
#' tmpf_rda <- tempfile(fileext = ".rda")
#' tmpf_rds <- tempfile(fileext = ".rds")
#' x <- 1:10
#' save_RDS_or_RDA(x, tmpf_rda)
#' save_RDS_or_RDA(x, tmpf_rds)
#' # Load as an RData file
#' x_rda <- read_single_object_from_rda(tmpf_rda)
#' # Load as an RDS file
#' x_rds <- readRDS(tmpf_rds)
#' # Both loaded objects are identical to the original
#' identical(x, x_rda)
#' identical(x, x_rds)
#'
#' @export
save_RDS_or_RDA <- function(object, file, ascii = FALSE, version = NULL,
                            compress = TRUE, savetype) {
    if (missing(savetype)) {
        if (is.character(file) && str_detect(file, regex("\\.rda(ta)?", ignore_case = TRUE))) {
            savetype <- "rda"
        } else {
            savetype <- "rds"
        }
    }
    savetype <- match_arg(savetype, choices = c("rda", "rds"))
    if (savetype == "rda") {
        save(list = "object", file = file, ascii = ascii, version = version, compress = compress)
    } else{
        saveRDS(object = object, file = file, ascii = ascii, version = version, compress = compress)
    }
}
