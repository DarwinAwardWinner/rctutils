# Make cairo_pdf use onefile = TRUE by default
## TODO: Rename to something unique

#' Same as `cairo_pdf()` but with a default of `onefile = TRUE`
#'
#' @export
cairo_pdf_onefile <- function(..., onefile = TRUE) {
    cairo_pdf(..., onefile = onefile)
}

## Note: relies on external command-line util. TODO: Use
## https://cran.r-project.org/web/packages/magick/vignettes/intro.html

#' @export
rasterpdf <- function(pdffile, outfile = pdffile, resolution = 600) {
    tempf <- tempfile(pattern = "raster", fileext = ".pdf")
    on.exit(unlink(tempf))
    exitcode <- system2("convert", args = c("-density", resolution, pdffile, tempf),
        stdout = FALSE, stderr = FALSE)
    assert_that(exitcode == 0)
    assert_that(file.exists(tempf))
    suppressWarnings(file.rename(tempf, outfile))
    # If file still exists, then the rename failed because it's a
    # cross-device move, so copy and delete instead.
    if (file.exists(tempf)) {
        file.copy(tempf, outfile)
    }
}

## Returns TRUE if x refers to the device number of a currently active
## graphics device.
#' @importFrom rlang is_scalar_integer
#' @export
is_dev <- function(x) {
    is_scalar_integer(x) && x %in% dev.list()
}

#' @export
with_dev <- function(dev, code, closedev) {
    orig.device <- dev.cur()
    new.device <- force(dev)
    # Functions that create devices don't generally return them, they
    # just set them as the new current device, so get the actual
    # device from dev.cur() instead.
    if (is.null(new.device)) {
        new.device <- dev.cur()
    }
    assert_that(is_dev(new.device) || new.device == 1)
    if (missing(closedev)) {
         closedev <- new.device != orig.device
    }
    on.exit({
        if (closedev) {
            dev.off(new.device)
        }
        if (is_dev(orig.device)) {
            dev.set(orig.device)
        }
    })
    force(code)
}

# Useful to wrap functions that both produce a plot and return a
# useful value, when you only want the return value and not the plot.
#' @export
suppressPlot <- function(arg) {
    png("/dev/null")
    result <- arg
    dev.off()
    result
}
