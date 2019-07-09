## Read a table from a R data file, csv, or xlsx file. Returns a data
## frame or throws an error.

#' General function for reading a table from any table-like file
#'
#' @export
read_table_general <- function(filename, read.table.args = NULL, read.xlsx.args = NULL,
                               dataframe.class = "data.frame") {
    req_ns("future", "openxlsx")
    suppressWarnings({
        read.table.args %<>% as.list
        read.table.args$file <- filename
        read.table.args$header <- TRUE
        read.xlsx.args %<>% as.list
        read.xlsx.args$xlsxFile <- filename
        lazy_results <- make_futures(
            .future.args = list(evaluator = future::sequential, lazy = TRUE),
            rdata = read_RDS_or_RDA(filename, dataframe.class),
            table = do.call(read.table, read.table.args),
            csv = do.call(read.csv, read.table.args),
            xlsx = do.call(openxlsx::read.xlsx, read.xlsx.args))
        for (lzresult in lazy_results) {
            result <- tryCatch({
                x <- as(future::value(lzresult), dataframe.class)
                assert_that(is(x, dataframe.class))
                x
            }, error = function(...) NULL)
            if (!is.null(result)) {
                return(result)
            }
        }
        stop(glue("Could not read a data frame from {deparse{filename}} as R data, csv, or xlsx"))
    })
}

#' @export
read_idr_table <- function(file) {
    idrcols <- c("chr", "start", "end", "name", "score", "strand",
                 "LocalIDR", "GlobalIDR", "startA", "endA", "scoreA", "startB", "endB", "scoreB")
    read.table(file, header = FALSE, sep = "\t", col.names = idrcols) %>%
        mutate_(LocalIDR = ~10^-LocalIDR, GlobalIDR = ~10^-GlobalIDR)
}

#' Read MotifMap-provided BED file into a GRanges object.
#'
#' Even though this function reads a BED file,
#' [rtracklayer::import.bed()] is not suitable because it chokes on
#' spaces in fields, which MotifMap contains. Hence the need for a
#' separate function.
#'
#' @param file The file name to read from.
#' @param parse_name If TRUE, separate the "name" column into
#'     "motif_ID" and "TF_name" fields.
#' @return a GRanges object with an mcol named "score". If
#'     `parse_name` is FALSE, it will also have an mcol named "name".
#'     If `parse_name` is TRUE, it will have mcols named "motif_ID",
#'     and "TF_name".
#'
#' @export
read_motifmap <- function(file, parse_name = TRUE) {
    req_ns("tidyr", "GenomicRanges", "readr")
    tab <- readr::read_tsv(file, col_names = c("chr", "start", "end", "name", "score", "strand"),
                    col_types = "ciicdc", progress = FALSE)
    if (parse_name) {
        tab %<>% tidyr::separate_(~name, into = c("motif_ID", "TF_name"), sep = " = ")
    }
    gr <- GenomicRanges::makeGRangesFromDataFrame(
        tab, starts.in.df.are.0based = TRUE,
        keep.extra.columns = TRUE)
    gr
}

#' Write a MotifMap-type GRanges into a BED file.
#'
#' This takes a GRanges object with the appropriate mcols (see
#' [read_motifmap()]) and saves it as a BED file.
#'
#' @param x A GRanges object with appropriate metadata columns.
#' @param file The file name to save.
#'
#' @export
write_motifmap <- function(x, file) {
    req_ns("S4Vectors", "rtracklayer", "tidyr")
    assert_that(is(x, "GRanges"))
    if (! "name" %in% names(S4Vectors::mcols(x))) {
        assert_that(all(c("motif_ID", "TF_name") %in% names(S4Vectors::mcols(x))))
        S4Vectors::mcols(x) %<>% as.data.frame %>%
            tidyr::unite_(~name, ~c(motif_ID, TF_name), sep = " = ") %>%
            as("DataFrame")
    }
    rtracklayer::export(x, file, format = "BED")
}

#' Read a narrowPeak format BED file.
#'
#' NarrowPeak is a sub-format of the BED file format used as the
#' output format of some peak callers. [rtracklayer::import.bed()]
#' often throws an error when reading these files, so a separate
#' function is needed.
#'
#' @param file The file name to read from.
#' @param ... Additional arguments are passed to `read.table()`.
#' @return A GRanges object containing the ranges from `file`.
#'
#' @export
read_narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep = "\t", row.names = NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$name <- as.character(peaks.df$name)
    GenomicRanges::makeGRangesFromDataFrame(
        peaks.df, starts.in.df.are.0based = TRUE,
        keep.extra.columns = TRUE)
}

#' Write a narrowPeak format BED file.
#'
#' This function takes a GRanges of the kind returned by
#' [read_narrowPeak()] and writes it back to a file.
#'
#' @param x A GRanges or data frame containing the all appropriate
#'     columns for a narrowPeak file. In addition to the normal BED
#'     columns, it must contain columns named "signalValue", "pValue",
#'     "qValue", and "summit".
#' @param file The file name to write to.
#' @param ... Additional arguments are passed to `write.table()`.
#'
#' If some of the required columns are not available (e.g. from a peak
#' caller that does not report a summit), you can fill in the
#' unavailable columns with NA.
#'
#' @export
write_narrowPeak <- function(x, file, ...) {
    x <- as(x, "data.frame")
    if("seqnames" %in% names(x))
        names(x)[names(x) == "seqnames"] <- "chr"
    x <- x[c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")]
    write.table(x, file, sep = "\t", row.names = FALSE, col.names = FALSE, ...)
}

#' Read a GRanges from a variety of possible file types.
#'
#' This function attempts to read in a set of genomic regions as a
#' GRanges object. It does its best to auto-detect the input file
#' type. If it cannot do so, it throws an error.
#'
#' @param filename The file to read regions from.
#' @return a GRanges object.
#'
#' Possible input file types include: RDA or RDS file (containing a
#' GRanges or data frame), narrowPeak, BED, GFF, and SAF.
#'
#' @export
read_regions <- function(filename) {
    req_ns("rtracklayer", "future")
    suppressWarnings({
        lazy_results <- make_futures(
            .future.args = list(evaluator = future::sequential, lazy = TRUE),
            rdata = read_RDS_or_RDA(filename),
            narrowPeak = read_narrowPeak(filename),
            bed = rtracklayer::import(filename, format = "bed"),
            gff = rtracklayer::import(filename, format = "gff"),
            saf = read_saf(filename),
            table = read_table_general(filename))
        for (lzresult in lazy_results) {
            result <- tryCatch({
                x <- future::value(lzresult)
                if (is(x, "List")) {
                    x <- unlist(x)
                }
                x <- as(x, "GRanges")
                x
            }, error = function(...) NULL)
            if (!is.null(result)) {
                return(result)
            }
        }
        stop(glue("Could not read genomic regions from {deparse(filename)} as R data, narrowPeak, bed, gff, SAF, or csv"))
    })
}

## TODO: Make import.saf and export.saf functions?

#' Read a SAF file into a GRangesList.
#'
#' The SAF format is described in [Rsubread::featureCounts()].
#'
#' @param filename The file name to read
#' @param ... Additional arguments are passed to
#'     [read_table_general()].
#' @return a GRangesList with one element for each GeneID, containing
#'     the ranges for that gene.
#'
#' @export
read_saf <- function(filename, ...) {
    saf <- read_table_general(filename, ...)
    assert_that("GeneID" %in% names(saf))
    gr <- as(saf, "GRanges")
    grl <- split(gr, gr$GeneID) %>%
        promote_common_mcols(delete_from_source = FALSE)
    return(grl)
}

## TODO: Give the below functions better names

#' Read a Salmon "geneMap" file into a data frame
#'
#' This function reads a file in the Salmon "geneMap" format into a
#' data frame. The format of this file is described in salmon's help
#' text for the "--geneMap" command-line option.
#'
#' @param fname The file to read from. This should generally be a file
#'     named "genemap.txt" in a Salmon index directory
#' @return A data frame with 2 columns named "TXNAME" and "GENEID".
#'
#' Since this uses [read_table_general()] and allows additional
#' columns to be present beyond the first two, it is actually quite
#' lenient about the format of the table.
#'
#' @seealso The "--geneMap" option to salmon, whose help text
#'     describes the format of the file.
#'
#' @export
read_tx2gene_from_genemap <- function(fname) {
    df <- read_table_general(
        fname,
        read.xlsx.args = list(colNames = FALSE),
        read.table.args = list(header = FALSE))
    df %<>% .[1:2]
    df[] %<>% lapply(as.character)
    names(df) <- c("TXNAME", "GENEID")
    df
}

#' @export
read_annotation_from_gff <- function(filename, format = "GFF3", ...) {
    req_ns("rtracklayer")
    gff <- NULL
    ## Allow the file to be an RDS file containing the GRanges
    ## resulting from import()
    gff <- tryCatch({
        read_RDS_or_RDA(filename, "GRanges")
    }, error = function(...) {
        req_ns("rtracklayer")
        rtracklayer::import(filename, format = format)
    })
    assert_that(is(gff, "GRanges"))
    grl <- gff_to_grl(gff, ...)
    return(grl)
}

#' @export
read_annotation_from_saf <- function(filename, ...) {
    saf <- read_table_general(filename, ...)
    assert_that("GeneID" %in% names(saf))
    gr <- as(saf, "GRanges")
    grl <- split(gr, gr$GeneID) %>% promote_common_mcols
    return(grl)
}

#' @export
read_annotation_from_rdata <- function(filename) {
    read_RDS_or_RDA(filename, "GRangesList")
}

#' @export
read_additional_gene_info <- function(filename, gff_format = "GFF3",
                                      geneFeatureType = "gene", ...) {
    req_ns("S4Vectors", "rtracklayer")
    df <- tryCatch({
        gff <- tryCatch({
            read_RDS_or_RDA(filename, "GRanges")
        }, error = function(...) {
            req_ns("rtracklayer")
            rtracklayer::import(filename, format = gff_format)
        })
        assert_that(is(gff, "GRanges"))
        if (!is.null(geneFeatureType)) {
            gff %<>% .[.$type %in% geneFeatureType]
        }
        gff %<>% .[!is.na(.$ID) & !duplicated(.$ID)]
        gff %>% S4Vectors::mcols %>% cleanup_mcols(mcols_df = .)
    }, error = function(...) {
        tab <- read_table_general(filename, ..., dataframe.class = "DataFrame")
        ## Nonexistent or automatic row names
        if (.row_names_info(tab) <= 0) {
            row.names(tab) <- tab[[1]]
        }
        tab
    })
    df %<>% S4Vectors::DataFrame
    assert_that(is(df, "DataFrame"))
    return(df)
}
