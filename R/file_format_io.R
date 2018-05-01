## Read a table from a R data file, csv, or xlsx file. Returns a data
## frame or throws an error.

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
#' We can't use rtracklayer::import.bed because it chokes on spaces in
#' fields, which MotifMap contains.
#'
#' @export
read_motifmap <- function(x, parse_name = TRUE) {
    req_ns("tidyr", "GenomicRanges", "readr")
    tab <- readr::read_tsv(x, col_names = c("chr", "start", "end", "name", "score", "strand"),
                    col_types = "ciicdc", progress = FALSE)
    if (parse_name) {
        tab %<>% tidyr::separate_(~name, into = c("motif_ID", "TF_name"), sep = " = ")
    }
    gr <- GenomicRanges::makeGRangesFromDataFrame(tab, starts.in.df.are.0based = TRUE)
    gr
}

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

#' @export
read_narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep = "\t", row.names = NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$name <- as.character(peaks.df$name)
    peaks.df
}

#' @export
write_narrowPeak <- function(x, file, ...) {
    x <- as(x, "data.frame")
    if("seqnames" %in% names(x))
        names(x)[names(x) == "seqnames"] <- "chr"
    x <- x[c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")]
    write.table(x, file, sep = "\t", row.names = FALSE, col.names = FALSE, ...)
}

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

#' @export
read_saf <- function(filename, ...) {
    saf <- read_table_general(filename, ...)
    assert_that("GeneID" %in% names(saf))
    gr <- as(saf, "GRanges")
    grl <- split(gr, gr$GeneID) %>% promote_common_mcols
    return(grl)
}

## TODO: Give the below functions better names

#' @export
read_tx2gene_from_genemap <- function(fname) {
    df <- read_table_general(fname)
    df %<>% .[1:2]
    df[] %<>% lapply(as.character)
    names(df) <- c("TXNAME", "GENEID")
    df
}

#' @export
read_annotation_from_gff <- function(filename, format = "GFF3", ...) {
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
    req_ns("S4Vectors")
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
