## Read a table from a R data file, csv, or xlsx file. Returns a data
## frame or throws an error.
read.table.general <- function(filename, read.table.args=NULL, read.xlsx.args=NULL,
                               dataframe.class="data.frame") {
    suppressWarnings({
        read.table.args %<>% as.list
        read.table.args$file <- filename
        read.table.args$header <- TRUE
        read.xlsx.args %<>% as.list
        read.xlsx.args$xlsxFile <- filename
        lazy.results <- list(
            rdata=future(read.RDS.or.RDA(filename, dataframe.class), lazy=TRUE),
            table=future(do.call(read.table, read.table.args), lazy=TRUE),
            csv=future(do.call(read.csv, read.table.args), lazy=TRUE),
            xlsx=future(do.call(read.xlsx, read.xlsx.args), lazy=TRUE))
        for (lzresult in lazy.results) {
            result <- tryCatch({
                x <- as(value(lzresult), dataframe.class)
                assert_that(is(x, dataframe.class))
                x
            }, error=function(...) NULL)
            if (!is.null(result)) {
                return(result)
            }
        }
        stop(glue("Could not read a data frame from {deparse{filename}} as R data, csv, or xlsx"))
    })
}

read.idr.table <- function(file) {
    idrcols <- c("chr", "start", "end", "name", "score", "strand",
                 "LocalIDR", "GlobalIDR", "startA", "endA", "scoreA", "startB", "endB", "scoreB")
    read.table(file, header=FALSE, sep="\t", col.names=idrcols) %>%
        mutate(LocalIDR=10^-LocalIDR, GlobalIDR=10^-GlobalIDR)
}

## Read MotifMap-provided BED file into a GRanges object. We can't use
## rtracklayer::import.bed because it chokes on spaces in fields,
## which MotifMap contains.
read.motifmap <- function(x, parse_name=TRUE) {
    tab <- read_tsv(x, col_names=c("chr", "start", "end", "name", "score", "strand"),
                    col_types="ciicdc", progress=FALSE)
    if (parse_name) {
        tab %<>% separate(name, into=c("motif_ID", "TF_name"), sep="=")
    }
    gr <- makeGRangesFromDataFrame(tab, starts.in.df.are.0based=TRUE)
    gr
}

## write.motifmap

read.narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep="\t", row.names=NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$name <- as.character(peaks.df$name)
    peaks.df
}

write.narrowPeak <- function(x, file, ...) {
    x <- as(x, "data.frame")
    if("seqnames" %in% names(x))
        names(x)[names(x) == "seqnames"] <- "chr"
    x <- x[c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")]
    write.table(x, file, sep="\t", row.names=FALSE, col.names=FALSE, ...)
}

read.regions <- function(filename) {
    suppressWarnings({
        lazy.results <- list(
            rdata=future(read.RDS.or.RDA(filename), lazy=TRUE),
            narrowPeak=future(read.narrowPeak(filename), lazy=TRUE),
            bed=future(import(filename, format="bed"), lazy=TRUE),
            gff=future(import(filename, format="gff"), lazy=TRUE),
            saf=future(read.saf(filename), lazy=TRUE),
            table=future(read.table.general(filename), lazy=TRUE))
        for (lzresult in lazy.results) {
            result <- tryCatch({
                x <- value(lzresult)
                if (is(x, "List")) {
                    x <- unlist(x)
                }
                x <- as(x, "GRanges")
                x
            }, error=function(...) NULL)
            if (!is.null(result)) {
                return(result)
            }
        }
        stop(glue("Could not read genomic regions from {deparse(filename)} as R data, narrowPeak, bed, gff, SAF, or csv"))
    })
}

read.saf <- function(filename, ...) {
    saf <- read.table.general(filename, ...)
    assert_that("GeneID" %in% names(saf))
    gr <- as(saf, "GRanges")
    grl <- split(gr, gr$GeneID) %>% promote.common.mcols
    return(grl)
}

## TODO: Give the below functions better names

read.tx2gene.from.genemap <- function(fname) {
    df <- read.table.general(fname)
    df %<>% .[1:2]
    df[] %<>% lapply(as.character)
    names(df) <- c("TXNAME", "GENEID")
    df
}

read.annotation.from.gff <- function(filename, format="GFF3", ...) {
    gff <- NULL
    ## Allow the file to be an RDS file containing the GRanges
    ## resulting from import()
    gff <- tryCatch({
        read.RDS.or.RDA(filename, "GRanges")
    }, error=function(...) {
        import(filename, format=format)
    })
    assert_that(is(gff, "GRanges"))
    grl <- gff.to.grl(gff, ...)
    return(grl)
}

read.annotation.from.saf <- function(filename, ...) {
    saf <- read.table.general(filename, ...)
    assert_that("GeneID" %in% names(saf))
    gr <- as(saf, "GRanges")
    grl <- split(gr, gr$GeneID) %>% promote.common.mcols
    return(grl)
}

read.annotation.from.rdata <- function(filename) {
    read.RDS.or.RDA(filename, "GRangesList")
}

read.additional.gene.info <- function(filename, gff_format="GFF3", geneFeatureType="gene",
 ...) {
    df <- tryCatch({
        gff <- tryCatch({
            read.RDS.or.RDA(filename, "GRanges")
        }, error=function(...) {
            import(filename, format=gff_format)
        })
        assert_that(is(gff, "GRanges"))
        if (!is.null(geneFeatureType)) {
            gff %<>% .[.$type %in% geneFeatureType]
        }
        gff %<>% .[!is.na(.$ID) & !duplicated(.$ID)]
        gff %>% mcols %>% cleanup.mcols(mcols_df=.)
    }, error=function(...) {
        tab <- read.table.general(filename, ..., dataframe.class="DataFrame")
        ## Nonexistent or automatic row names
        if (.row_names_info(tab) <= 0) {
            row.names(tab) <- tab[[1]]
        }
        tab
    })
    df %<>% DataFrame
    assert_that(is(df, "DataFrame"))
    return(df)
}
