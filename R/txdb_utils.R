#' Get a TxDb from either a package name or file name
#'
#' @export
get.txdb <- function(txdbname) {
    tryCatch({
        library(txdbname, character.only=TRUE)
        pos <- str_c("package:", txdbname)
        get(txdbname, pos)
    }, error=function(...) {
        library(GenomicFeatures)
        loadDb(txdbname)
    })
}

#' @export
get.tx2gene.from.txdb <- function(txdb) {
    k <- keys(txdb, keytype = "GENEID")
    suppressMessages(AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")) %>%
        .[c("TXNAME", "GENEID")]
}