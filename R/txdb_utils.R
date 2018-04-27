#' Get a TxDb from either a package name or file name
#'
#' @export
get_txdb <- function(txdbname) {
    tryCatch({
        return(get(txdbname, loadNamespace(txdbname)))
        ## pos <- str_c("package:", txdbname)
        ## get(txdbname, pos)
    }, error=function(...) {
        req_ns("AnnotationDbi")
        AnnotationDbi::loadDb(txdbname)
    })
}

#' Generate a tx2gene table from a TxDb
#'
#' @export
get_tx2gene_from_txdb <- function(txdb) {
    req_ns("AnnotationDbi")
    k <- AnnotationDbi::keys(txdb, keytype = "GENEID")
    suppressMessages(AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")) %>%
        .[c("TXNAME", "GENEID")]
}
