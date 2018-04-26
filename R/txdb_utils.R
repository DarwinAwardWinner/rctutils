#' Get a TxDb from either a package name or file name
#'
#' @export
get.txdb <- function(txdbname) {
    tryCatch({
        return(get(txdbname, loadNamespace(txdbname)))
        ## pos <- str_c("package:", txdbname)
        ## get(txdbname, pos)
    }, error=function(...) {
        req_ns("AnnotationDbi")
        AnnotationDbi::loadDb(txdbname)
    })
}

#' @export
get.tx2gene.from.txdb <- function(txdb) {
    req_ns("AnnotationDbi")
    k <- AnnotationDbi::keys(txdb, keytype = "GENEID")
    suppressMessages(AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")) %>%
        .[c("TXNAME", "GENEID")]
}
