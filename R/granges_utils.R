#' Compute the midpoints of a collection of read fragments
#'
#' @param reads A [GenomicAlignments::GAlignments],
#'     [GenomicAlignments::GAlignmentPairs], or
#'     [GenomicAlignments::GAlignmentsList] object representing a set
#'     of reads or read pairs.
#' @param fraglen The estimated fragment length. This is only required
#'     if `reads` contains any un-paired reads, since the far end of
#'     the fragment represented by a single-end read is unknown.
#' @param extend_fragment_upstream If FALSE (the default), each
#'     single-end read is assumed to represent the 5-prime end of a
#'     fragment, so the fragment will be extended `fraglen` base pairs
#'     in the 3-prime direction starting from the 5-prime end of the
#'     read. If reads represent the 3-prime ends of fragments, set
#'     this argument to TRUE, and the reverse will happen: fragments
#'     will be extended `fraglen` base pairs in the 5-prime direction
#'     starting from the 3-prime end of the read.
#'
#' This should really be an S4 method, but writing S4 methods is a
#' pain.
#'
#' @include internal.R
#' @importFrom magrittr %>% %<>%
#' @export
readsToFragmentMidpoints <- function(reads, fraglen, extend_fragment_upstream=FALSE) {
    req_ns("GenomicRanges")
    ## Depending on whether reads are paired or single, get the full
    ## fragment length
    if (is(reads, "GAlignmentsList")) {
        ## Split into single and paired reads
        isSingle <- lengths(reads) == 1
        mated.frags <- GenomicRanges::granges(reads[!isSingle], ignore.strand=TRUE)
        if (any(isSingle)) {
            single.frags <- GenomicRanges::granges(reads[isSingle]) %>%
                ## Extend smaller single reads to fraglen, but don't
                ## shrink longer reads
                GenomicRanges::resize(
                    width=pmax(fraglen, width(.)),
                    fix=ifelse(extend_fragment_upstream, "end", "start"))
        } else {
            single.frags <- GRanges()
        }
        frags <- c(mated.frags, single.frags)
    } else if (is(reads, "GAlignmentPairs")) {
        ## Only paired reads
        frags <- GenomicRanges::granges(reads)
    } else if (is(reads, "GAlignments")) {
        ## Only single reads
        frags <- reads %>% GenomicRanges::granges %>%
            ## Extend smaller reads to fraglen, but don't shrink
            ## longer reads
            GenomicRanges::resize(
                width=pmax(fraglen, GenomicRanges::width(.)),
                fix=ifelse(extend_fragment_upstream, "end", "start"))
    } else {
        warning(glue("Unknown reads type: {class(reads)}. Attempting to coerce to GRanges."))
        reads %<>% as("GRanges")
        ## If all ranges are the same width, assume they represent
        ## single-end reads and resize them to fraglen
        if (all(GenomicRanges::width(reads) == GenomicRanges::width(reads[1])) &&
            GenomicRanges::width(reads[1]) < fraglen) {
            warning(glue("All reads from unknown class have the same length, {width(reads[1])}, and are therefore assumed to be single-end reads, which will be resized to {fraglen}."))
            frags <- GenomicRanges::resize(reads, width=fraglen, fix="start")
        } else {
            frags <- reads
        }
    }
    ## Finally, get the center of each fragment
    GenomicRanges::resize(frags, width=1, fix="center")
}

## Convert strand to -1, 0, or 1
#' @export
strand_sign <- function(x, allow.unstranded=FALSE) {
    s <- strand(x)
    ss <- (s == "+") - (s == "-")
    if (allow.unstranded) {
        ss[ss == 0] <- NA
    } else if (any(unlist(ss == 0))) {
        stop("Strand must be either '+' or '-'")
    }
    ss
}

## This merges exons into genes (GRanges to GRangesList)
#' @export
gff.to.grl <- function(gr, exonFeatureType="exon", geneIdAttr="gene_id", geneFeatureType="gene") {
    exon.gr <- gr[gr$type %in% exonFeatureType]
    exon.gr %<>% cleanup.mcols
    grl <- split(exon.gr, as.character(mcols(exon.gr)[[geneIdAttr]])) %>%
        promote.common.mcols
    if (!is.null(geneFeatureType)) {
        gene.meta <- gr[gr$type %in% geneFeatureType] %>%
            mcols %>% cleanup.mcols(mcols_df=.) %>% .[match(names(grl), .[[geneIdAttr]]),]
        for (i in names(gene.meta)) {
            if (i %in% names(mcols(grl))) {
                value <- ifelse(is.na(gene.meta[[i]]), mcols(grl)[[i]], gene.meta[[i]])
            } else {
                value <- gene.meta[[i]]
            }
            mcols(grl)[[i]] <- value
        }
    }
    return(grl)
}

## This converts a GRangesList into the SAF ("Simplified annotation
## format")
#' @export
grl.to.saf <- function(grl) {
    gr <- unlist(grl)
    data.frame(Chr=as.vector(seqnames(gr)),
               Start=start(gr),
               End=end(gr),
               Strand=as.vector(strand(gr)),
               GeneID=rep(names(grl), lengths(grl)))
}

## Get column names that are always the same for all elements of a
## gene. Used for extracting only the gene metadata from exon
## metadata.
get.gene.common.colnames <- function(df, geneids, blacklist=c("type", "Parent")) {
    if (nrow(df) < 1) {
        return(character(0))
    }
    if (any(is.na(geneids))) {
        stop("Gene IDs cannot be undefined")
    }
    if (any(lengths(geneids) > 1)) {
        stop("Gene IDs must not be a list")
    }
    if (!anyDuplicated(geneids)) {
        return(names(df))
    }
    ## Forget blacklisted columns
    df <- df[setdiff(names(df), blacklist)]
    ## Forget list columns
    df <- df[sapply(df, . %>% lengths %>% max) == 1]
    ## Forget empty columns
    df <- df[!sapply(df, is.empty)]
    if (ncol(df) < 1) {
        return(character(0))
    }
    ## Convert to Rle
    df <- DataFrame(lapply(df, . %>% unlist %>% Rle))
    geneids %<>% Rle
    genecols <- sapply(df, . %>% split(geneids) %>% runLength %>% lengths %>% max %>% is_weakly_less_than(1))
    names(which(genecols))
}

#' Promote common mcols from the GRanges inside a GRangesList
#'
#' Given a GRangesList whose underlying ranges have mcols, this finds
#' mcols of the ranges that are constant within each gene and promotes
#' them to mcols of the GRangesList. For example, you have a
#' GRangesList of exons grouped by gene, and the exons are annotated
#' with gene IDs, then the gene ID column will be promoted to the
#' GRangesList object itself.
#'
#' @include internal.R
#' @importFrom magrittr %<>%
#' @export
promote.common.mcols <- function(grl, delete.from.source=FALSE, ...) {
    req_ns("S4Vectors", "GenomicRanges")
    colnames.to.promote <- get.gene.common.colnames(S4Vectors::mcols(unlist(grl)), rep(names(grl), lengths(grl)), ...)
    promoted.df <- S4Vectors::mcols(unlist(grl))[cumsum(lengths(grl)),colnames.to.promote, drop=FALSE]
    if (delete.from.source) {
        mcols(grl@unlistData) %<>% .[setdiff(names(.), colnames.to.promote)]
    }
    S4Vectors::mcols(grl) %<>% cbind(promoted.df)
    grl
}

#' Like rtracklayer::liftOver but "fills in" small gaps.
#'
#' @param x,chain,... These arguments have the same meaning as in
#'     [rtracklayer::liftOver()].
#' @param allow.gap Maximum gap size to "fill in". If set to zero,
#'     this is equivalent to [rtracklayer::liftOver()].
#'
#' @include internal.R
#' @export
liftOverLax <- function(x, chain, ..., allow.gap=0) {
    req_ns("rtracklayer", "S4Vectors")

    newx <- rtracklayer::liftOver(x, chain)
    if (allow.gap > 0) {
        gapped <- which(lengths(newx) > 1)
        newx.gapped.reduced <- reduce(newx[gapped], min.gapwidth = allow.gap + 1, with.revmap=TRUE)
        S4Vectors::mcols(newx.gapped.reduced@unlistData) <-
            rep(S4Vectors::mcols(x[gapped]), lengths(newx.gapped.reduced))
        newx[gapped] <- newx.gapped.reduced
    }
    return(newx)
}

#' Convenience function for running liftOver on a MotifMap BED file
#'
#' @include internal.R file_format_io.R
#' @export
liftOver_motifMap <- function(infile, chainfile, outfile, allow.gap=2) {
    req_ns("rtracklayer")
    gr <- read.motifmap(infile, parse_name=TRUE)
    chain <- rtracklayer::import.chain(chainfile)
    gr2 <- rtracklayer::liftOverLax(gr, chain, allow.gap=allow.gap)
    gr2 <- unlist(gr2[lengths(gr2) == 1])
    write.motifmap(gr2, outfile)
}
