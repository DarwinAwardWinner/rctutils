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
#' @export
readsToFragmentMidpoints <- function(reads, fraglen, extend_fragment_upstream = FALSE) {
    req_ns("GenomicRanges")
    ## Depending on whether reads are paired or single, get the full
    ## fragment length
    if (is(reads, "GAlignmentsList")) {
        ## Split into single and paired reads
        isSingle <- lengths(reads) == 1
        mated.frags <- GenomicRanges::granges(reads[!isSingle], ignore.strand = TRUE)
        if (any(isSingle)) {
            single.frags <- GenomicRanges::granges(reads[isSingle]) %>%
                ## Extend smaller single reads to fraglen, but don't
                ## shrink longer reads
                GenomicRanges::resize(
                    width = pmax(fraglen, GenomicRanges::width(.)),
                    fix = ifelse(extend_fragment_upstream, "end", "start"))
        } else {
            single.frags <- GenomicRanges::GRanges()
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
                width = pmax(fraglen, GenomicRanges::width(.)),
                fix = ifelse(extend_fragment_upstream, "end", "start"))
    } else {
        warning(glue("Unknown reads type: {class(reads)}. Attempting to coerce to GRanges."))
        reads %<>% as("GRanges")
        ## If all ranges are the same width, assume they represent
        ## single-end reads and resize them to fraglen
        if (all(GenomicRanges::width(reads) == GenomicRanges::width(reads[1])) &&
            GenomicRanges::width(reads[1]) < fraglen) {
            warning(glue("All reads from unknown class have the same length, {width(reads[1])}, and are therefore assumed to be single-end reads, which will be resized to {fraglen}."))
            frags <- GenomicRanges::resize(reads, width = fraglen, fix = "start")
        } else {
            frags <- reads
        }
    }
    ## Finally, get the center of each fragment
    GenomicRanges::resize(frags, width = 1, fix = "center")
}

## Convert strand to -1, 0, or 1
#' @export
strand_sign <- function(x, allow.unstranded = FALSE) {
    req_ns("GenomicRanges")
    s <- GenomicRanges::strand(x)
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
gff_to_grl <- function(gr, exonFeatureType = "exon", geneIdAttr = "gene_id", geneFeatureType = "gene") {
    req_ns("S4Vectors")
    exon.gr <- gr[gr$type %in% exonFeatureType]
    exon.gr %<>% cleanup_mcols
    grl <- split(exon.gr, as.character(S4Vectors::mcols(exon.gr)[[geneIdAttr]])) %>%
        promote_common_mcols
    if (!is.null(geneFeatureType)) {
        gene.meta <- gr[gr$type %in% geneFeatureType] %>%
            S4Vectors::mcols %>% cleanup_mcols(mcols_df = .) %>%
            .[match(names(grl), .[[geneIdAttr]]),]
        for (i in names(gene.meta)) {
            if (i %in% names(S4Vectors::mcols(grl))) {
                value <- ifelse(is.na(gene.meta[[i]]), S4Vectors::mcols(grl)[[i]], gene.meta[[i]])
            } else {
                value <- gene.meta[[i]]
            }
            S4Vectors::mcols(grl)[[i]] <- value
        }
    }
    return(grl)
}

## This converts a GRangesList into the SAF ("Simplified annotation
## format")
#' @export
grl_to_saf <- function(grl) {
    req_ns("GenomicRanges")
    gr <- unlist(grl)
    data.frame(Chr = as.vector(GenomicRanges::seqnames(gr)),
               Start = GenomicRanges::start(gr),
               End = GenomicRanges::end(gr),
               Strand = as.vector(GenomicRanges::strand(gr)),
               GeneID = rep(names(grl), lengths(grl)))
}

## Get column names that are always the same for all elements of a
## gene. Used for extracting only the gene metadata from exon
## metadata.

get_gene_common_colnames <- function(df, geneids) {
    req_ns("S4Vectors")
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
    ## Forget list columns
    df <- df[sapply(df, . %>% lengths %>% max) == 1]
    ## Forget empty columns
    df <- df[!sapply(df, is_valueless)]
    if (ncol(df) < 1) {
        return(character(0))
    }
    ## Convert to Rle
    df <- S4Vectors::DataFrame(lapply(df, function(x) S4Vectors::Rle(unlist(x))))
    geneids <- S4Vectors::Rle(geneids)
    genecols <- sapply(df, . %>% split(geneids) %>% as("List") %>% (S4Vectors::runLength) %>% lengths %>% max %>% is_weakly_less_than(1))
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
#' @param grl A GRangesList
#' @param delete_from_source If TRUE, delete any promoted mcols from
#'     the underlying GRanges, so that the GRangesList and underlying
#'     GRanges will not have any mcols in common (unless they already
#'     had some before).
#' @param blacklist A vector of column names that should never be
#'     promoted. The default is to never promote "type" or "Parent",
#'     which is appropriate for typical GFF files.
#' @return The same GRangesList, possibly with additional mcols
#'     derived from the underlying GRanges.
#'
#' Note that if the GRangesList's mcols already contain some of the
#' same colnames as the mcols of the inner GRanges, these will be
#' replaced by the promoted mcols.
#'
#' @examples
#'
#' gr <- GRanges("chr1", IRanges(start = (1:10) * 100, width = 50),
#'               GeneID = rep(c("geneA", "geneB"), each=5),
#'               GeneName = rep(c("Gene A", "Gene B"), each=5),
#'               ExonID = rep(1:5, 2))
#' grl <- split(gr, gr$GeneID)
#' # All the mcols are still in the underlying GRanges, with none in
#' # the GRangesList.
#' names(mcols(grl))
#' names(mcols(unlist(grl)))
#'
#' grl <- promote_common_mcols(grl)
#' # Now grl has received the common mcols
#' names(mcols(grl))
#' names(mcols(unlist(grl)))
#'
#' grl <- promote_common_mcols(grl, delete_from_source = TRUE)
#' # Now grl and only grl has the common mcols
#' names(mcols(grl))
#' names(mcols(unlist(grl)))
#'
#' @export
promote_common_mcols <- function(grl, delete_from_source = FALSE, blacklist = c("type", "Parent")) {
    req_ns("S4Vectors", "GenomicRanges")
    colnames.to.promote <- get_gene_common_colnames(S4Vectors::mcols(unlist(grl)), rep(names(grl), lengths(grl))) %>%
        setdiff(blacklist)
    promoted.df <- S4Vectors::mcols(unlist(grl))[cumsum(lengths(grl)),colnames.to.promote, drop = FALSE]
    if (delete_from_source) {
        S4Vectors::mcols(grl@unlistData) %<>%
            .[setdiff(names(.), colnames.to.promote)]
    }
    S4Vectors::mcols(grl)[colnames(promoted.df)] <- promoted.df
    grl
}

#' Like `rtracklayer::liftOver()` but "fills in" small gaps.
#'
#' This is identical to [rtracklayer::liftOver()], except that gaps
#' smaller than a specified size will be filled in. This allows the
#' liftover process to tolerate small indels in the old genome
#' relative to the new one without breaking up contiguous regions.
#'
#' @param x,chain,... These arguments have the same meaning as in
#'     [rtracklayer::liftOver()].
#' @param allow.gap Maximum gap size to "fill in". If set to zero,
#'     this is equivalent to [rtracklayer::liftOver()].
#' @return See [rtracklayer::liftOver()].
#'
#' @export
liftOverLax <- function(x, chain, ..., allow.gap = 0) {
    req_ns("rtracklayer", "S4Vectors", "GenomicRanges")
    newx <- rtracklayer::liftOver(x, chain, ...)
    if (allow.gap > 0) {
        gapped <- which(lengths(newx) > 1)
        newx.gapped.reduced <- GenomicRanges::reduce(
            newx[gapped], min.gapwidth = allow.gap + 1, with.revmap = TRUE)
        S4Vectors::mcols(newx.gapped.reduced@unlistData) <-
            rep(S4Vectors::mcols(x[gapped]), lengths(newx.gapped.reduced))
        newx[gapped] <- newx.gapped.reduced
    }
    return(newx)
}

#' Convenience function for running liftOver on a MotifMap BED file
#'
#' This function is a shortcut for calling [read_motifmap()], [liftOverLax()],
#' and then [write_motifmap()].
#'
#' @param infile The input file to read with `read_motifmap()`.
#' @param chaingile The chain file to use for the liftover process, to
#'     be read with [rtracklayer::import.chain()].
#' @param outfile The output file to write with `write_motifmap()`.
#' @param allow.gap,... These arguments are passed to `liftOverLax()`.
#'
#' Note that only features that remain contiguous after lifting over
#' are written to the output file. The rest are discarded.
#'
#' @export
liftOver_motifMap <- function(infile, chainfile, outfile, allow.gap = 2, ...) {
    req_ns("rtracklayer")
    gr <- read_motifmap(infile, parse_name = TRUE)
    chain <- rtracklayer::import.chain(chainfile)
    gr2 <- rtracklayer::liftOverLax(gr, chain, allow.gap = allow.gap, ...)
    gr2 <- unlist(gr2[lengths(gr2) == 1])
    write_motifmap(gr2, outfile)
}
