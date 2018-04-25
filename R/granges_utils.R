## Should really be an S4 method, but writing S4 methods is a pain
readsToFragmentMidpoints <- function(reads, fraglen) {
    if (is(reads, "GAlignmentsList")) {
        isSingle <- lengths(reads) == 1
        mated.frags <- granges(reads[!isSingle], ignore.strand=TRUE)
        ## Extend smaller reads to fraglen, but don't shrink longer
        ## reads
        single.frags <- granges(reads[isSingle]) %>% resize(width=pmax(fraglen, width(.)), fix="start")
        frags <- c(mated.frags, single.frags)
    } else if (is(reads, "GAlignmentPairs")) {
        frags <- granges(reads)
    } else if (is(reads, "GAlignments")) {
        ## Extend smaller reads to fraglen, but don't shrink longer
        ## reads
        frags <- reads %>% granges %>% resize(width=pmax(fraglen, width(.)), fix="start")
    } else {
        warning(glue("Unknown reads type: {class(reads)}. Attempting to coerce to GRanges."))
        reads %<>% as("GRanges")
        ## If all ranges are the same width, assume they represent
        ## single-end reads and resize them to fraglen
        if (all(width(reads) == width(reads[1])) && width(reads[1]) < fraglen) {
            warning(glue("All reads from unknown class have the same length, {width(reads[1])}, and are therefore assumed to be single-end reads, which will be resized to {fraglen}."))
            frags <- resize(reads, width=fraglen, fix="start")
        } else {
            frags <- reads
        }
    }
    ## Finally, get the center of each fragment
    resize(frags, width=1, fix="center")
}

## Convert strand to -1, 0, or 1
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

## Given a GRangesList whose underlying ranges have mcols, find mcols
## of the ranges that are constant within each gene and promote them
## to mcols of the GRangesList. For example, if exons are annotated with
promote.common.mcols <- function(grl, delete.from.source=FALSE, ...) {
    colnames.to.promote <- get.gene.common.colnames(mcols(unlist(grl)), rep(names(grl), lengths(grl)), ...)
    promoted.df <- mcols(unlist(grl))[cumsum(lengths(grl)),colnames.to.promote, drop=FALSE]
    if (delete.from.source) {
        mcols(grl@unlistData) %<>% .[setdiff(names(.), colnames.to.promote)]
    }
    mcols(grl) %<>% cbind(promoted.df)
    grl
}

## Like rtracklayer::liftOver but "fills in" small gaps induced by the
## liftOver process (i.e. no larger than allow.gap). If allow.gap is
## zero, this is equivalent to liftOver.
#' @importFrom rtracklayer liftOver
liftOverLax <- function(x, chain, ..., allow.gap=0) {
    newx <- liftOver(x, chain)
    if (allow.gap > 0) {
        gapped <- which(lengths(newx) > 1)
        newx.gapped.reduced <- reduce(newx[gapped], min.gapwidth = allow.gap + 1, with.revmap=TRUE)
        mcols(newx.gapped.reduced@unlistData) <- rep(mcols(x[gapped]), lengths(newx.gapped.reduced))
        newx[gapped] <- newx.gapped.reduced
    }
    return(newx)
}

liftOver_motifMap <- function(infile, chainfile, outfile, allow.gap=2) {
    gr <- read.motifmap(infile, parse_name=TRUE)
    chain <- import.chain(chainfile)
    gr2 <- liftOverLax(gr, chain, allow.gap=allow.gap)
    gr2 <- unlist(gr2[lengths(gr2) == 1])
    write.motifmap(gr2, outfile)
}
