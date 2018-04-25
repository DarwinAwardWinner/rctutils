## readsToFragmentMidpoints
## strand.sign
## grl.to.saf
## gff.to.grl
## get.gene.common.colnames
## liftOverLax

liftOver_motifMap <- function(infile, chainfile, outfile, allow.gap=2) {
    gr <- read.motifmap(infile, parse_name=TRUE)
    chain <- import.chain(chainfile)
    gr2 <- liftOverLax(gr, chain, allow.gap=allow.gap)
    gr2 <- unlist(gr2[lengths(gr2) == 1])
    write.motifmap(gr2, outfile)
}

## get.tx2gene.from.txdb
