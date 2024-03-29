% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/granges_utils.R
\name{liftOver_motifMap}
\alias{liftOver_motifMap}
\title{Convenience function for running liftOver on a MotifMap BED file}
\usage{
liftOver_motifMap(infile, chainfile, outfile, allow.gap = 2, ...)
}
\arguments{
\item{infile}{The input file to read with \code{read_motifmap()}.}

\item{chainfile}{The chain file to use for the liftover process, to
be read with \code{\link[rtracklayer:Chain-class]{rtracklayer::import.chain()}}.}

\item{outfile}{The output file to write with \code{write_motifmap()}.}

\item{allow.gap, ...}{These arguments are passed to \code{liftOverLax()}.

Note that only features that remain contiguous after lifting over
are written to the output file. The rest are discarded.}
}
\description{
This function is a shortcut for calling \code{\link[=read_motifmap]{read_motifmap()}}, \code{\link[=liftOverLax]{liftOverLax()}},
and then \code{\link[=write_motifmap]{write_motifmap()}}.
}
