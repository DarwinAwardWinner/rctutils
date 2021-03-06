% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_format_io.R
\name{write_motifmap}
\alias{write_motifmap}
\title{Write a MotifMap-type GRanges into a BED file.}
\usage{
write_motifmap(x, file)
}
\arguments{
\item{x}{A GRanges object with appropriate metadata columns.}

\item{file}{The file name to save.}
}
\description{
This takes a GRanges object with the appropriate mcols (see
\code{\link[=read_motifmap]{read_motifmap()}}) and saves it as a BED file.
}
