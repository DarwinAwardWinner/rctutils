% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_format_io.R
\name{read_saf}
\alias{read_saf}
\title{Read a SAF file into a GRangesList.}
\usage{
read_saf(filename, ...)
}
\arguments{
\item{filename}{The file name to read}

\item{...}{Additional arguments are passed to
\code{\link[=read_table_general]{read_table_general()}}.}
}
\value{
a GRangesList with one element for each GeneID, containing
the ranges for that gene.
}
\description{
The SAF format is described in \code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}.
}
