% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallel_utils.R
\name{featureCountsQuiet}
\alias{featureCountsQuiet}
\title{\code{Rsubread::featureCounts()} with output suppressed}
\usage{
featureCountsQuiet(...)
}
\arguments{
\item{...}{See \code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}.}
}
\value{
See \code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}.
}
\description{
This is a helper function for \code{featureCountsParallel()}.
}
\details{
This redirects the output to \verb{/dev/null}, so it assumes a UNIX-like
system.
}
