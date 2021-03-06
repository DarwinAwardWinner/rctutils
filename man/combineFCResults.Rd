% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallel_utils.R
\name{combineFCResults}
\alias{combineFCResults}
\title{Combine multiple featureCounts results}
\usage{
combineFCResults(fcreslist)
}
\arguments{
\item{fcreslist}{A list of values returned by multiple calls to
\code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}}. They should all contain the same
set of features but different samples.}
}
\value{
A single value of the same type that is returned by
\code{Rsubread::featureCounts()}. It will contain the same features
as the inputs and the union of all samples.
}
\description{
This is a helper function for \code{featureCountsParallel()}.
}
