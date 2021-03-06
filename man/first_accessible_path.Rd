% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_utils.R
\name{first_accessible_path}
\alias{first_accessible_path}
\title{Find the first accessible file path from a vector of paths.}
\usage{
first_accessible_path(paths, mode = 0)
}
\arguments{
\item{paths}{"A character vector of paths to check for accessibility, in order."}

\item{mode}{"Has the same meaning as in \code{file.access()}."}
}
\value{
The first element of \code{paths} for which
\code{file.access()} returns TRUE. If none of \code{paths} is
accessible, NA is returned.
}
\description{
This function is useful for finding a file that may be in one of
several different locations.
}
\examples{
exec_paths <- strsplit(Sys.getenv("PATH"), ":")[[1]]
first_accessible_path(exec_paths)
}
