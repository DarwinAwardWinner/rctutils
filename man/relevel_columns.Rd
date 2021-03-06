% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_manip_utils.R
\name{relevel_columns}
\alias{relevel_columns}
\title{Relevel many factor columns in a data frame at once}
\usage{
relevel_columns(df, ...)
}
\arguments{
\item{df}{The data frame in which to re-level some factor columns.}

\item{...}{Additional arguments control the re-leveling of factors.
The name of each argument indicates the column of \code{df} to
re-level (this column should be a factor). The value should be
the new levels in the order that they should appear (or any
valid list of arguments to \code{\link[forcats:fct_relevel]{forcats::fct_relevel()}}).}
}
\value{
\code{df}, with the specified factor columns re-leveled
according to the provided specifications.
}
\description{
This is a shortcut for calling \code{\link[forcats:fct_relevel]{forcats::fct_relevel()}} on multiple
columns of a data frame.
}
\examples{

x <- data.frame(a = letters[1:2], b = LETTERS[1:6])
sapply(x, levels)
x2 <- relevel_columns(x, a = c("b", "a"), b = list("B", "A", after = 2))
sapply(x2, levels)

}
