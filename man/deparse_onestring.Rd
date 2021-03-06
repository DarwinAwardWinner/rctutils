% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixes_for_builtins.R
\name{deparse_onestring}
\alias{deparse_onestring}
\title{Deparse and then concatenate into a single string}
\usage{
deparse_onestring(...)
}
\arguments{
\item{...}{See \code{\link[=deparse]{deparse()}}.}
}
\value{
A character vector of length 1.

Note that depending on the different deparse options, the resulting
string may still not be a complete R expression that recreates the
original object. Also note that the returned string may contain
newlines.
}
\description{
This is like the base R \code{deparse()}, except that it always returns
a single string, rather than a character vector of lines of text.
}
\examples{

# Make an object whose representation spans multiple lines
x <- seq(100)[-50]
length(deparse(x))
length(deparse_onestring(x))
identical(x, eval(parse(text=deparse_onestring(x))))

}
