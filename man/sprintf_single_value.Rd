% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_manip_utils.R
\name{sprintf_single_value}
\alias{sprintf_single_value}
\title{Like \code{sprintf()} but using the same value in each specification}
\usage{
sprintf_single_value(fmt, value)
}
\arguments{
\item{fmt}{This has the same meaning as in \code{sprintf()}.}

\item{value}{The single value use for every conversion
specification in \code{fmt}.}
}
\value{
See \code{\link[=sprintf]{sprintf()}}.

Note that vector arguments are still recycled normally.
}
\description{
This is equivalent to calling \code{\link[=sprintf]{sprintf()}} with \code{value} passed after
\code{fmt} as many times as there are conversion specifications in
\code{fmt}.
}
\examples{

fmt <- "Hello \%s. Your name is \%s. Isn't it good to be named \%s?"
# With normal sprintf, you need to specify the right number of
# arguments:
name <- c("John Doe", "Jane Doe")
sprintf(fmt, name, name, name)
# Instead, this figures out the right number of times to pass it:
sprintf_single_value(fmt, name)

}
