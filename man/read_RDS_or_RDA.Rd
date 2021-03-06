% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdata_io.R
\name{read_RDS_or_RDA}
\alias{read_RDS_or_RDA}
\title{Read a single object from an RDS or RDA file}
\usage{
read_RDS_or_RDA(filename, expected.class = "ANY")
}
\arguments{
\item{filename}{The file to read an R object from.}

\item{expected.class}{If specified, the object will be coerced into
this class using \code{as()}, which will throw an error as normal if
the coercion is not possible. This allows you to restrict what
kind of objects you will accept. This can also be a function
that accepts a single argument and performs the proper coercion
itself.}
}
\value{
The object read from the file, possibly after coercing it
into another class.
}
\description{
This is like \code{\link[=readRDS]{readRDS()}} except that it can also read from an RData
file containing a single object (i.e. the kind of file that is read
using \code{\link[=load]{load()}}). Use this in place of \code{readRDS()} if you want to be
slightly more forgiving about what kind of R data file you accept.
}
\examples{

tmpf <- tempfile()
saveRDS(1:10, tmpf)
read_RDS_or_RDA(tmpf)
read_RDS_or_RDA(tmpf, "character")
# Using a function instead of a class name.
read_RDS_or_RDA(tmpf, as.character)
read_RDS_or_RDA(tmpf, "factor")
read_RDS_or_RDA(tmpf, "data.frame")

\dontrun{
# This will throw an error because the coercion to "lm" is not
# possible.
read_RDS_or_RDA(tmpf, "lm")
}

}
\seealso{
\code{\link[=readRDS]{readRDS()}}, \code{\link[=read_single_object_from_rda]{read_single_object_from_rda()}}, \code{\link[=as]{as()}}
}
