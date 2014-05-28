\name{pointMatch}
\alias{pointMatch}
\title{
Fuzzy Matching of Genomic Locations
}

\description{
  Finds the nearest locations in a given table of genomic targets
}

\usage{
pointMatch(object1, object2, col1=2, col2=2, verbose)
}

\arguments{
  \item{object1}{A \code{data.frame} holding the query positions in column
\code{col1}, on the chromosomes given in the \code{chr} column.}
  \item{object2}{As \code{object1}, but holding the target positions in
the \code{col2} and \code{chr} columns.}
  \item{col1}{The column in \code{object1} holding the query positions.}
  \item{col2}{The column in \code{object2} holding the target positions.}
  \item{verbose}{For backward compatibility; ignored.}
}

\details{
Similar to \code{annotateNearest}, but for points instead of ranges.
}

\value{
An integer matrix with columns \code{c("dist", "matchIndex")}.
\code{dist} is the signed distance from a query to its nearest target,
which is given by its index \code{matchIndex} in \code{object2}.  The
distance is negative for queries past their nearest target.  A row of
\code{NA} is returned for queries on chromosomes not in \code{object2}.
}

\author{
Harris Jaffee, Peter Murakami and Rafael A. Irizarry
}
\seealso{
\code{annotateNearest}, \code{nearest}, and \code{distanceToNearest}
}

\examples{
object1 <- data.frame(3, chr="chr1")
object2 <- data.frame(c(1,4), chr="chr1")
pointMatch(object1, object2, col1=1, col2=1)
object2 <- data.frame(c(2,4), chr="chr1")
pointMatch(object1, object2, col1=1, col2=1)
}
