\name{annotateNearest}
\alias{annotateNearest}
\title{
Annotate Nearest
}

\description{
  Finds the nearest region in a table of target regions and annotates
  the relationship.
}

\usage{
annotateNearest(x, subject, pointMatch=FALSE, ...)
}

\arguments{
  \item{x}{An \code{IRanges} or \code{GenomicRanges} object, or a
    \code{data.frame} with columns for \code{start}, \code{end},
    and, optionally, \code{chr} or \code{seqnames}.}
  \item{subject}{Object of the same class as \code{x}.}
  \item{pointMatch}{FIXME}
  \item{\dots}{Options passed to \code{nearest}.}
}

\details{
A wrapper for \code{nearest(select="arbitrary")} with annotation of
the result.
}

\value{
  
A data frame with columns \code{c("dist", "matchIndex", "type",
  "amountOverlap", "insideDist", "size1", "size2")} unless
  \code{annotate} is \code{FALSE}, in which case only the first two
  columns are returned, as an integer matrix. 
\item{item}{The signed distance to the nearest target, with queries
downstream from (past) their nearest target given a negative distance.}
\item{matchIndex}{The index of the nearest target.}
\item{type}{one of \code{c("inside", "cover", "disjoint", "overlap")}.}
\item{amountOverlap}{The width of the overlap region (when there is one).}
\code{insideDist}{The distance when one region is inside another. dist
  will be 0 in this case.}
\item{size1}{equals \code{width(x)}.}
\item{size2}{equals \code{width(subject)}.}

}

\author{
Harris Jaffee, Peter Murakami and Rafael A. Irizarry
}
\seealso{
\code{nearest} and \code{distanceToNearest}
}

\examples{
query <- IRanges(c(1, 4, 9), c(5, 7, 10))
subject <- IRanges(c(2, 2, 10), c(2, 3, 12))
nearest(query, subject)
distanceToNearest(query, subject)
annotateNearest(query, subject)
}
