\name{annotateNearest}
\alias{annotateNearest}
\title{
Annotate the results of nearest
}
\description{
Annotate the results of nearest with more information about the type of match.
}
\usage{
annotateNearest(x, subject, annotate = TRUE, ...)
}
\arguments{
  \item{x}{
The query. An \code{IRanges} or \code{GenomicRanges} object, or a
    \code{data.frame} with columns for \code{start}, \code{end},
    and, optionally, \code{chr} or \code{seqnames}.
}
  \item{subject}{
The subject.  An \code{IRanges} or \code{GenomicRanges} object, or a
    \code{data.frame} with columns for \code{start}, \code{end},
    and, optionally, \code{chr} or \code{seqnames}.
}
  \item{annotate}{
Whether to annotate the result.
}
  \item{\dots}{
Arguments passed along to \code{nearest}. 
}
}
\details{
This function runs \code{\link[IRanges:nearest-methods]{IRanges::nearest()}} and then annotates the
nearest hit. Note that the nearest subject range to a given query may not be
unique and we arbitrarily chose one as done by default by
\code{\link[IRanges:nearest-methods]{IRanges::nearest()}}.
}
\value{
A data frame with columns \code{c("distance", "subjectHits", "type",
  "amountOverlap", "insideDistamce", "size1", "size2")} unless
  \code{annotate} is \code{FALSE}, in which case only the first two
  columns are returned as an integer matrix.

 \item{dist}{Signed distance to the nearest target. Queries downstream
  from (i.e. past) their nearest target are given a negative distance.}
\item{subjectHits}{The index of the nearest target.}
\item{type}{one of \code{c("inside", "cover", "disjoint", "overlap")}.}
\item{amountOverlap}{The width of the overlap region, if any.}
\item{insideDistance}{When a query is contained in its nearest target, the
  signed minimum of the two distances target-start-to-query-start and
  query-end-to-target-end. The former is taken positive, and the latter,
  which wins in ties, negative. \code{dist} will be \code{0} in this case.}
\item{size1}{equals \code{width(x)}.}
\item{size2}{equals \code{width(subject)}.}
}
\author{
Harris Jaffee, Peter Murakami and Rafael A. Irizarry
}

\seealso{
\code{\link[IRanges:nearest-methods]{IRanges::nearest()}}, \code{\link{matchGenes}} 
}

\examples{
query <- GRanges(seqnames = 'chr1', IRanges(c(1, 4, 9), c(5, 7, 10)))
subject <- GRanges('chr1', IRanges(c(2, 2, 10), c(2, 3, 12)))
nearest(query, subject)
distanceToNearest(query, subject)

## showing 'cover' and 'disjoint', and 'amountOverlap'
annotateNearest(query, subject)

## showing 'inside' and 'insideDist', and 'amountOverlap'
annotateNearest(subject, query)
annotateNearest(GRanges('chr1', IRanges(3,3)), GRanges('chr1', IRanges(2,5)))
annotateNearest(GRanges('chr1', IRanges(3,4)), GRanges('chr1', IRanges(2,5)))
annotateNearest(GRanges('chr1', IRanges(4,4)), GRanges('chr1', IRanges(2,5)))

}
