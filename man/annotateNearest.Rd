\name{annotateNearest}
\alias{annotateNearest}
\alias{regionMatch}
\alias{matchGenes}
\title{
Annotate Nearest
}

\description{
  Finds the nearest region in a table of target regions and annotates
  the relationship.
}

\usage{
annotateNearest(x, subject, annotate=TRUE, ...)

matchGenes(x, build, ...)
}

\arguments{
  \item{x}{An \code{IRanges} or \code{GenomicRanges} object, or a
    \code{data.frame} with columns for \code{start}, \code{end},
    and, optionally, \code{chr} or \code{seqnames}.}
  \item{subject}{For \code{annotateNearest} to operate in
    \code{matchGenes} mode, a single word indicating the desired
    species/build; this is currently limited to \code{hg19}.
    Otherwise, an object of the same class as \code{x}.}
  \item{annotate}{Whether to annotate the result.}
  \item{build}{For example, \code{hg19}.}
  \item{\dots}{Options passed to \code{nearest}.}
}

\details{
The name \code{matchGenes} is provided for backward compatibility.
\code{annotateNearest} is a wrapper for \code{nearest(select="arbitrary")}
with annotation of the result.
}

\value{
  
For \code{matchGenes}, a data frame with columns \code{c("name",
  "annotation", "description", "region", "distance", "subregion",
  "insidedistance", "exonnumber", "nexons", "UTR", "strand", "geneL",
  "codingL")} unless \code{all} is \code{TRUE}, in which case there
  are two additional columns at the beginning holding the query index
  and nearest transcript name.
\item{name}{}
\item{annotation}{}
\item{description}{}
\item{region}{}
\item{distance}{}
\item{subregion}{}
\item{insidedistance}{}
\item{exonnumber}{}
\item{nexons}{}
\item{UTR}{}
\item{strand}{}
\item{geneL}{}
\item{codingL}{}

Otherwise, a data frame with columns \code{c("dist", "matchIndex", "type",
  "amountOverlap", "insideDist", "size1", "size2")} unless
  \code{annotate} is \code{FALSE}, in which case the first two columns,
  only, are returned as an integer matrix. 
\item{dist}{Signed distance to the nearest target. Queries downstream
  from (i.e. past) their nearest target are given a negative distance.}
\item{matchIndex}{The index of the nearest target.}
\item{type}{one of \code{c("inside", "cover", "disjoint", "overlap")}.}
\item{amountOverlap}{The width of the overlap region, if any.}
\code{insideDist}{When a query is contained in its nearest target, the
  signed minimum of the distances "target start to query start" and
  "query end to target end". The former is taken positive, and the latter,
  which wins in ties, negative. \code{dist} will be \code{0} in this case.}
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
#
## showing 'cover' and 'disjoint', and 'amountOverlap'
annotateNearest(query, subject)
#
## showing 'inside' and 'insideDist', and 'amountOverlap'
annotateNearest(subject, query)
annotateNearest(IRanges(3,3), IRanges(2,5))
annotateNearest(IRanges(3,4), IRanges(2,5))
annotateNearest(IRanges(4,4), IRanges(2,5))
}
