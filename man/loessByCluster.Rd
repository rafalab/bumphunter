\name{loessByCluster}
\alias{loessByCluster}
\title{
  Apply loess smoothing to values within each spatially-defined cluster.
}

\description{
Loess smoothing is applied independently to each cluster of genomic
locations. Locations within the same cluster are close together to
warrant smoothing across neighbouring locations.
}

\usage{
loessByCluster(y, x = NULL, cluster, weights = NULL, bpSpan = 1000,
               minNum = 7, minInSpan = 5, maxSpan = 1, verbose = TRUE)
}

\arguments{
  \item{y}{A vector or matrix of values to be smoothed. If a matrix,
    each column represents a sample.} 
  \item{x}{The genomic location of the values in y}
  \item{cluster}{A vector indicating clusters of locations. A cluster is
    typically defined as a region that is small enough that it makes sense
    to smooth across neighbouring locations. Smoothing will only be applied
    within a cluster, not across locations from different clusters.} 
  \item{weights}{weights used by the loess smoother}
  \item{bpSpan}{The span used when loess smoothing. (Expressed in base
    pairs.)}
  \item{minNum}{Clusters with fewer than \code{minNum} locations will
    not be smoothed}
  \item{minInSpan}{Only smooth the region if there are at least this
    many locations in the span.}
  \item{maxSpan}{The maximum span. Spans greater than this value will be
    capped.} 
  \item{verbose}{Boolean. Should progress be reported?}
}
\details{
This function is typically called by \code{\link{smoother}}, which is in
  turn called by \code{\link{bumphunter}}.
}
\value{
  \item{fitted }{The smoothed data values}
  \item{smoothed }{A boolean vector indicating whether a given position was smoothed}
  \item{smoother }{always set to \sQuote{loess}.}
}
\author{
Rafael A. Irizarry
}
\seealso{
\code{\link{smoother}}, \code{\link{runmedByCluster}}, \code{\link{locfitByCluster}}
}
\examples{
dat <- dummyData()
smoothed <- loessByCluster(y=dat$mat[,1], cluster=dat$cluster, bpSpan = 1000,
                         minNum=7, minInSpan=5, maxSpan=1)
}
