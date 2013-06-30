\name{locfitByCluster}
\alias{locfitByCluster}
\title{
  Apply local regression smoothing to values within each spatially-defined cluster.
}

\description{
  
Local regression smoothing with a gaussian kernal, is applied
independently to each cluster of genomic locations. Locations within the
same cluster are close together to warrant smoothing across neighbouring
locations.

}

\usage{
locfitByCluster(y, x = NULL, cluster, weights = NULL, minNum = 7,
               bpSpan = 1000, minInSpan = 0, verbose = TRUE)
}

\arguments{
  \item{y}{A vector or matrix of values to be smoothed. If a matrix,
    each column represents a sample.} 
  \item{x}{The genomic location of the values in y}
  \item{cluster}{A vector indicating clusters of locations. A cluster is
    typically defined as a region that is small enough that it makes sense
    to smooth across neighbouring locations. Smoothing will only be applied
    within a cluster, not across locations from different clusters.} 
  \item{weights}{weights used by the locfit smoother}
  \item{minNum}{Clusters with fewer than \code{minNum} locations will
    not be smoothed}
  \item{bpSpan}{The span used when locfit smoothing. (Expressed in base
    pairs.)}
  \item{minInSpan}{Only smooth the region if there are at least this
    many locations in the span.}
  \item{verbose}{Boolean. Should progress be reported?}
}
\details{
  This function is typically called by \code{\link{smoother}}, which is in
  turn called by \code{\link{bumphunter}}.
}
\value{
  \item{fitted }{The smoothed data values}
  \item{smoothed }{A boolean vector indicating whether a given position was smoothed}
  \item{smoother }{always set to \sQuote{locfit}.}
}
\author{
Rafael A. Irizarry and Kasper D. Hansen
}
\seealso{
\code{\link{smoother}}, \code{\link{runmedByCluster}}, \code{\link{loessByCluster}}
}
\examples{
dat <- dummyData()
smoothed <- locfitByCluster(y=dat$mat[,1], cluster=dat$cluster, bpSpan = 1000,
                         minNum=7, minInSpan=5)
}
