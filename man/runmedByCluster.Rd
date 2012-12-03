\name{runmedByCluster}
\alias{runmedByCluster}

\title{
Apply running median smoothing to values within each spatially-defined
cluster
}

\description{
Running median smoothing is applied independently to each cluster of
genomic locations. Locations within the same cluster are close together
to warrant smoothing across neighbouring locations. 
}

\usage{
runmedByCluster(y, x = NULL, cluster, weights = NULL, k = 5,
                endrule = "constant", verbose = TRUE) 
}

\arguments{
  \item{y}{A vector or matrix of values to be smoothed. If a matrix,
  each column represents a sample.}
  \item{x}{The genomic location of the values in y.}
  \item{cluster}{A vector indicating clusters of locations. A cluster is
    typically defined as a region that is small enough that it makes sense
    to smooth across neighbouring locations. Smoothing will only be applied
    within a cluster, not across locations from different clusters.} 
  \item{weights}{weights used by the smoother.}
  \item{k}{integer width of median window; must be odd. See \code{\link{runmed}}}
  \item{endrule}{character string indicating how the values at the
    beginning and the end (of the data) should be treated. See
    \code{\link{runmed}}.} 
  \item{verbose}{Boolean. Should progress be reported?}
}

\details{
This function is typically called by \code{\link{smoother}}, which is in
turn called by \code{\link{bumphunter}}.  
}

\value{
  \item{fitted }{The smoothed data values}
  \item{smoothed }{A boolean vector indicating whether a given position was smoothed}
  \item{spans }{The span used by the loess smoother. One per cluster.}
  \item{clusterL }{The number of locations in each cluster.}
  \item{smoother }{always set to \sQuote{runmed}.}
}

\author{
  Rafael A. Irizarry
}

\seealso{
  \code{\link{smoother}}, \code{\link{loessByCluster}}. Also see \code{\link{runmed}}.
}

\examples{
dat <- dummyData()
smoothed <- runmedByCluster(y=dat$mat[,1], cluster=dat$cluster,
                            k=5, endrule="constant")
}
