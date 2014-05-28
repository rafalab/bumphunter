\name{smoother}
\alias{smoother}
\title{
Smooth genomic profiles
}

\description{
Apply smoothing to values typically representing the difference between
two populations across genomic regions.
}

\usage{
smoother(y, x = NULL, cluster, weights = NULL, smoothFunction,
         verbose = TRUE, ...) 
}

\arguments{
  \item{y}{A vector or matrix of values to be smoothed. If a matrix,
  each column represents a sample.}
  \item{x}{The genomic location of the values in y}
  \item{cluster}{A vector indicating clusters of locations. A cluster
    is typically defined as a region that is small enough that it makes
    sense to smooth across neighbouring  locations. Smoothing will only be
    applied within a cluster, not across locations from different
    clusters} 
  \item{weights}{weights used by the smoother.}
  \item{smoothFunction}{A function to be used for smoothing the estimate
    of the genomic profile. Two functions are provided by the package:
    \code{loessByCluster} and \code{runmedByCluster}.}
  \item{verbose}{Boolean. Should progress be reported?}
  \item{\dots}{Further arguments to be passed to \code{smoothFunction}}
}

\details{
This function is typically called by bumphunter prior to identifying
candidate bump regions. Smoothing is carried out within regions defined
by the \code{cluster} argument.
}

\value{
  \item{fitted }{The smoothed data values}
  \item{smoothed }{A boolean vector indicating whether a given position was smoothed}
  \item{spans }{The span used by the loess smoother. One per cluster.}
  \item{clusterL }{The number of locations in each cluster.}
  \item{smoother }{The name of the smoother used}
}

\author{
Rafael A. Irizarry and Martin J. Aryee
}

\seealso{
\code{\link{loessByCluster}}, \code{\link{runmedByCluster}}
}

\examples{
dat <- dummyData()

# Enable parallelization
require(doParallel)
registerDoParallel(cores = 2)

## loessByCluster
smoothed <- smoother(y=dat$mat[,1], cluster=dat$cluster, smoothFunction=loessByCluster, 
                       bpSpan = 1000, minNum=7, minInSpan=5, maxSpan=1)
## runmedByCluster
smoothed <- smoother(y=dat$mat[,1], cluster=dat$cluster, smoothFunction=runmedByCluster, 
                     k=5, endrule="constant")
# cleanup, for Windows
bumphunter:::foreachCleanup()
}
