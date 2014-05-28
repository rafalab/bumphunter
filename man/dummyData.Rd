\name{dummyData}
\alias{dummyData}
\title{
  Generate dummy data for use with bumphunter functions
}

\description{
  This function generates a small dummy dataset representing
  samples from two different groups (cases and controls) that is used in
  bumphunter examples.
}

\usage{
dummyData(n1 = 5, n2 = 5, sd = 0.2, l = 100, spacing = 100,
          clusterSpacing=1e5, numClusters=5)
}

\arguments{
  \item{n1}{Number of samples in group 1 (controls)}
  \item{n2}{Number of samples in group 2 (cases)}
  \item{sd}{Within group standard deviation to be used when simulating data}
  \item{l}{The number of genomic locations for which to simulate data}
  \item{spacing}{The average spacing between locations. The actual
    locations have a random component so the actual spacing will be
    non-uniform}
  \item{clusterSpacing}{The spacing between clusters. (Specifically, the
    spacing between the first location in each cluster.)} 
  \item{numClusters}{Divide the genomic locations into this number of
    clusters, each of which will contain locations spaced \code{spacing}
    bp apart.} 
}

\value{
  A list containing data that can be used with various bumphunter functions.   
  \item{mat}{ A simulated data matrix with rows representing genomic locations and
  columns representing samples.}
  \item{design}{ Design matrix with rows representing samples and columns
  representing covariates.}  
  \item{chr}{A character vector with the chromosomes of each location.}
  \item{pos}{A numeric vector representing the chromosomal position.}
  \item{cluster}{A vector representing the cluster of each locations}
  \item{n1}{Number of samples in cluster 1}
  \item{n2}{Number of samples in cluster 2}
}

\author{
  Martin J. Aryee
}


\examples{
dat <- dummyData()
names(dat)
head(dat$pos)
}
