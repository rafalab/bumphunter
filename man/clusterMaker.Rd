\name{clusterMaker}
\alias{clusterMaker}
\title{
  Make clusters of genomic locations based on distance
}

\description{
  Genomic locations are grouped into clusters based on distance:
  locations that are close to each other are assigned to the same
  cluster. The operation is performed on each chromosome independently.
}

\usage{
clusterMaker(chr, pos, assumeSorted = FALSE, maxGap = 300)
}

\arguments{
  \item{chr}{A vector representing chromosomes.  This is usually a
    character vector, but may be a factor or an integer.}
  \item{pos}{A numeric vector with genomic locations.}
  \item{assumeSorted}{This is a statement that the function may assume
    that the vector \code{pos} is sorted (within each \code{chr}).
    Allowing the function to make this assumption may increase the speed
    of the function slightly.} 
  \item{maxGap}{An integer. Genomic locations within maxGap from each other are
    placed into the same cluster.}
}

\details{
  
  The main purpose of the function is to genomic location into clusters
  that are close enough to perform operations such as smoothing.  A
  genomic location is a combination of a chromosome (\code{chr}) and an
  integer position (\code{pos}).  Specifically, genomic intervals are
  not handled by this function.

  Each chromosome is clustered independently from each other.  Within
  each chromosome, clusters are formed in such a way that two positions
  belong to the same cluster if they are within \code{maxGap} of each
  other.
}

\value{
  A vector of integers to be interpreted as IDs for the clusters, such
  that two genomic positions with the same cluster ID is in the same
  cluster.  Each genomic position receives one integer ID.
}

\author{
  Rafael A. Irizarry 
}

\examples{
N <- 1000
chr <- sample(1:5, N, replace=TRUE)
pos <- round(runif(N, 1, 10^5))
o <- order(chr, pos)
chr <- chr[o]
pos <- pos[o]
regionID <- clusterMaker(chr, pos)
}
