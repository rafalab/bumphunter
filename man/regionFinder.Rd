\name{regionFinder}
\alias{regionFinder}
\title{
  Find non-zero regions in vector
}

\description{
  Find regions for which a numeric vector is above (or below) predefined thresholds.
}

\usage{
regionFinder(x, chr, pos, cluster = NULL, y = x, summary = mean,
             ind = seq(along = x), order = TRUE, oneTable = TRUE,
             maxGap = 300, cutoff=quantile(abs(x), 0.99),
             assumeSorted = FALSE, verbose = TRUE)
}

\arguments{
  \item{x}{A numeric vector.}
  \item{chr}{A character vector with the chromosomes of each location.}
  \item{pos}{A numeric vector representing the genomic location.}
  \item{cluster}{The clusters of locations that are to be analyzed
    together. In the case of microarrays, the cluster is many times
    supplied by the manufacturer. If not available the function
    \code{\link{clusterMaker}} can be used.}
  \item{y}{
    A numeric vector with same length as \code{x} containing values to
    be averaged for the region summary. See details for more.}
  \item{summary}{
    The function to be used to construct a summary of the \code{y}
    values for each region.}
  \item{ind}{
    an optional vector specifying a subset of observations to be used
    when finding regions.}
  \item{order}{
    if \code{TRUE} then the resulting tables are ordered
    based on area of each region. Area is defined as the absolute value
    of the summarized \code{y} times the number of features in the
    regions.}
  \item{oneTable}{
    if \code{TRUE} only one results table is returned. Otherwise, two
    tables are returned: one for the regions with positive values and
    one for the negative values.}
  \item{maxGap}{
    If cluster is not provided this number will be used to define clusters
    via the \code{\link{clusterMaker}} function.}
  \item{cutoff}{
    This argument is passed to \code{\link{getSegments}}.  It represents
    the upper (and optionally the lower) cutoff for \code{x}. }
  \item{assumeSorted}{
    This argument is passed to \code{\link{getSegments}} and
    \code{\link{clusterMaker}}.} 
  \item{verbose}{
    Should the function be verbose?}
}

\details{

  This function is used in the final steps of
  \code{\link{bumphunter}}. While \code{bumphunter} does many things,
  such as regression and permutation, \code{regionFinder} simply finds
  regions that are above a certain threshold (using
  \code{\link{getSegments}}) and summarizes them. The regions are found
  based on \code{x} and the summarized values are based on \code{y}
  (which by default equals \code{x}). The summary is used for the
  ranking so one might, for example, use t-tests to find regions but
  summarize using effect sizes.

}

\value{
  
 If \code{oneTable} is \code{FALSE} it returns two tables otherwise it
 returns one table. The rows of the table are regions. Information on
 the regions is included in the columns.

}

\author{
  Rafael A Irizarry
}

\seealso{

  \code{\link{bumphunter}} for the main usage of this function,
  \code{\link{clusterMaker}} for the typical input to the \code{cluster}
  argument and \code{\link{getSegments}} for a function used within
  \code{regionFinder}.

}

\examples{
x <- seq(1:1000)
y <- sin(8*pi*x/1000) + rnorm(1000, 0, 0.2)
chr <- rep(c(1,2), each=length(x)/2)
tab <- regionFinder(y, chr, x, cutoff=0.8)
print(tab[tab$L>10,])
}
