\name{getSegments}
\alias{getSegments}
\title{
  Segment a vector into positive, zero, and negative regions
}
\description{
  Given two cutoffs, L and U, this function divides a numerical vector
  into contiguous parts that are above U, between L and U, and below L.
}
\usage{
getSegments(x, f = NULL, cutoff = quantile(abs(x), 0.99),
            assumeSorted = FALSE, verbose = FALSE)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{f}{A factor used to pre-divide \code{x} into pieces. Each
    piece is then segmented based on the cutoff.  Setting this to
    \code{NULL} assumes that there are no pre-division.  Often,
    \code{\link{clusterMaker}} is used to define this factor.}
  \item{cutoff}{a numeric vector of length either 1 or 2.  If length is
    1, U (see details) will be cutoff and L will be -cutoff.  Otherwise
    it specifies L and U.  The function will furthermore always use the
    minimum of cutoff for L and the maximum for U.}
  \item{assumeSorted}{This is a statement that the function may assume
    that the vector \code{f} is sorted.  Allowing the function to make
    this assumption may increase the speed of the function slightly.}
  \item{verbose}{Should the function be verbose?}
}
\details{
  This function is used to find the indexes of the \sQuote{bumps} in functions
  such as \code{\link{bumphunter}}. 

  \code{x} is a numeric vector, which is converted into three levels
  depending on whether x>=U (\sQuote{up}), L<x<U (\sQuote{zero}) or x<=L
  (\sQuote{down}), with L and U coming from \code{cutoff}.  We assume
  that adjacent entries in \code{x} are next to each other in some
  sense.  Segments, consisting of consecutive indices into \code{x}
  (ie. values between 1 and \code{length(x)}), are formed such that all
  indices in the same segment belong to the same level of \code{f} and
  have the same discretized value of \code{x}.

  In other words, we can use \code{getSegments} to find runs of \code{x}
  belonging to the same level of \code{f} and with all of the values of
  x either above U, between L and U, or below L.
}

\value{
  A list with three components, each a list of indices.  Each component
  of these lists represents a segment and this segment is represented by a
  vector of indices into the original vectors \code{x} and \code{f}.
  \item{upIndex:}{a list with each entry an index of contiguous values
    in the same segment.  These segments have values of \code{x} above U.}
  \item{dnIndex:}{a list with each entry an index of contiguous values
    in the same segment.  These segments have values of \code{x} below L.}
  \item{zeroIndex:}{a list with each entry an index of contiguous
    values in the same segment. These segments have values of \code{x}
    between L and U.}
}

\author{
  Rafael A Irizarry and Kasper Daniel Hansen
}

\seealso{
  \code{\link{clusterMaker}}
}
\examples{
x <- 1:100
y <- sin(8*pi*x/100)
chr <- rep(1, length(x))
indexes <- getSegments(y, chr, cutoff=0.8)
plot(x, y, type="n")
for(i in 1:3){
   ind <- indexes[[i]]
   for(j in seq(along=ind)) {
      k <- ind[[j]]
      text(x[k], y[k], j, col=i)
   }
}
abline(h=c(-0.8,0.8))
}
