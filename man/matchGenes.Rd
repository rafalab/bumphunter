\name{matchGenes}
\alias{matchGenes}
\title{
Human Gene Matcher
}

\description{
  Finds the nearest transcript(s) and annotates the relationship.
}

\usage{
matchGenes(x, subject, annotate=TRUE, ...)
}

\arguments{
  \item{object}{}
}

\details{
}

\value{
A data frame with columns \code{c("name", "annotation", "description",
  "levels", "region", "distance", "subregion", "insidedistance",
  "exonnumber", "nexons", "UTR", "strand", "geneL", "codingL")} unless
  \code{all} is \code{TRUE}, in which case there are two additional
  columns at the beginning holding the query indices and the nearest
  transcript names.
}

\author{
Peter Murakami and Rafael A. Irizarry
}
\seealso{
\link{annotateNearest}
}

\examples{
}
