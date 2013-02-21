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
  \item{\dots}{For \code{matchGenes} or \code{matchGenes} mode,
	\code{promoterDist} (default \code{2500}), \code{verbose}
	(\code{TRUE}), \code{all} (\code{FALSE}), and \code{mc.cores}
	(\code{1}).  \code{mc.cores} is crucial for parallelizing
	the annotation process, assuming the hardware facilities.
	But note that the default is NOT \code{getOption("mc.cores", 2L)},
	yet, which leaves the explicit decision to the user.
	\code{up} and \code{down} limits on the nearest gene location
	are planned but not implemented yet.  Otherwise, options for
	\code{nearest}, except that \code{select="all"} is not
	implemented yet and is currently rejected.}
}

\details{
The name \code{matchGenes} is provided for backward compatibility.
The preferred call is to \code{nearest} with \code{subject} set to
the build name.
\code{annotateNearest} and \code{matchGenes} are basically wrappers
for \code{nearest(select="arbitrary")} with specialized annotation
of the result.
}

\value{

Normally, a data frame with columns \code{c("dist", "matchIndex", "type",
  "amountOverlap", "insideDist", "size1", "size2")} unless
  \code{annotate} is \code{FALSE}, in which case the first two columns,
  only, are returned as an integer matrix.

For matchGenes or matchGenes mode, with \code{all} unset or \code{FALSE},
  a data frame with one row for each query and with columns \code{c("name",
  "annotation", "description", "region", "distance", "subregion",
  "insidedistance", "exonnumber", "nexons", "UTR", "strand", "geneL", "codingL")}.
  The first column is the _gene_ nearest the query, by virtue of owning the
  transcript nearest the query.  Note that this nearest transcript (name),
  per se, does not exist in the result.  If \code{all} is \code{TRUE}, then each
  row of the previous result is expanded into perhaps multiple rows, reflecting
  multiple (equally) nearest transcripts for a given query.  In this case, two
  new columns at the beginning hold the query index and nearest transcript name.
  The nearest gene to a given query (in column 3) may or may not be unique.
  
\item{dist}{Signed distance to the nearest target. Queries downstream
  from (i.e. past) their nearest target are given a negative distance.}
\item{matchIndex}{The index of the nearest target.}
\item{type}{one of \code{c("inside", "cover", "disjoint", "overlap")}.}
\item{amountOverlap}{The width of the overlap region, if any.}
\item{insideDist}{When a query is contained in its nearest target, the
  signed minimum of the two distances target-start-to-query-start and
  query-end-to-target-end. The former is taken positive, and the latter,
  which wins in ties, negative. \code{dist} will be \code{0} in this case.}
\item{size1}{equals \code{width(x)}.}
\item{size2}{equals \code{width(subject)}.}

\item{name}{nearest gene}
\item{annotation}{RefSeq ID}
\item{description}{a factor with levels \code{c("upstream", "promoter",
	"overlaps 5'", "inside intron", "inside exon", "covers exon(s)",
	"overlaps exon upstream", "overlaps exon downstream",
	"overlaps two exons", "overlaps 3'", "close to 3'", "downstream",
	"covers")}}
\item{region}{a factor with levels \code{c("upstream", "promoter", "overlaps 5'",
	"inside", "overlaps 3'", "close to 3'", "downstream", "covers")}}
\item{distance}{distance before 5' end of gene}
\item{subregion}{a factor with levels \code{c("inside intron", "inside exon",
	"covers exon(s)", "overlaps exon upstream", "overlaps exon downstream",
	"overlaps two exons")}}
\item{insidedistance}{distance past 5' end of gene}
\item{exonnumber}{which exon}
\item{nexons}{number of exons}
\item{UTR}{a factor with levels \code{c("inside transcription region", "5' UTR",
	"overlaps 5' UTR", "3'UTR", "overlaps 3'UTR", "covers transcription region")}}
\item{strand}{"+" or "-"}
\item{geneL}{the gene length}
\item{codingL}{the coding length}

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
