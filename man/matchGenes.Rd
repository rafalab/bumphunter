\name{matchGenes}
\alias{matchGenes}
\title{
Find and annotate closest genes to genomic regions
}
\description{
Find and annotate closest genes to genomic regions
}
\usage{
matchGenes(x, subject, type = c("any", "fiveprime"), promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{
    An \code{IRanges} or \code{GenomicRanges} object, or a
    \code{data.frame} with columns for \code{start}, \code{end},
    and, optionally, \code{chr} or \code{seqnames}.
  }
  \item{subject}{
    An \code{GenomicRanges} object containing transcripts or genes that
    have been annotated by the function
    \code{\link{annotateTranscripts}}. 
  }
  \item{promoterDist}{
    Anything within this distance to the transcription start site (TSE)
    will be considered a promoter.
  }
  \item{type}{
    Should the distance be computed to any part of the transcript or the
    five prime end.
}
\item{skipExons}{
    Should the annotation of exons be skipped. Skipping this part makes the
    code slightly faster.
  }
  \item{verbose}{
    logical value. If 'TRUE', it writes out some messages
    indicating progress. If 'FALSE' nothing should be printed.
  }
}
\details{
 This function runs \code{\link[IRanges:nearest-methods]{IRanges::nearest()}} and then annotates
 the the relationship between the region and the transcript/gene that is
 closest. Many details are provided on this relationship as described in
 the next section. 

}
\value{
  A data frame with one row for each query and with columns \code{c("name",
  "annotation", "description", "region", "distance", "subregion",
  "insideDistance", "exonnumber", "nexons", "UTR", "strand", "geneL", "codingL","Entrez","subhectHits")}.
  The first column is the _gene_ nearest the query, by virtue of it
  owning the transcript determined (or chosen by \code{nearest}) to be
  nearest the query. Note that  the nearest gene to a given
  query, in column 3, may not be unique and we arbitrarily chose one as
  done by default by \code{\link[IRanges:nearest-methods]{IRanges::nearest()}}.
  
  The "distance" column is the distance from the query to the 5' end of the
  nearest transcript, so may be different from the distance computed by
  \code{nearest} to that transcript, as a range.

\item{name}{Symbol of nearest gene}
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
\item{insideDistance}{distance past 5' end of gene}
\item{exonnumber}{which exon}
\item{nexons}{number of exons}
\item{UTR}{a factor with levels \code{c("inside transcription region", "5' UTR",
	"overlaps 5' UTR", "3'UTR", "overlaps 3'UTR", "covers transcription region")}}
\item{strand}{"+" or "-"}
\item{geneL}{the gene length}
\item{codingL}{the coding length}
\item{Entrez}{Entrez ID of closest gene}
\item{subjectHits}{Index in subject of hit}
}
\author{
Harris Jaffee, Peter Murakami and Rafael A. Irizarry
}
\seealso{
\code{annotateNearest}, \code{annotateTranscripts}
}
\examples{
\dontrun{
    islands=makeGRangesFromDataFrame(read.delim("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg19.txt")[1:100,])
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
    tab<- matchGenes(islands,genes)
}
}
