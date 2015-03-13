\name{annotateTranscripts}
\alias{annotateTranscripts}
\title{
Annotate transcripts
}
\description{
Annotate transcripts
}
\usage{
annotateTranscripts(txdb, annotationPackage = NULL, by = c("tx","gene"), codingOnly=FALSE, verbose = TRUE, requireAnnotation = FALSE)
}
\arguments{
  \item{txdb}{
    A \code{TxDb} database object such as \code{TxDb.Hsapiens.UCSC.hg19.knownGene}
  }
  \item{annotationPackage}{
    An annotation data package from which to obtain gene/transcript
    annotation. For example \code{org.Hs.eg.db}. If none is provided the
    function tries to infer it from \code{species(txdb)} and if it can't it
    proceeds without annotation unless \code{ requireAnnotation = TRUE}.
  }
  \item{by}{
    Should we create a \code{GRanges} of transcripts (\code{tx}) or genes (\code{gene}).
  }
  \item{codingOnly}{
    Should we exclude all the non-coding transcripts.
  }
  
  \item{verbose}{
    logical value. If 'TRUE', it writes out some messages
    indicating progress. If 'FALSE' nothing should be printed.
  }
  \item{requireAnnotation}{
    logical value. If 'TRUE' function will stop if no annotation package
    is successfully loaded.
  }
}
\details{
  This function prepares a \code{GRanges} for the \code{\link{matchGenes}}
  function. It adds information and in particular adds exons information
  to each gene/transcript.
}
\value{
  A \code{GRanges} object with an attribute \code{description} set to
  \code{annotatedTranscripts}. The following columns are added.
  \code{seqinfo} is the information returned by
  \code{\link[GenomeInfoDb]{seqinfo}}, \code{CSS} is the coding region
  start, \code{CSE} is the coding region end, \code{Tx} is the transcript
  ID used in TxDb, \code{Entrez} is the Entrez ID, \code{Gene} is the gene
  symbol, \code{Refseq} is the RefSeq annotation, \code{Nexons} is the
  number of exons, \code{Exons} is an \code{IRanges} with the exon information.
}
\author{
  Harris Jaffee and Rafael A. Irizarry
}

\seealso{
  \code{\link{matchGenes}}
}
\examples{
\dontrun{
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    ##and to avoid guessing the annotation package:
    genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene,annotation="org.Hs.eg.db")
    
}
}