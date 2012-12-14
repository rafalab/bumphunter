\name{TT}
\docType{data}
\alias{TT}
\title{HG19 Transcripts}
\description{
Database of transcripts associated with "known" hg19 genes,
namely those with Entrez ID gene identifiers associated.
}
\usage{data(TT)}
\format{
  A \code{\link[GenomicRanges]{GRanges}} object with 8 metadata
  columns, "CSS" (coding start), "CSE" (coding end), "Tx" (transcript  
  name), "Entrez" (Entrez ID), "Gene" (gene name), "Refseq" (Refseq
  number), "Nexons" (number of exons), and "Exons" (the exons themselves,
  as \code{IRanges} objects).  Note that CSS is always less than CSE,
  even on minus strands where their "start" and "end" meanings are
  reversed.  Similarly with the exons.
}
\source{the result of \code{known_transcripts()}
}
\keyword{hg19, known transcripts}
