# helper for nearestgene()
shrink_to_midpoints <- function(x) {
	width <- ifelse(width(x)==0L, 0L, 1L+(end(x)-start(x))%%2L)
	resize(x, fix="center", width=width)
}

# hg19 only, up/down/r not supported
nearestgene <- function(queries, ignore.strand=TRUE, all=FALSE, useMidpts=FALSE){
  q = dataframe_to_GRanges(queries)

  # shrink queries to their mid*point*, if integral, otherwise an interval of length 1
  if (useMidpts)
	q <- shrink_to_midpoints(q)

  # the reference transcript database, TT
cat("nearestgene: loading bumphunter hg19 transcript database\n")
  tt = data(TT, package="bumphunter")
  TT = get(tt)

  # toss null-coding genes ?
  ##annotation = elementMetadata(TT)
  ##no_coding_region = which(is.na(annotation$CSS))
  ##TT = TT[-no_coding_region,]

  # nearest transcripts
cat("finding nearest transcripts...\n")

  if (all) {
    n <- nearest(x=q, subject=TT, ignore.strand=ignore.strand, select="all")
    N <- TT[subjectHits(n)]
    query.hits <- queryHits(n)
    q <- q[query.hits]
    d <- distance(q, N, ignore.strand=ignore.strand)
    strand <- as.vector(strand(N))
    negatives <- (strand=="+" & end(N) < start(q)) | (strand=="-" & end(q) < start(N))
    d[negatives] <- -d[negatives]
    list(q=query.hits, N=N, d=d)
  }
  else {
    n <- nearest(x=q, subject=TT, ignore.strand=ignore.strand)
    N <- TT[n]
    # normal distance (unsigned)
    d <- distance(q, N, ignore.strand=ignore.strand)
    # signed distance, reflecting strand
    strand <- as.vector(strand(N))
    negatives <- (strand=="+" & end(N) < start(q)) | (strand=="-" & end(q) < start(N))
    d[negatives] <- -d[negatives]
    list(N=N, d=d)
  }
}
