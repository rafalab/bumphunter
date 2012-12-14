# defaults to hg19, i.e.
#	TxDb.Hsapiens.UCSC.hg19.knownGene
# and
#	org.Hs.eg.db
known_transcripts <- function() {

	# the org package should be implied by the tx_db package, e.g.
	#	tx_db = "TxDb.Hsapiens.UCSC.hg19.knownGene"
	#	species = strsplit(tx_db, split=".", fixed=TRUE)[[1]][2]
	#	paste("org", substr(species, 1, 2), "eg", "db", sep=".")

	package = "TxDb.Hsapiens.UCSC.hg19.knownGene"

cat("loading the knownGene database\n")
	require(package=package, character.only=TRUE) || stop('need package: ', package)
	the_db = get(package)

cat("... ")
	# transcripts-by-gene
	tt = transcriptsBy(the_db, by="gene")
	seqinfo = seqinfo(tt)

	chr = unlist(CharacterList(seqnames(tt)))
	strand = unlist(CharacterList(strand(tt)))
	#
	RR = ranges(tt)
	entrez = names(RR)
	# starts / ends (left/right endpoints):
	TSS = unlist(start(RR))
	TSE = unlist(end(RR))

	# the transcript names
	LL = tt@unlistData@elementMetadata@listData
	txx = LL$tx_name
cat(length(txx), "transcripts in", length(tt), "genes\n")

cat("exons by transcript\n")
	# exons-by-transcript
	xx = exonsBy(the_db, by="tx", use.names=TRUE)
	# exon groups of interest in transcript order,
	# and sorted "left to right" within each group
	Exons = reduce(ranges(xx)[txx])
	Nexons = elementLengths(Exons)
	
cat("coding by transcript\n")
	# coding-regions-by-transcript
	cds = cdsBy(the_db, by="tx", use.names=TRUE)
	with_coding = which(txx %in% names(cds))
	
	# take only the left- and right-most endpts:
	cds_ranges = range(cds[txx[with_coding]])
	
	# starts and ends (really left / right endpts):
	CSS = CSE = integer(length(txx))
	CSS[with_coding] = unlist(start(cds_ranges))
	CSE[with_coding] = unlist(end(cds_ranges))
	CSS[-with_coding] = NA
	CSE[-with_coding] = NA
	
cat("entrez ids to gene symbol and refseq id\n")
	# Get the symbols associated with the entrez ids
	require(org.Hs.eg.db) || stop("need package org.Hs.eg.db")

	map = org.Hs.egSYMBOL
	symbols = unlist(as.list(org.Hs.egSYMBOL[mappedkeys(map)]))
	genes = symbols[entrez]	# a handful of these are NA

	map = org.Hs.egREFSEQ
	symbols = unlist(as.list(org.Hs.egREFSEQ[mappedkeys(map)]))
	refseq = symbols[entrez]	# a handful of these are NA

	TT = elementLengths(tt)
	Entrez = Rle(entrez, TT)
	Gene = Rle(genes, TT)
	Refseq = Rle(refseq, TT)
	Tx = txx
	
cat("formatting the table\n")
	# the reference transcripts, with annotation
	GRanges(ranges=IRanges(start=TSS, end=TSE),
		seqnames=chr, strand=strand,
		# general sequence info
		seqinfo=seqinfo,
		# transcript metadata
		CSS, CSE, Tx, Entrez, Gene, Refseq, Nexons, Exons)
}
