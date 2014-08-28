known_transcripts <- function(species="Hsapiens", build="hg19") {
	#package <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
	package <- paste("TxDb", species, "UCSC", build, "knownGene", sep=".")
	Org <- paste("org", substr(species, 1, 2), sep=".")

cat("building known_transcripts for:\t")
cat(species, "/", build, "\n")

cat("loading the knownGene database\n")
	library(package=package, character.only=TRUE)
	txdb_package_description <- packageDescription(package)
	the_db <- get(package)

cat("... ")
	# transcripts-by-gene
	tt <- transcriptsBy(the_db, by="gene")
	seqinfo <- seqinfo(tt)

	chr <- unlist(CharacterList(seqnames(tt)))
	strand <- unlist(CharacterList(strand(tt)))
	#
	RR <- ranges(tt)
	entrez <- names(RR)
	# starts / ends (left/right endpoints):
	TSS <- unlist(start(RR))
	TSE <- unlist(end(RR))

	# the transcript names
	LL <- tt@unlistData@elementMetadata@listData
	txx <- LL$tx_name
cat(length(txx), "transcripts in", length(tt), "genes\n")

cat("exons by transcript\n")
	# exons-by-transcript
	xx <- exonsBy(the_db, by="tx", use.names=TRUE)
	# exon groups of interest in transcript order,
	# and sorted "left to right" within each group
	Exons <- reduce(ranges(xx)[txx])
	Nexons <- elementLengths(Exons)
	
cat("coding by transcript\n")
	# coding-regions-by-transcript
	cds <- cdsBy(the_db, by="tx", use.names=TRUE)
	with_coding <- which(txx %in% names(cds))
	
	# take only the left- and right-most endpts:
	cds_ranges <- range(cds[txx[with_coding]])
	
	# starts and ends (really left / right endpts):
	CSS <- CSE <- integer(length(txx))
	CSS[with_coding] <- unlist(start(cds_ranges))
	CSE[with_coding] <- unlist(end(cds_ranges))
	CSS[-with_coding] <- NA
	CSE[-with_coding] <- NA
	
cat("entrez ids to gene symbol and refseq id\n")
	# Get the symbols associated with the entrez ids
	org_package <- paste(Org, "eg", "db", sep=".")
	#library(org.Hs.eg.db)
	library(org_package, character.only=TRUE)
	org_package_description <- packageDescription(org_package)

	#map <- org.Hs.egSYMBOL
	map <- get(paste(Org, "egSYMBOL", sep="."))

	which <- mappedkeys(map)
	#symbols <- unlist(as.list(org.Hs.egSYMBOL[which]))
	symbols <- unlist(as.list(map[which]))
	genes <- symbols[entrez]		# a handful of NA

	#map <- org.Hs.egREFSEQ
	map <- get(paste(Org, "egREFSEQ", sep="."))
	which <- mappedkeys(map)
	#symbols <- sapply(as.list(org.Hs.egREFSEQ[which]),
	symbols <- sapply(as.list(map[which]), paste, collapse=" ")
	refseq <- symbols[entrez]	# a handful of NA

	TT <- elementLengths(tt)
	Entrez <- Rle(entrez, TT)
	Gene <- Rle(genes, TT)
	Refseq <- Rle(refseq, TT)
	Tx <- txx
	
	# return a list containing the two package descriptions and
	# the associated known transcripts object
cat("formatting the table\n")
	list(txdb=txdb_package_description, org=org_package_description,
	# the reference transcripts, with annotation
	transcripts=GRanges(ranges=IRanges(start=TSS, end=TSE),
		seqnames=chr, strand=strand,
		# general sequence info
		seqinfo=seqinfo,
		# transcript metadata
		CSS, CSE, Tx, Entrez, Gene, Refseq, Nexons, Exons)
	)
}
