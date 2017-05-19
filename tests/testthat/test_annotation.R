context('Annotation-related functions')


## Create toy data with all cases
x <- data.frame(start = c(2, 4, 7, 10, 13, 17, 19, 22, 24), end = c(3, 5, 9, 12, 16, 18, 20, 23, 26), chr = 'chr1')
subject <- data.frame(start = c(1, 4, 7, 11, 14, 21, 25), end = c(3, 6, 9, 12, 15, 22, 27), chr = 'chr1')

## Test cases and types of input
test_that('Annotate Nearest', {
    expect_that(annotateNearest(x, subject, select = 'all'), throws_error())
    expect_that(annotateNearest(x, subject), equals(annotateNearest(x, makeGRangesFromDataFrame(subject))))
    expect_that(annotateNearest(x, subject), equals(annotateNearest(makeGRangesFromDataFrame(x), subject)))
    expect_that(annotateNearest(x, subject)$distance, equals(rep(c(0, -2, 1, 0), c(5, 1, 1, 2))))
    expect_that(annotateNearest(x, subject)$subjectHits, equals(c(1:5, 5:6, 6:7)))
    expect_that(annotateNearest(x, subject)$type, equals(rep(c('inside', 'cover', 'disjoint', 'overlap'), c(2, 3, 2, 2))))
    expect_that(annotateNearest(x, subject)$amountOverlap, equals(rep(c(NA, -1, 2), c(7, 1, 1))))
    expect_that(annotateNearest(x, subject)$insideDist, equals(rep(c(0, NA), c(2, 7))))
    expect_that(annotateNearest(x, subject)$size1, equals(rep(c(2, 3, 4, 2, 3), c(2, 2, 1, 3, 1))))
    expect_that(annotateNearest(x, subject)$size2, equals(rep(c(3, 2, 3), c(3, 5, 1))))
})




## Human (UCSC hg19 knownGene) case
# library('TxDb.Hsapiens.UCSC.hg19.knownGene')
# genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

## Canis Familiaris case
library('GenomicFeatures')
system.time(can_txdb <- makeTxDbFromUCSC("canFam3", "refGene"))
#    user  system elapsed 
#  70.588   1.002 123.785 

## Using a smaller set is not significantly faster
# system.time(can_txdb <- makeTxDbFromUCSC('canFam3', 'refGene', transcript_ids = c('NM_001193298', 'NM_001002949', 'NM_001080724', 'NM_001003268', 'NM_001252259', 'NM_001012395', 'NM_001286958', 'NM_001284493', 'NM_001122602', 'NM_001003021')))
#   user  system elapsed 
# 61.193   1.035 117.035 

can_genes <- annotateTranscripts(can_txdb)
can_genes_NA <- annotateTranscripts(can_txdb, NA)

test_that('Annotate Transcripts', {
    expect_that(attributes(can_genes)$description, equals('annotatedTranscripts'))
    expect_that(attributes(can_genes_NA)$description, equals('annotatedTranscripts'))
    expect_that(colnames(mcols(can_genes)), equals(colnames((mcols(can_genes_NA)))))
    expect_that(ranges(can_genes), equals(ranges(can_genes_NA)))
    expect_that(can_genes_NA$Gene, equals(Rle(NA, length(can_genes_NA))))
    expect_that(can_genes_NA$Refseq, equals(Rle(NA, length(can_genes_NA))))
    expect_that(annotateTranscripts(can_txdb, by = 'gene', annotationPackage = 'DoesNotExist', requireAnnotation = TRUE), throws_error())
    expect_that(sum(can_genes$CSE < can_genes$CSS, na.rm = TRUE), equals(0))
})



## Test matchGenes()

test <- head(can_genes)
strand(test) <- ifelse(strand(test) == "+", "-", "+")
test2 <- head(can_genes)
strand(test2) <- '*'
test3 <- resize(head(can_genes), 10, fix = 'center')
matchGenes(test, can_genes)


matched <- matchGenes(head(can_genes), can_genes)
matched_rev <- matchGenes(test, can_genes)
test_that('Match genes', {
    expect_that(matched, equals(matchGenes(test2, can_genes)))
    expect_that(matchGenes(test3, can_genes)$distance, equals(matched$distance + floor((width(test) - 10) / 2)))
    expect_that(matchGenes(head(can_genes), can_genes, type = 'fiveprime')$subjectHits, equals(1:6))
    expect_that(matchGenes(head(can_genes_NA), can_genes_NA)[, 3:ncol(matched)], equals(matched[, 3:ncol(matched)]))
    expect_that(is.na(matched_rev$codingL[which(matched_rev$Geneid == 104797479)]), equals(TRUE))
    expect_that(sum(matched_rev$description %in% c('downstream', 'upstream')), equals(6))
    expect_that(matchGenes(data.frame(start = start(head(can_genes)), end = end(head(can_genes)), chr = seqnames(head(can_genes)), strand = strand(head(can_genes))), can_genes), equals(matched))
})


## Testing https://github.com/ririzarr/bumphunter/issues/4
library('GenomicRanges')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
query <- GRanges('chr21', IRanges(22115534, 22115894), strand = "*")
genes <- annotateTranscripts(txdb = txdb)
res <- matchGenes(x = query, subject = genes)

test_that('NA bug', {
    expect_that(as.character(res$UTR), equals('inside transcription region'))
    expect_that(res$insideDistance, equals(0))
    expect_that(as.character(res$description), equals('inside exon'))
})


## Test using a GENCODE v25 txdb object
## Details from https://github.com/rafalab/bumphunter/issues/15

library('rtracklayer')
## Get the raw data
gr <- import('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz')

## Subset and add the chromosome length info
gr_small <- keepSeqlevels(gr, c('chrY', 'chrM'), pruning.mode = 'tidy')
hg38_chrominfo <- fetchExtendedChromInfoFromUCSC("hg38")
new_info <- hg38_chrominfo$UCSC_seqlength[match(seqlevels(gr_small),
    hg38_chrominfo$UCSC_seqlevel)]

## Create a small TxDb object
txdb <- makeTxDbFromGRanges(gr_small)

ann <- annotateTranscripts(txdb, annotationPackage = 'org.Hs.eg.db',
    mappingInfo = list('column' = 'ENTREZID', 'keytype' = 'ENSEMBL',
    'multiVals' = 'first'), simplifyGeneID = TRUE)

## Annotate some dummy regions
genes_gencode <- matchGenes(ann[which(ann$Gene == 'CD99')], ann)

## Check the reuslt
test_that('Gencode v25 genes', {
    expect_equal(unique(gr_small[grepl(genes_gencode$Geneid[1], gr_small$gene_id)]$gene_name), unique(genes_gencode$name))
    expect_equal('CD99', unique(genes_gencode$name))
})
