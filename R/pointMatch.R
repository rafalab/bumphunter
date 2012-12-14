pointMatch <- function(object1, object2, col1=2, col2=2, verbose) {
	ans <- matrix(NA, nrow=dim(object1)[1], ncol=2)
	colnames(ans) <- c("dist", "matchIndex")
	x.chrs <- object1$chr
	if (is.numeric(x.chrs))
		x.chrs <- paste("chr", x.chrs, sep="")
	y.chrs <- object2$chr
	if (is.numeric(y.chrs))
		y.chrs <- paste("chr", y.chrs, sep="")
	II <- which(x.chrs %in% y.chrs)
	if (length(II) == 0)
		return(ans)
	x <- object1[II, col1]
	x.chrs <- x.chrs[II]
	x <- GRanges(seqnames=x.chrs, ranges=IRanges(x, x))
	y <- object2[, col2]
	y <- GRanges(seqnames=y.chrs, ranges=IRanges(y, y))
	ans[II, ] <- regionMatch(x, y, annotate=FALSE)
	ans
}
