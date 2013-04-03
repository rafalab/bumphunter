# helper for nearestgene()
dataframe_to_GRanges <- function(D) {
	if (is(D, "GRanges"))
		return(D)
	names <- names(D)
	seqnames <- if ("chr" %in% names) {
		tmp <- D$chr
		if (is.factor(tmp) || is.character(tmp))
			tmp
		else
			# allow numeric
			paste("chr", tmp, sep="")
	    }
	    else
		D$seqnames
	ranges <- if ("start" %in% names)
		IRanges(start=D$start, end=D$end)

	strands <- if("strand" %in% names) D$strand else "*"
	GRanges(seqnames=seqnames, ranges=ranges, strand=strands)
}

## this doesn't work yet
#DataFrameToGRanges <- function(D) {
	#if (is(D, "GRanges"))
		#return(D)
	#suppressWarnings(as(as(as.data.frame(D), "RangedData"), "GRanges"))
#}

dataframe_to_IRanges <- function(D) {
	if (is(D, "IRanges"))
		return(D)
	IRanges(start=D$start, end=D$end)
}


