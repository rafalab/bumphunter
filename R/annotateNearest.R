dataframe_to_GRanges <- function(D) {
        require(GenomicRanges)
        if (is(D, "GRanges"))
                return(D)
        names <- names(D)
        seqnames <- if ("chr" %in% names) {
                tmp <- D$chr
                if (is.character(tmp))
                        tmp
                else
                        # the 'chr' variable may be numeric
                        paste("chr", tmp, sep="")
            }
            else
                D$seqnames
        ranges <- if ("start" %in% names)
                IRanges(start=D$start, end=D$end)

        strands <- if("strand" %in% names) D$strand else "*"
        GRanges(seqnames=seqnames, ranges=ranges, strand=strands)
}

dataframe_to_IRanges <- function(D) {
        require(IRanges)
        if (is(D, "IRanges"))
                return(D)
        IRanges(start=D$start, end=D$end)
}

annotateNearest <- function(x, subject, pointMatch=FALSE, ...) {

    if (class(x) == "data.frame") {
	names <- names(x)
	if ("chr" %in% names || "seqnames" %in% names) {
		x <- dataframe_to_GRanges(x)
		subject <- dataframe_to_GRanges(subject)
	}
	else {
		x <- dataframe_to_IRanges(x)
		subject <- dataframe_to_IRanges(subject)
	}
    }

    # validate argument types
    type <- rep(NA, 2)
    if (is(x, "IRanges"))
	type[1] <- 1
    else if (is(x, "GRanges"))
	type[1] <- 2
    if (is(subject, "IRanges"))
	type[2] <- 1
    else if (is(subject, "GRanges"))
	type[2] <- 2
    type <- prod(type)
    if(is.na(type) || type==2)
	stop("annotatNearest: arguments must be both IRanges or GRanges")
    if (type==1)
	require(IRanges)
    else
	require(GenomicRanges)

    ret <- matrix(NA, nrow=length(x), ncol=7)
    colnames(ret) <- c("dist", "matchIndex", "type", "amountOverlap",
	"insideDist", "size1", "size2")

    dots <- list(...)
    names.dots <- names(dots)
    if ("select" %in% names.dots && dots$select == "all")
	stop("regionMatch: select=\"all\" not supported (yet)")
    if (type==1 && "ignore.strand" %in% names.dots)
	stop("regionMatch: cannot specify 'ignore.strand' with IRanges")

#cat("calling nearest... ")
    NN <- nearest(x, subject, ...)
    if (any(is.na(NN)))
	stop("regionMatch: nearest returned NAs")
#cat("Done\n")

    ret[,"matchIndex"] <- NN
    ret <- data.frame(ret, stringsAsFactors=FALSE)
    y <- subject[NN,]

    x.start <- start(x)
    y.start <- start(y)
    x.end <- end(x)
    y.end <- end(y)

    ret$type[(x.start > y.start & x.end <= y.end) |
             (x.start >= y.start & x.end < y.end)] <- "inside"
    ret$type[x.start <= y.start & x.end >= y.end] <- "cover"
    ret$type[x.start > y.end] <- "disjointR"
    ret$type[x.end < y.start] <- "disjointL"
    ret$type[x.start > y.start & x.start <= y.end & x.end > y.end] <- "overlapR"    
    ret$type[x.start < y.start & x.end >= y.start & x.end < y.end] <- "overlapL"

    ret$dist <- 0L

    i <- ret$type == "disjointR"
    ret$dist[i] <- y.end[i] - x.start[i]

    i <- ret$type == "disjointL"
    ret$dist[i] <- y.start[i] - x.end[i]

    if (pointMatch)
	return(as.matrix(ret[, c("dist", "matchIndex")]))

    i <- ret$type == "overlapR"
    ret$amountOverlap[i] <- -1L * (y.end[i] - x.start[i] + 1L)

    i <- ret$type == "overlapL"
    ret$amountOverlap[i] <- x.end[i] - y.start[i] + 1L

    ret$type[ret$type%in%c("disjointR", "disjointL")] <- "disjoint"
    ret$type[ret$type%in%c("overlapR", "overlapL")] <- "overlap"

    ## insideDist column:
    i <- ret$type == "inside"	#no missing ret$type at this point
    tmp0 <- cbind(x.end[i] - y.end[i], x.start[i] - y.start[i])
    tmp <- apply(abs(tmp0), 1, which.min)
    tmpinsidedist <- tmp0[, 1]
    tmpIndex <- tmp == 2
    tmpinsidedist[tmpIndex] <- tmp0[tmpIndex, 2]
    ret$insideDist[i] <- tmpinsidedist

    ## size1 and size2 columns:
    ret$size1 <- x.end - x.start + 1L
    ret$size2 <- y.end - y.start + 1L
    
    ret
}
