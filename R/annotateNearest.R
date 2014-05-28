annotateNearest <- function(x, subject, annotate=TRUE, ...) {

    # matchGenes
    if (class(subject) == "character") {
	if (subject != "hg19")
		stop("matchGenes: only 'hg19' supported for now")
	return(.matchGenes(object=x, build=subject, ...))
    }

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
	stop("arguments must be both IRanges or GRanges")

    dots <- list(...)
    names.dots <- names(dots)
    if ("select" %in% names.dots && dots$select == "all")
	stop("select=\"all\" not supported (yet)")
    if (type==1 && "ignore.strand" %in% names.dots)
	stop("cannot specify 'ignore.strand' with IRanges")

    # find nearests
    NN <- nearest(x, subject, ...)
    nearest_exists <- !is.na(NN)

    # full return value
    Ret <- matrix(NA, nrow=length(NN), ncol=7)
    colnames(Ret) <- c("dist", "matchIndex", "type", "amountOverlap",
	"insideDist", "size1", "size2")
    Ret <- data.frame(Ret, stringsAsFactors=FALSE)

    # queries with a nearest target
    x <- x[nearest_exists]
    # corresponding nearest targets
    y <- subject[NN[nearest_exists],]

    # partial value for when nearest exists
    ret <- matrix(NA, nrow=length(x), ncol=7)
    colnames(ret) <- c("dist", "matchIndex", "type", "amountOverlap",
	"insideDist", "size1", "size2")
    ret[,"matchIndex"] <- NN[nearest_exists]
    ret <- data.frame(ret, stringsAsFactors=FALSE)

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

    if (!annotate) {
	Ret[nearest_exists,] <- ret
	return(as.matrix(Ret[, c("dist", "matchIndex")]))
    }

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
    
    Ret[nearest_exists,] <- ret
    Ret
}

# for use by pointMatch
regionMatch <- annotateNearest

# for backward compatibility
matchGenes <- function(x, build="hg19", ...) annotateNearest(x=x, subject=build, ...)
