.getSegments_old <- function(x, factor, cutoff=quantile(abs(x),0.99), verbose=TRUE){
    Indexes <- split(seq(along=x),factor)
    regionID <- vector("numeric", length(x))
    LAST <- 0
  
    segmentation <- vector("numeric", length(x))
    type <- vector("numeric", length(x))

    for (i in seq(along = Indexes)) {
        if (verbose) if (i%%1000 == 0) cat(".")
        Index <- Indexes[[i]]
        y <- x[Index]
        z <- sign(y) * as.numeric(abs(y) > cutoff)
        w <- cumsum(c(1, diff(z) != 0)) + LAST
        segmentation[Index] <- w
        type[Index] <- z
        LAST <- max(w)
    }
    ##add a vector of the pns
    res <- list(upIndex = split(which(type>0), segmentation[type>0]),
                dnIndex = split(which(type<0), segmentation[type<0]),
                zeroIndex = split(which(type==0), segmentation[type==0]))
    names(res[[1]]) <- NULL
    names(res[[2]]) <- NULL
    names(res[[3]]) <- NULL
    return(res)
}

getSegments <- function(x, f = NULL, cutoff=quantile(abs(x), 0.99), assumeSorted = FALSE, verbose=FALSE){
    if(is.null(f))
        f <- rep(1L, length(x))
    stopifnot(length(x) == length(f))
    stopifnot(length(cutoff) <= 2)
    if(is.character(f))
        f <- as.factor(f)
    if(is.numeric(f))
        f <- as.integer(f)
    stopifnot(is.factor(f) || is.integer(f))
    if(length(cutoff) == 1)
        cutoff <- c(-cutoff, cutoff)
    cutoff <- sort(cutoff)

    if(assumeSorted) {
        needReorder <- is.unsorted(f)
        if(needReorder) {
            od <- order(f)
            x <- x[od]
            f <- f[od]
        }
    }
        
    if(verbose) message("getSegments: segmenting")
    Indexes <- split(seq(along=x), f)
    direction <- as.integer(x >= cutoff[2])
    direction[x <= cutoff[1]] <- -1L

    ## We now need to convert f into cid
    if(verbose) message("getSegments: splitting")
    segments <- cumsum(c(1, diff(direction) != 0)) +
        cumsum(c(1, diff(f) != 0))
    names(segments) <- NULL

    res <- list(upIndex = split(which(direction>0), segments[direction>0]),
                dnIndex = split(which(direction<0), segments[direction<0]),
                zeroIndex = split(which(direction==0), segments[direction==0]))
    names(res[[1]]) <- NULL
    names(res[[2]]) <- NULL
    names(res[[3]]) <- NULL

    if(assumeSorted && needReorder) {
        res <- lapply(res, function(sp) lapply(sp, function(xx) od[xx]))
    }
    res
}

clusterMaker <- function(chr, pos, assumeSorted = FALSE, maxGap=300){
    
    nonaIndex <- which(!is.na(chr) & !is.na(pos))
    Indexes <- split(nonaIndex, chr[nonaIndex])
    clusterIDs <- rep(NA, length(chr))
    LAST <- 0
    for(i in seq(along = Indexes)){
        Index <- Indexes[[i]]
        x <- pos[Index]
        if(!assumeSorted){
            Index <- Index[order(x)]
            x <- pos[Index]
        }
        y <- as.numeric(diff(x) > maxGap)
        z <- cumsum(c(1, y))
        clusterIDs[Index] <- z + LAST
        LAST <- max(z) + LAST
    }
    clusterIDs
}

##you can pass cutoff through the ...
regionFinder <- function(x, chr, pos, cluster=NULL, y=x, summary=mean, ind=seq(along=x),
                         order=TRUE, oneTable=TRUE, maxGap=300, cutoff=quantile(abs(x), 0.99),
                         assumeSorted = FALSE, verbose = TRUE){
    if(any(is.na(x[ind]))){
        warning("NAs found and removed. ind changed.")
        ind <- intersect(which(!is.na(x)),ind)
    } 
    if(is.null(cluster))
        cluster <- clusterMaker(chr, pos, maxGap=maxGap, assumeSorted = assumeSorted)
    Indexes <- getSegments(x = x[ind], f = cluster[ind], cutoff = cutoff,
                           assumeSorted = assumeSorted, verbose = verbose)
    clusterN <- table(cluster)[as.character(cluster)]
    
    res <- vector("list",2)
    for(i in 1:2){
        res[[i]] <- data.frame(chr = sapply(Indexes[[i]], function(Index) chr[ind[Index[1]]]),
                               start = sapply(Indexes[[i]], function(Index) min(pos[ind[Index]])),
                               end = sapply(Indexes[[i]], function(Index) max(pos[ind[Index]])),
                               value = sapply(Indexes[[i]], function(Index) summary(y[ind[Index]])),
                               area = sapply(Indexes[[i]], function(Index) abs(sum(y[ind[Index]]))),
                               cluster = sapply(Indexes[[i]], function(Index) cluster[ind[Index]][1]),
                               indexStart = sapply(Indexes[[i]], function(Index) min(ind[Index])),
                               indexEnd = sapply(Indexes[[i]], function(Index) max(ind[Index])))
        res[[i]]$L <- res[[i]]$indexEnd - res[[i]]$indexStart+1
        res[[i]]$clusterL <- sapply(Indexes[[i]], function(Index) clusterN[ind[Index]][1])
    }
    names(res) <- c("up","dn")
    if(order & !oneTable){
        if(nrow(res$up)>0) res$up <- res$up[order(-res$up$area),]
        if(nrow(res$dn)>0) res$dn <- res$dn[order(-res$dn$area),]
    }
    if(oneTable){
        res <- rbind(res$up,res$dn)
        if(order & nrow(res)>0) res <- res[order(-res$area),]
    }
    return(res)
}
