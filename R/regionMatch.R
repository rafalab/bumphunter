regionMatch <- function(x, subject, verbose=TRUE) {
    ret <- matrix(NA,nrow=length(x),ncol=7)
    colnames(ret) <- c("dist", "matchIndex", "type", "amountOverlap",
                       "insideDist", "size1", "size2")
    sp1 <- split(1:length(x), as.character(seqnames(x)))
    sp2 <- split(1:length(subject), as.character(seqnames(subject)))
    for(i in names(sp1)) {
        if(verbose) cat(i," ")
        inds1 <- sp1[[i]]
        if(i %in% names(sp2)) {
          inds2 <- sp2[[i]]
          ret[inds1,"matchIndex"] <- inds2[
               nearest(ranges(x[inds1,]),ranges(subject[inds2,]))]
        } else ret[inds1,"matchIndex"] = NA
      }
    if(verbose) cat("\n")
    y <- subject[ret[,"matchIndex"],]
    ret <- data.frame(ret, stringsAsFactors=FALSE)
    
    ret$type[(start(x)>start(y) & end(x)<=end(y)) |
             (start(x)>=start(y) & end(x)<end(y))] <- "inside"
    ret$type[start(x)<=start(y) & end(x)>=end(y)] <- "cover"
    ret$type[start(x)>end(y)] <- "disjointR"
    ret$type[end(x)<start(y)] <- "disjointL"
    ret$type[is.na(ret$matchIndex)] <- "disjoint"
    ret$type[(start(x)>start(y) & start(x)<=end(y)) & end(x)>end(y)] <- "overlapR"    
    ret$type[start(x)<start(y) & (end(x)>=start(y) & end(x)<end(y))] <- "overlapL"

    ret$dist <- 0
    ret$dist[ret$type=="disjoint"] <- NA
    ret$dist[ret$type=="disjointR"] <- end(y)[ret$type=="disjointR"] - start(x)[ret$type=="disjointR"]
    ret$dist[ret$type=="disjointL"] <- start(y)[ret$type=="disjointL"] - end(x)[ret$type=="disjointL"]
    ret$amountOverlap[ret$type=="overlapR"] <- -1*(end(y)[ret$type=="overlapR"]-start(x)[ret$type=="overlapR"]+1)
    ret$amountOverlap[ret$type=="overlapL"] <- end(x)[ret$type=="overlapL"]-start(y)[ret$type=="overlapL"]+1
    ret$type[ret$type%in%c("disjointR","disjointL")] <- "disjoint"
    ret$type[ret$type%in%c("overlapR","overlapL")] <- "overlap"

    ## insideDist column:
    insideIndex <- ret$type=="inside" #no missing ret$type at this point
    tmp0 <- cbind(end(x)[insideIndex] - end(y)[insideIndex],
                  start(x)[insideIndex] - start(y)[insideIndex])
    tmp <- apply(abs(tmp0),1,which.min)
    tmpinsidedist <- tmp0[,1]
    tmpIndex <- tmp==2
    tmpinsidedist[tmpIndex] <- tmp0[tmpIndex,2]
    ret$insideDist[insideIndex] <- tmpinsidedist

    ## size1 and size2 columns:
    ret$size1 <- end(x) -start(x) +1
    ret$size2 <- end(y)-start(y)+1
    
    ret
}
