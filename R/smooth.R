smoother <- function(y, x=NULL, cluster, weights=NULL,
                     smoothFunction, verbose=TRUE, ...) {
    if(is.null(dim(y)))
        y <- matrix(y, ncol=1) ##need to change this to be more robust
    if(!is.null(weights) && is.null(dim(weights)))
        weights <- matrix(weights, ncol=1)
    cores <- getDoParWorkers()
    Indexes <- split(seq(along=cluster), cluster)
    ## Split Indexes into chunks that will be distributed among the workers
    IndexesChunks <- lapply(as.list(isplitVector(Indexes, chunks=cores)), unlist) # not super-elegant

    # make binding of 'idx' explicit
    foreach_object <- foreach(idx=IndexesChunks, .packages = "bumphunter")
    foreach_argname <- foreach_object$argnames[1]
    ret <- foreach_object %dorng% {
	idx <- get(foreach_argname)
        smoothFunction(y=y[idx,], x=x[idx], cluster=cluster[idx],
                       weights=weights[idx,], verbose=verbose, ...)
    }

    attributes(ret)[["rng"]] <- NULL
    ## Paste together results from different workers
    ret <- reduceIt(ret)
    ret$smoother <- ret$smoother[1]
    return(ret)  
}

locfitByCluster <- function(y, x = NULL, cluster, weights = NULL,
                            minNum = 7, bpSpan = 1000, minInSpan = 5,
                            verbose = TRUE) {
    ## Depends on limma
    ## bpSpan is in basepairs
    ## assumed x are ordered
    ## if y is vector change to matrix
    if(is.null(dim(y)))
        y <- matrix(y, ncol = 1) 
    if(!is.null(weights) && is.null(dim(weights)))
        weights <- matrix(weights, ncol = 1)
    if(is.null(x))
        x <- seq(along = y)
    if(is.null(weights))
        weights <- matrix(1, nrow = nrow(y), ncol = ncol(y))
    Indexes <- split(seq(along = cluster), cluster)
    clusterL <- sapply(Indexes, length)
    smoothed <- rep(TRUE, nrow(y))
    
    for(i in seq(along=Indexes)) {
        if(verbose) if(i %% 1e4 == 0) cat(".")
        Index <- Indexes[[i]]
        if(clusterL[i] >= minNum) {
            nn <- minInSpan / length(Index)
            for(j in 1:ncol(y)) {
                sdata <- data.frame(pos = x[Index],
                                    y = y[Index, j],
                                    weights = weights[Index,j])
                fit <- locfit(y ~ lp(pos, nn = nn, h = bpSpan), data = sdata,
                              weights = weights, family = "gaussian", maxk = 10000)
                pp <- preplot(fit, where = "data", band = "local",
                              newdata = data.frame(pos = x[Index]))
                y[Index, j] <- pp$trans(pp$fit)
            }
        } else {
            y[Index,] <- NA
            smoothed[Index] <- FALSE
        }
    }                
    return(list(fitted=y, smoothed=smoothed, clusterL=clusterL, smoother="locfit"))
}

loessByCluster <- function(y, x=NULL, cluster, weights= NULL,
                         bpSpan = 1000, minNum=7, minInSpan=5,
                         maxSpan=1, verbose=TRUE) {
    ## Depends on limma
    ## bpSpan is in basepairs
    ## assumed x are ordered
    ## if y is vector change to matrix
    if(is.null(dim(y)))
        y <- matrix(y,ncol=1) ##need to change this to be more robust
    if(!is.null(weights) && is.null(dim(weights)))
        weights <- matrix(weights,ncol=1)
    
    if(is.null(x))
        x <- seq(along=y)
    if(is.null(weights))
        weights <- matrix(1, nrow=nrow(y), ncol=ncol(y))
    Indexes <- split(seq(along=cluster), cluster)
    clusterL <- sapply(Indexes, length)
    spans <- rep(NA, length(Indexes))
    smoothed <- rep(TRUE, nrow(y))
    
    for(i in seq(along=Indexes)) {
        if(verbose) if(i %% 1e4 == 0) cat(".")
        Index <- Indexes[[i]]
        if(clusterL[i] >= minNum) {
            span = bpSpan/median(diff(x[Index]))/length(Index) # make a span
            if(span > maxSpan) span <- maxSpan
            spans[i] <- span
            if(span*length(Index)>minInSpan){
                ##this can be parallelized
                for(j in 1:ncol(y)){
                    y[Index,j] <- limma::loessFit(y[Index,j], x[Index], span = span,
                                                  weights = weights[Index,j])$fitted
                }
            } else{
                y[Index,] <- NA
                smoothed[Index] <- FALSE
            }
        } else{
            y[Index,] <- NA
            spans[i] <- NA
            smoothed[Index] <- FALSE
        }
    }
    return(list(fitted=y, smoothed=smoothed, spans=spans, clusterL=clusterL, smoother="loess"))
}


runmedByCluster <- function(y, x=NULL, cluster, weights=NULL,
                            k=5, endrule="constant", verbose=TRUE) {
    ## we dont use chr and pos byt
    ## Depends on limma
    ## bpSpan is in basepairs
    ## assumed y are ordered
    ## x is ingored (for now... should be used for ties)
    ## weight is ingored too. x and weight are included for consistency
    ## in parameters
    if(is.null(dim(y)))
        y <- matrix(y, ncol=1) ##need to change this to be
    
    Indexes <- split(seq(along=cluster), cluster)
    clusterL <- sapply(Indexes, length)
    spans <- rep(k, length(Indexes)) ##using spans for consistency with loessByClusters
    smoothed <- rep(TRUE, nrow(y))
    ## We can parallelize here
    for(i in seq(along=Indexes)) {
        if(verbose) if(i %% 1e4 == 0) cat(".")
        Index <- Indexes[[i]]
        if(clusterL[i] >= k) {
            for(j in 1:ncol(y)){
                y[Index,j] <- runmed(y[Index,j], k=k, endrule=endrule)
            }
        } else {
            y[Index] <- NA
            smoothed[Index] <- FALSE
        }
    }
    return(list(fitted=y, smoothed=smoothed, spans=spans, clusterL=clusterL, smoother="runmed"))
}


reduceIt <- function(x, elem, bind=rbind) {
    if(missing(elem)) elem <- names(x[[1]])
    ret <- lapply(elem, function(el) {
        xx <- lapply(x, "[[", el)
        if (is.matrix(xx[[1]])) return(do.call(bind, xx))
        else if (is.vector(xx[[1]])) return(do.call("c", xx))
        else stop("reduce can only handle matrices or vectors")
    })
    names(ret) <- elem
    ret
}
