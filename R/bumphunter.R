setMethod("bumphunter", signature(object = "matrix"),
          function(object, design, chr=NULL, pos, cluster=NULL,
                   coef=2, cutoff=NULL, pickCutoff=FALSE, cutoffQ=0.99,
                   maxGap=500, smooth=FALSE, smoothFunction=locfitByCluster,
                   useWeights=FALSE, B=1000, verbose=TRUE, ...){
            if(missing(design)) stop("design must be specified")
            if(missing(pos)) stop("If object is a matrix, pos must be specified")
            bumphunterEngine(object, design=design, chr=chr, pos,
                             cluster=cluster,
                             coef=coef,cutoff=cutoff,pickCutoff=pickCutoff,
                             cutoffQ=cutoffQ,
                             maxGap=maxGap,smooth=smooth,
                             smoothFunction=smoothFunction,
                             useWeights=useWeights,
                             B=B, verbose=verbose, ...)
          })

bumphunterEngine <- function(mat, design, chr=NULL, pos, cluster=NULL,
                             coef=2, cutoff=NULL, pickCutoff=FALSE,
                             cutoffQ=0.99,
                             maxGap=500, smooth=FALSE,
                             smoothFunction=loessByCluster,
                             useWeights=FALSE, B=1000, verbose=TRUE, ...){
  if(!is.matrix(mat))
    stop("'mat' must be a matrix.")
  if(ncol(mat) != nrow(design))
    stop("Number of columns of data matrix must match number of rows of design matrix")
  if(is.null(cutoff) & !pickCutoff) stop("Must pick a cutoff")
  if(is.null(cutoff) & pickCutoff & B<1) stop("Must do at least one permution to pick a cutoff")
  if(!is.null(cutoff) & pickCutoff)
                             warning("cutoff was supplied so ignoring pickCutoff=TRUE")
  
  if (!getDoParRegistered())
    registerDoSEQ()
  workers <- getDoParWorkers()
  backend <- getDoParName()
  version <- getDoParVersion()
  subverbose <- max(as.integer(verbose) - 1L, 0)
  if (verbose) {
    if (workers == 1) {
      mes <- "bumphunterEngine: Using a single core (backend: %s, version: %s)"
      message(sprintf(mes, backend, version))
    } else {
      mes <- "bumphunterEngine: Parallelizing using %s workers/cores (backend: %s, version: %s)"
      message(sprintf(mes, workers, backend, version))
    }
  }
  
  ##B is the number of random samples we take
  if (is.null(chr))
    chr <- rep("Unspecified", length(pos))
  if(is.null(cluster))
    cluster <- clusterMaker(chr, pos, maxGap=maxGap)
  
  if(B>0 & B< workers) stop("B must be bigger than workers (or 0)")
  
  if(verbose) message("bumphunterEngine: Computing coefficients.")
  
  if(useWeights & smooth){
    tmp <- .getEstimate(mat=mat, design=design, coef=coef, full=TRUE)
    rawBeta <- tmp$coef
    weights <- tmp$sigma
    rm(tmp)
  } else{
    rawBeta <- .getEstimate(mat=mat, design=design, coef=coef, full=FALSE)
    weights <- NULL
  }
  
  if(smooth){
    if(verbose) message("bumphunterEngine: Smoothing coefficients.")
    beta <- smoother(y=rawBeta, x=pos, cluster=cluster, weights=weights,
                     smoothFunction=smoothFunction, verbose=subverbose,...) 
    Index <- which(beta$smoothed)
    beta <- beta$fitted
  } else {
    beta <- rawBeta
    Index <- seq(along=beta)
  }
  
  if(B>0){
    if (verbose) message("bumphunterEngine: Performing ", B, " permutations.")
  
  
    if(useWeights & smooth){
      tmp <- .getEstimate(mat, design, coef, B, full=TRUE)
      permRawBeta <- tmp$coef
      weights <- tmp$sigma
      rm(tmp)
    } else{
      permRawBeta <- .getEstimate(mat, design, coef, B, full=FALSE)
      weights <- NULL
    }
    
    ## Get individual p-values based on permutation of samples
    ## For each permutation we consider whether the absolute value of
    ## the observed beta is 
    ##
    if(verbose) message("bumphunterEngine: Computing marginal permutation p-values.")
    
    precision <- sqrt(.Machine$double.eps)
    sumGreaterOrEqual <- rowSums(abs(permRawBeta) >= abs(as.vector(rawBeta)) |
                                 abs(abs(permRawBeta)-abs(as.vector(rawBeta)))<=precision)
    pvs <- (sumGreaterOrEqual + 1L) / (B + 1L)
    
    if(smooth){
      if(verbose) message("bumphunterEngine: Smoothing permutation coefficients.")
      permBeta <- smoother(y=permRawBeta, x=pos, cluster=cluster, weights=weights,
                           smoothFunction=smoothFunction,
                           verbose=subverbose,...)$fitted
    } else permBeta <- permRawBeta
    
    if(is.null(cutoff))
      cutoff <- quantile(abs(permBeta), cutoffQ, na.rm=TRUE)
    if(verbose) message(sprintf("bumphunterEngine cutoff: %s", round(cutoff,3)))

  }
  ###Done with permutations

  if(verbose) message("bumphunterEngine: Finding regions.")
  tab <- regionFinder(x=beta, chr=chr, pos=pos, cluster=cluster,
                      cutoff=cutoff, ind=Index, verbose=FALSE)
  if (nrow(tab)==0) {
    if (verbose) message ("bumphunterEngine: No bumps found!")
    return(list(table=NA, coef=rawBeta, fitted=beta, pvaluesMarginal=NA))
  }
  
  if (B<1) {
    return(list(table=tab, coef=rawBeta, fitted=beta, pvaluesMarginal=NA))
  }
  
  if (verbose) message("bumphunterEngine: Computing regions for each permutation.")
  chunksize <- ceiling(B/workers)
  #nulltabs <- foreach(subMat=iter(permBeta, by="col", chunksize=chunksize),
                      #.combine="c", .packages = "bumphunter") %dorng% {
                        #apply(subMat, 2, regionFinder, chr=chr,
                              #pos=pos,
                              #cluster=cluster,
                              #cutoff=cutoff, ind=Index, verbose=FALSE)
                      #}
  # make binding of 'subMat' explicit
  foreach_object <- foreach(subMat=iter(permBeta, by="col", chunksize=chunksize),
	.combine="c", .packages = "bumphunter")
  foreach_argname <- foreach_object$argnames[1]
  nulltabs <- foreach_object %dorng% { 
	subMat <- get(foreach_argname)
	apply(subMat, 2, regionFinder, chr=chr, pos=pos, cluster=cluster,
		cutoff=cutoff, ind=Index, verbose=FALSE)
  }
  attributes(nulltabs)[["rng"]] <- NULL

  if (verbose) message("bumphunterEngine: Estimating p-values and FWER.")

  L <- V <- A <- as.list(rep(0, B))
  for(i in 1:B) {
    nulltab <- nulltabs[[i]]
    if (nrow(nulltab)>0) { 
      L[[i]] <- nulltab$L
      V[[i]] <- nulltab$value
      A[[i]] <- nulltab$area
    }
  }
  ###for observed length and height
  ###compute the total compute total number of times
  ###it is seen in permutations
  tots <- sapply(seq(along=V), function(i) {
    apply(cbind(tab$L,abs(tab$value)), 1, function(x) {
      sum(L[[i]]>=x[1] & abs(V[[i]])>=x[2])
    })
  })
  if (is.vector(tots)) {
    tots <- matrix(tots, nrow=1)
  }
  ##This is like a FWER
  rate1 <- rowMeans(tots>0)
  
  ##Now compute pvalues by assuming everything is exchangeable
  pvalues1 <- rowSums(tots)/sum(sapply(nulltabs,nrow))

  tots2 <- sapply(seq(along=A), function(i) {
    sapply(tab$area, function(x) { sum(A[[i]]>=x[1]) })
  })
  if (is.vector(tots2)) {
    tots2 <- matrix(tots2, nrow=1)
  }
  rate2 <- rowMeans(tots2>0)
  pvalues2 <- rowSums(tots2)/sum(sapply(nulltabs,nrow))
 
  tab$p.value <- pvalues1
  tab$fwer <- rate1

  tab$p.valueArea <- pvalues2
  tab$fwerArea <- rate2

  tab <- tab[order(tab$fwer,-tab$area),]
  return(list(table=tab, coef=rawBeta, fitted=beta, pvaluesMarginal=pvs, null=list(value=V,length=L)))
}

