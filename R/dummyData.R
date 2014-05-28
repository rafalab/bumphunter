## Generate dummy data that can be used for unit tests and simple examples

dummyData <- function(n1=5, n2=5, sd=0.2, l=100, spacing=100, clusterSpacing=1e5, numClusters=5) {
    ## removed a set seed
    cluster <- rep(1:numClusters, each=l/numClusters)  
    chr <- c(rep(1, l*.8), rep(2, l*.2))
    pos <- cluster*clusterSpacing + seq(0, by=spacing, length=l) + sample(l/10, l, replace=TRUE)
    ## Prepare data
    x <- c(rep(0, 2*l/numClusters), seq(0, 1, length.out=l/numClusters) * 2 * pi, rep(0, 2*l/numClusters))
    grp1 <- replicate(n1, sin(x)) + rnorm(n1*l, sd=sd)
    grp2 <- replicate(n2, rep(0,l)) + rnorm(n1*l, sd=sd)
    mat <- cbind(grp1, grp2)
    dd <- data.frame(x=factor(c(rep(1, n1), rep(2,n2))))
    design <- model.matrix(~ x, dd, contrasts = list(x="contr.helmert"))
    return(list(mat=mat, design=design, chr=chr, pos=pos, cluster=cluster, n1=n1, n2=n2))
}

