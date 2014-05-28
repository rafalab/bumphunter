test_loessByCluster <- function () {
    set.seed(123)
    dat <- dummyData()
    smoothed <- loessByCluster(y=dat$mat[,1], cluster=dat$cluster, bpSpan = 1000,
                               minNum=7, minInSpan=5, maxSpan=1)
    load(file.path(path.package("bumphunter"), "unitTests/loessSmoothed.rda"))
    checkEquals(smoothed, loessSmoothed)
}

test_runmedByCluster <- function () {
    set.seed(123)
    dat <- dummyData()
    smoothed <- runmedByCluster(y=dat$mat[,1], cluster=dat$cluster, k=5, endrule="constant")
    load(file.path(path.package("bumphunter"), "unitTests/runmedSmoothed.rda"))
    checkEquals(smoothed, runmedSmoothed)
}

test_smoother <- function() {
    set.seed(123)
    dat <- dummyData()
    ## loessByCluster
    smoothed <- smoother(y=dat$mat[,1], cluster=dat$cluster,
                         smoothFunction=loessByCluster, 
                         bpSpan = 1000, minNum=7, minInSpan=5, maxSpan=1)
    load(file.path(path.package("bumphunter"), "unitTests/loessSmoothed.rda"))
    checkEquals(smoothed, loessSmoothed)
    ## runmedByCluster
    smoothed <- smoother(y=dat$mat[,1], cluster=dat$cluster,
                         smoothFunction=runmedByCluster, 
                         k=5, endrule="constant")
    load(file.path(path.package("bumphunter"), "unitTests/runmedSmoothed.rda"))
    checkEquals(smoothed, runmedSmoothed)  
}
