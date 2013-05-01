test_bumphunter <- function () {
    ## Generate dummy data
    set.seed(123)
    dat <- dummyData()
    
    ## Test 1
    set.seed(123)
    bumps <- bumphunter (dat$mat, design=dat$design, chr=dat$chr, pos=dat$pos,
                         cluster=dat$cluster, coef=2, cutoff=0.28,
                         smooth=TRUE, B=500, verbose=TRUE, smoothFunction=loessByCluster)
    load(file.path(path.package("bumphunter"), "unitTests", "bumpsTest1.rda"))
    nam <- setdiff(names(bumps), "algorithm")
    checkEquals(bumps[nam], bumpsTest1[nam])
    
    ## Test 2: chr unspecified
    set.seed(123)
    bumps <- bumphunter (dat$mat, design=dat$design, pos=dat$pos,
                         cluster=dat$cluster, coef=2, cutoff=0.28,
                         maxGap=500, smooth=TRUE, B=500, verbose=TRUE, smoothFunction=loessByCluster)
    load(file.path(path.package("bumphunter"), "unitTests", "bumpsTest2.rda"))
    nam <- setdiff(names(bumps), "algorithm")
    checkEquals(bumps[nam], bumpsTest2[nam])
}

test_bumphunterParallel <- function() {
    set.seed(123)
    dat <- dummyData()

    if(require(doParallel)) {
        registerDoParallel(cores = 2)
        set.seed(123)
        bumps <- bumphunter (dat$mat, design=dat$design, chr=dat$chr, pos=dat$pos,
                             cluster=dat$cluster, coef=2, cutoff=0.28,
                             smooth=TRUE, B=500, verbose=TRUE, smoothFunction=loessByCluster)
        load(file.path(path.package("bumphunter"), "unitTests", "bumpsTest1.rda"))
        nam <- setdiff(names(bumps), "algorithm")
        checkEquals(bumps[nam], bumpsTest1[nam])
        cl <- makeCluster(2)
        registerDoParallel(cl = cl)
        set.seed(123)
        bumps <- bumphunter (dat$mat, design=dat$design, chr=dat$chr, pos=dat$pos,
                             cluster=dat$cluster, coef=2, cutoff=0.28,
                             smooth=TRUE, B=500, verbose=TRUE, smoothFunction=loessByCluster)
        load(file.path(path.package("bumphunter"), "unitTests", "bumpsTest1.rda"))
        nam <- setdiff(names(bumps), "algorithm")
        checkEquals(bumps[nam], bumpsTest1[nam])
    }
}
    


    
