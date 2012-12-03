test_bumphunter <- function () {
    ## Generate dummy data
    set.seed(123)
    dat <- dummyData()
    
    ## Test 1
    set.seed(123)
    bumps <- bumphunter (dat$mat, design=dat$design, chr=dat$chr, pos=dat$pos,
                         cluster=dat$cluster, coef=2, cutoffQ=0.9,
                         smooth=TRUE, B=500, verbose=TRUE, smoothFunction=loessByCluster)
    load(file.path(path.package("bumphunter"), "unitTests", "bumpsTest1.rda"))
    checkEquals(bumps, bumpsTest1)
    
    ## Test 2: chr unspecified
    set.seed(123)
    bumps <- bumphunter (dat$mat, design=dat$design, pos=dat$pos,
                         cluster=dat$cluster, coef=2, cutoffQ=0.9,
                         maxGap=500, smooth=TRUE, B=500, verbose=TRUE, smoothFunction=loessByCluster)
    load(file.path(path.package("bumphunter"), "unitTests", "bumpsTest2.rda"))
    checkEquals(bumps, bumpsTest2)
}

test_bumphunterParallel <- function() {
    set.seed(123)
    dat <- dummyData()

    if(require(doParallel)) {
        registerDoParallel(cores = 2)
        set.seed(123)
        bumps <- bumphunter (dat$mat, design=dat$design, chr=dat$chr, pos=dat$pos,
                             cluster=dat$cluster, coef=2, cutoffQ=0.9,
                             smooth=TRUE, B=500, verbose=TRUE, smoothFunction=loessByCluster)
        load(file.path(path.package("bumphunter"), "unitTests", "bumpsTest1.rda"))
        checkEquals(bumps, bumpsTest1)
        cl <- makeCluster(2)
        registerDoParallel(cl = cl)
        set.seed(123)
        bumps <- bumphunter (dat$mat, design=dat$design, chr=dat$chr, pos=dat$pos,
                             cluster=dat$cluster, coef=2, cutoffQ=0.9,
                             smooth=TRUE, B=500, verbose=TRUE, smoothFunction=loessByCluster)
        load(file.path(path.package("bumphunter"), "unitTests", "bumpsTest1.rda"))
        checkEquals(bumps, bumpsTest1)
    }
}
    


    
