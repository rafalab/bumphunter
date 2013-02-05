library(bumphunter)
  
## Generate dummy data
set.seed(123)
dat <- dummyData()

## loessByCluster      
loessSmoothed <- loessByCluster(cluster=dat$cluster, y=dat$mat[,1],
                              bpSpan = 1000, minNum=7, minInSpan=5, maxSpan=1)
save(loessSmoothed, file="inst/unitTests/loessSmoothed.rda")

## runmedByCluster
runmedSmoothed <- runmedByCluster(cluster=dat$cluster, y=dat$mat[,1],  k=5, endrule="constant")
save(runmedSmoothed, file="inst/unitTests/runmedSmoothed.rda")

## bumphunter test 1
set.seed(123)
bumpsTest1 <- bumphunter (dat$mat, design=dat$design, chr=dat$chr, pos=dat$pos, cluster=dat$cluster,
                          coef=2, cutoffQ=0.9, smooth=TRUE, B=500,
                          verbose=TRUE, smoothFunction=loessByCluster)
save(bumpsTest1, file="inst/unitTests/bumpsTest1.rda")

## bumphunter test 2: chr unspecified
set.seed(123)
bumpsTest2 <- bumphunter (dat$mat, design=dat$design, pos=dat$pos, cluster=dat$cluster,
                          coef=2, cutoffQ=0.9, maxGap=500, smooth=TRUE, B=500,
                          verbose=TRUE, smoothFunction=loessByCluster)
save(bumpsTest2, file="inst/unitTests/bumpsTest2.rda")


