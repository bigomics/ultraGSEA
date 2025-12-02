##--------------------------------------------------------
## benchmarking for runtime
##--------------------------------------------------------

source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")
source("../R/goat.R")

library(fgsea)
library(goat)

data(examplePathways)
data(exampleRanks)
fc = exampleRanks
length(fc)
gmt = examplePathways
range(sapply(gmt,length))
gmt <- lapply(gmt, function(s) intersect(s,names(fc)))
gmt <- gmt[sapply(gmt,length)>10]
G <- gmt2mat(gmt)

library(peakRAM)
gmtx <- rep(gmt,150)
length(gmtx)
names(gmtx) <- make.unique(names(gmtx))
Gx <- do.call(cbind, rep(list(G),150))
colnames(Gx) <- make.unique(colnames(Gx))
dim(Gx)

tt <- list()
nn <- c(100,1000)
nn <- c(100,1000,10000)
nn <- c(100,1000,10000,50000,100000)
nt <- length(nn)
for(i in 1:1) {
  n=100
  for(n in nn) {
    t0 <- peakRAM::peakRAM(
      f0 <- fgsea::fgsea(gmtx[1:n], fc),
      g1 <- goat(gmtx[1:n], fc, filter=FALSE),
      f1 <- ultragsea(fc, Gx[,1:n], method="cor"),
      c1 <- gset.cor(fc, Gx[,1:n], compute.p=TRUE),
      c2 <- gset.cor(fc, Gx[,1:n], compute.p=TRUE, zbias=10),
      z1 <- fc_ztest(fc, Gx[,1:n], zmat=FALSE),
      m1 <- limma::cameraPR(fc, gmt, use.ranks=FALSE)
    )
    t0[,1] <- paste0(c("fgsea","goat","ultragsea","gsetcor","gsetcorB",
      "ztest","cameraPR"),".n",n)
    tt <- c(tt, list(t0))
  }
}
ttx <- do.call(rbind,tt)
ttx

tmean <- tapply(ttx[,2], ttx[,1], median)
tmean

ii <- grep("fgsea",names(tmean))
nm <- length(tmean) / length(ii)
acc = round(rep(tmean[ii], nm) / tmean, 2)
acc = matrix(acc, nrow=nt)
colnames(acc) <- unique(sub("[.]n.*","",names(tmean)))

tmat <- matrix(tmean, nrow=nt)
colnames(tmat) <- unique(sub("[.]n.*","",names(tmean)))
B <- data.frame(time=tmat, speedup = acc)
rownames(B) <- paste0("N=",nn)
t(B)
write.csv(t(B), file="timings.csv")



##----------------------------------------------
##----------------------------------------------
##----------------------------------------------
