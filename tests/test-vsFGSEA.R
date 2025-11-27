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

f0 <- fgsea::fgsea(gmt, fc)
f0 <- fgsea::fgsea(gmt, fc, eps=0)
f1 <- fgsea::fgseaSimple(gmt, fc, nperm=10000)
f1 <- fgsea::fgseaSimple(gmt, fc, nperm=20)
f2 <- fgsea::fgseaMultilevel(gmt, fc, nPermSimple=100, eps=1)

f0 <- data.frame(f0, row.names=f0$pathway)
head(f0)
head(f1)
head(f2)
S <- cbind(f0$NES, f1$NES, f2$NES)
pairs(S)

z1 <- ultragsea(fc, G, format='simple')
z2 <- ultragsea(fc, G, format='as.gsea')
z3 <- ultragsea(fc, G, format='as.gsea2')
head(z1)
S <- cbind(z1$score, z2$NES, z3$NES)
pairs(S)

range(z1$ES)
sd(z1$NES)

source("../R/goat.R")
g1 <- goat(gmt, fc, method="goat", filter=FALSE)
g2 <- goat(gmt, fc, method="goat_bootstrap", filter=FALSE)
g3 <- goat(gmt, fc, method="gsea", filter=FALSE)
rownames(g1) <- g1$pathway
rownames(g2) <- g2$pathway
head(g1)
g2 <- g2[match(g1$pathway,g2$pathway),]
g3 <- g3[match(g1$pathway,g3$pathway),]

ss <- g1$pathway
S <- cbind(
  fgsea = f0[ss,]$NES,
  ultragsea = z1[ss,]$score,
  goat = g1$score,
  goat.bs = g2$score,
  goat.gsea = g3$score
)
pairs(S)

P <- cbind( fgsea = f0[ss,]$pval,
  ultragsea = z1[ss,]$pval,
  goat = g1$pval, goat.bs = g2$pval,
  goat.gsea = g3$pval
)
pairs(-log10(P))


##--------------------------------------------------------
##--------------------------------------------------------
##--------------------------------------------------------

f0 <- fgsea::fgsea(gmt, fc, eps=0, gseaParam=0)
f0 <- f0[match(names(gmt),f0$pathway),]

f1 <- fgsea::fgsea(gmt, fc, eps=0, gseaParam=1)
f1 <- f1[match(names(gmt),f1$pathway),]

z0 <- ultragsea(fc, G, method="ztest", format='as.gsea')
z0 <- z0[match(names(gmt),z0$pathway),]
z1 <- ultragsea(fc, G, method="cor", format='as.gsea')
z1 <- z1[match(names(gmt),z1$pathway),]

g1 <- goat(gmt, fc)
g2 <- goat(gmt, fc, method="goat_bootstrap")
g1 <- g1[match(names(gmt),g1$pathway),]
g2 <- g2[match(names(gmt),g2$pathway),]

ungap <- function(x) ifelse(x<0, x+1, x-1)
S <- cbind(
  ES0=f0$ES, NES0=f0$NES, 
  ES1=f1$ES, NES1=f1$NES,
  ungapped.NES0=ungap(f0$NES),
  ungapped.NES1=ungap(f1$NES),
  ultragsea.ztest=z0$NES,
  ultragsea.cor=z1$NES, 
  goat=g1$score,
  goat.bs=g2$score
)
gsize <- Matrix::colSums(G!=0)[names(gmt)]
pairs(S, pch=20, cex = 0.2*log10(1+gsize)**2 )


cor.test( ungap(f0$NES), g1$score )
cor.test( ungap(f1$NES), z0$NES )
cor.test( g1$score, z0$NES )

P <- cbind(
  fgsea0 = f0$pval,fgsea1 = f1$pval,
  ultragsea.ztest = z0$pval, ultragsea.cor = z1$pval, 
  goat = g1$pval, goat.bs = g2$pval)
pairs(-log(P), pch=20, cex=1)

##--------------------------------------------------------
## benchmarking for runtime
##--------------------------------------------------------

library(peakRAM)
gmtx <- rep(gmt,150)
length(gmtx)
names(gmtx) <- make.unique(names(gmtx))
Gx <- do.call(cbind, rep(list(G),150))
colnames(Gx) <- make.unique(colnames(Gx))
dim(Gx)
system.time(G2<-gmt2mat(gmtx))

tt <- list()
nn <- c(100,1000,10000,50000,100000)
nt <- length(nn)
for(i in 1:10) {
  n=100
  for(n in nn) {
    t0 <- peakRAM(
      f0 <- fgsea::fgsea(gmtx[1:n], fc),
      g1 <- goat(gmtx[1:n], fc, filter=FALSE),
      f1 <- ultragsea(fc, Gx[,1:n], method="cor"),
      c1 <- gset.cor(fc, Gx[,1:n], compute.p=TRUE)
    )
    t0[,1] <- paste0(c("fgsea","goat","ultragsea","gsetcor"),".n",n)
    tt <- c(tt, list(t0))
  }
}
ttx <- do.call(rbind,tt)
ttx

tmean <- tapply(ttx[,2], ttx[,1], median)
tmean
acc = round(rep(tmean[1:nt], 4) / tmean, 2)
acc = matrix(acc, nrow=nt)[,2:4]
colnames(acc) <- c("goat","gsetcor","ultragsea")

tmat <- matrix(tmean, nrow=nt)
colnames(tmat) <- c("fgsea","goat","gsetcor","ultragsea")
B <- data.frame(time=tmat, speedup = acc)
rownames(B) <- paste0("N=",nn)
B
write.csv(B, file="timings.csv")



##----------------------------------------------
##----------------------------------------------
##----------------------------------------------
