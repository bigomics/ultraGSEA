source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")

library(fgsea)
library(goat)

data(examplePathways)
data(exampleRanks)
fc = exampleRanks
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
g1 <- goat_test(fc, gmt, method="goat")
g2 <- goat_test(fc, gmt, method="gsea")
rownames(g1) <- g1$pathway
rownames(g2) <- g2$pathway
head(g1)
g2 <- g2[match(g1$pathway,g2$pathway),]

ss <- g1$pathway
S <- cbind(
  fgsea = f0[ss,]$NES,
  ultragsea = z1[ss,]$score,
  goat = g1$score, goat.gsea = g2$score)
pairs(S)

P <- cbind( fgsea = f0[ss,]$pval,
  ultragsea = z1[ss,]$pval,
  goat = g1$pval)
pairs(-log10(P))


##--------------------------------------------------------
##--------------------------------------------------------
##--------------------------------------------------------

f0 <- fgsea::fgsea(gmt, fc, eps=0, gseaParam=0)
f0 <- f0[match(names(gmt),f0$pathway),]

f1 <- fgsea::fgsea(gmt, fc, eps=0, gseaParam=1)
f1 <- f1[match(names(gmt),f1$pathway),]

z0 <- ultragsea(fc, G, format='as.gsea')
z0$NES <- z0$NES / sd(z0$NES,na.rm=TRUE) * 1.0
z0 <- z0[match(names(gmt),z0$pathway),]

g1 <- goat_test(fc, gmt)
g2 <- goat_test(fc, gmt, method="gsea")
g1 <- g1[match(names(gmt),g1$pathway),]
g2 <- g2[match(names(gmt),g2$pathway),]

ungap <- function(x) ifelse(x<0, x+1, x-1)
S <- cbind(
  ES0=f0$ES, NES0=f0$NES, 
  ES1=f1$ES, NES1=f1$NES,
  ungapped.NES0=ungap(f0$NES),
  ungapped.NES1=ungap(f1$NES),
  ultragsea=z0$NES, 
  goat=g1$score,
  goat.gsea=g2$score
)
gsize <- Matrix::colSums(G!=0)[names(gmt)]
pairs(S, pch=20, cex = 0.2*log10(1+gsize)**2 )

P <- cbind( fgsea0 = f0$pval,fgsea1 = f1$pval,
           ultragsea = z0$pval, 
           goat = g1$pval, goat.rnk = g2$pval)
pairs(-log(P), pch=20, cex=1)

##--------------------------------------------------------
## benchmarking for runtime
##--------------------------------------------------------

library(peakRAM)
gmtx <- rep(gmt,100)
length(gmtx)
names(gmtx) <- make.unique(names(gmtx))
Gx <- do.call(cbind, rep(list(G),100))
colnames(Gx) <- make.unique(colnames(Gx))
dim(Gx)
system.time(G2<-gmt2mat(gmtx))

tt <- list()
for(i in 1:3) {
    n=100
    for(n in c(100,1000,10000)) {
        t0 <- peakRAM(
            f0 <- fgsea::fgsea(gmtx[1:n], fc),
            f1 <- ultragsea(fc, Gx[,1:n], method="cor"),
            g1 <- goat_test(fc, gmtx[1:n])
        )
        t0[,1] <- paste0(c("fgsea","ultragsea","goat"),".n",n)
        tt <- c(tt, list(t0))
    }
}
ttx <- do.call(rbind,tt)
ttx
tapply(ttx[,4], ttx[,1], mean)

tmean <- tapply(ttx[,2], ttx[,1], mean)
acc = round(rep(tmean[1:3],3) / tmean, 2)
acc = matrix(acc, nrow=3)[,2:3]
colnames(acc) <- c("goat","ultragsea")
B <- data.frame(
  fgsea = tmean[1:3],
  goat = tmean[4:6],
  ultragsea = tmean[7:9],
  speedup = acc
)
rownames(B) <- paste0("N=",c(100,1000,10000))
B
write.csv(B, file="timings.csv")



##----------------------------------------------
##----------------------------------------------
##----------------------------------------------
