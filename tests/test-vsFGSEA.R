source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")

library(fgsea)
data(examplePathways)
data(exampleRanks)
fc = exampleRanks
gmt = examplePathways
G <- gmt2mat(gmt)

f0 <- fgsea::fgsea(gmt, fc)
f0 <- fgsea::fgsea(gmt, fc, eps=0)
f1 <- fgsea::fgseaSimple(gmt, fc, nperm=10000)
f1 <- fgsea::fgseaSimple(gmt, fc, nperm=20)
f2 <- fgsea::fgseaMultilevel(gmt, fc, nPermSimple=100, eps=1)

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

z1 <- ultragsea(fc, G, format='as.gsea2')
z1$NES <- z1$NES / sd(z1$NES,na.rm=TRUE) * 1.0
z1 <- z1[match(names(gmt),z1$pathway),]

c1 <- limma::cameraPR(fc, gmt, use.ranks=FALSE)
c1 <- c1[names(gmt),]

ungap <- function(x) ifelse(x<0, x+1, x-1)
gsize <- Matrix::colSums(G!=0)[names(ks)]
S <- cbind(
  ES0=f0$ES, NES0=f0$NES, 
  ES1=f1$ES, NES1=f1$NES,
  ungapped.NES1=ungap(f1$NES),
  z0=z0$NES, z1=z1$NES
)
pairs(S, pch=20, cex = 0.2*log10(1+gsize)**2 )

P <- cbind( fgsea0=f0$pval,fgsea1=f1$pval,
           z0=z0$pval, z1=z1$pval,
           camera=c1$PValue)
pairs(-log(P), pch=20, cex=1)


## benchmarking for runtime
library(peakRAM)
gmtx <- rep(gmt,80)
names(gmtx) <- make.unique(names(gmtx))
Gx <- do.call(cbind, rep(list(G),80))
colnames(Gx) <- make.unique(colnames(Gx))
dim(Gx)
system.time(G2<-gmt2mat(gmtx))



tt <- list()
for(i in 1:10) {
    n=100
    for(n in c(100,1000,10000,100000)) {
        t0 <- peakRAM(
            f0 <- fgsea::fgsea(gmtx[1:n], fc),
            f1 <- ultragsea(fc, Gx[,1:n], method="cor")
        )
        t0[,1] <- paste0(c("fgsea","ultragsea"),".n",n)
        tt <- c(tt, list(t0))
    }
}
ttx <- do.call(rbind,tt)
ttx
tapply(ttx[,4], ttx[,1], mean)

tmean <- tapply(ttx[,2], ttx[,1], mean)
acc = round(tmean[1:4] / tmean[5:8],2)
B <- data.frame( fgsea=tmean[1:4], ultragsea=tmean[5:8], speedup=acc)
rownames(B) <- paste0("N=",c(100,1000,10000,100000))
B
write.csv(B, file="timings-vs-fgsea2.csv")



##----------------------------------------------
##----------------------------------------------
##----------------------------------------------
