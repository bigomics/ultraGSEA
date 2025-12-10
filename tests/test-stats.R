library(playbase)
library(limma)

source("~/Playground/playbase/dev/include.R",chdir=TRUE)
source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/goat.R")
source("../R/utils.R")

library(fgsea)
library(goat)
library(ultragsea)
data(examplePathways)
data(exampleRanks)
fc = exampleRanks
gmt = examplePathways
range(sapply(gmt,length))
gmt <- lapply(gmt, function(s) intersect(s,names(fc)))
gmt <- gmt[sapply(gmt,length)>=10]
G <- gmt2mat(gmt)

gg <- intersect(rownames(G),names(fc))
G <- G[gg,]
fc <- fc[gg]

length(fc)
length(gmt)
range(fc)
range(sapply(gmt,length))

plot(sort(fc))
fc1 = scale(fc, center=FALSE)
fc1 = fc1 - mean(fc1)
plot(sort(fc1)); abline(h=0)

range(fc1)
table(sign(fc))
table(sign(fc1))
F <- cbind(fc1, fc1+0.5)
colMeans(F)

source("../R/ultragsea.R")
res1  <- gset.cor(G, F, compute.p=TRUE)
res2a <- gset.ztest(G, F, alpha=0)
res2b <- gset.ztest(G, F, alpha=1)
res3a <- gset.ztest(G, F, pdist="t", alpha=0)
res3b <- gset.ztest(G, F, pdist="t", alpha=1)
res4a <- gset.ztest(G, F, center=FALSE, alpha=0)
res4b <- gset.ztest(G, F, center=FALSE, alpha=1)
res5a <- gset.ztest(G, F, pdist="t", center=FALSE, alpha=0)
res5b <- gset.ztest(G, F, pdist="t", center=FALSE, alpha=1)

cx <- log10(1+Matrix::colSums(G!=0))
range(cx)

k=2
k=1
S <- cbind(
  res1$rho[,k], res2a$stat[,k], res2b$stat[,k],
  res3a$stat[,k], res3b$stat[,k],
  res4a$stat[,k], res4b$stat[,k],
  res5a$stat[,k], res5b$stat[,k])
dim(S)
colnames(S) <- c(
  "cor","ztest.a0","ztest.a1",
  "ttest.a0", "ttest.a1",
  "ztest0.a0","ztest0.a1",
  "ttest0.a0","ttest0.a1")
pairs(S, panel=function(x,y) {
  abline(h=0,v=0,lty=2,lwd=0.7)
  points(x,y,cex=0.3*cx,pch=19)
})

P <- cbind( res1$p.value[,k], res2a$p[,k], res2b$p[,k], 
  res3a$p[,k], res3b$p[,k], res4a$p[,k], res4b$p[,k],
  res5a$p[,k], res5b$p[,k])
dim(P)
colnames(P) <- c("cor","ztest.a0","ztest.a1",
  "ttest.a0","ttest.a1",
  "ztest0.a0","ztest0.a1",
  "ttest0.a0","ttest0.a1")
pairs(-log10(P), cex=0.3*cx,pch=19)

system.time(pnorm(rnorm(1e6)))
system.time(pt(rnorm(1e6),df=5))
system.time(pt(rnorm(1e6),df=100))

par(mfrow=c(4,4),mar=c(4,4,2,2))
kk <- round(seq(3,200,length.out=16))
tstat <- rnorm(1e3,sd=1) 
for(k in kk) {
  p1 = -log10(pnorm(tstat))
  p2 = -log10(pt(tstat,df=k))
  xl=yl= c(0, max(c(p1,p2)))
  plot( p1, p2, xlim=xl, ylim=yl,
    main=paste("K=",k))
}

## test matrix var
dim(G)
gset_size <- Matrix::colSums(G!=0)
gset_var <- apply(F, 2, function(fc) {
  gfc <- (G[gg,]!=0) * fc[gg]
  sparseMatrixStats::colVars(gfc) * length(gg) / gset_size
})

fsumsq <- Matrix::t(G!=0) %*% F**2 / gset_size
dim(fsumsq)
head(fsumsq)
fmeansq <- ((Matrix::t(G[gg,]!=0) %*% F[gg,,drop=FALSE]) / gset_size)**2
head(fmeansq)
gsvar <- (fsumsq - fmeansq) * ( gset_size / (gset_size - 1))
dim(gset_var)
dim(gsvar)
head(gset_var)
head(gsvar)

## manual
sum(fc==0)
fvar <- apply(gfc, 2, function(x) var(x[which(x!=0)]))
head(fvar,3)

cx <- log10(1+gset_size)**2
V <- cbind(gset_var[,1], gsvar[,1], fvar)
pairs(V, pch=20, cex=0.2*cx)
