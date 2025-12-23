##--------------------------------------------------------
## Comparison of ORA
##--------------------------------------------------------

##remotes::install_github("MaayanLab/genesetr")

source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")
source("../R/goat.R")
source("../R/fastFET.R")

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
fc <- scale(fc,center=FALSE)[,1]

gg <- intersect(rownames(G),names(fc))
G <- G[gg,]
fc <- fc[gg]
G <- G[,names(gmt)]
dim(G)

sig.genes <- names(which(abs(fc) > 1))
genes=sig.genes
gs <- gmt[[6]]
bg = names(fc)

manual.ft <- function(sig.genes, gmt, bg) {
  pv <- rep(NA,length(gmt))
  or <- rep(NA,length(gmt))  
  i=1
  for(i in 1:length(gmt)) {
    gs <- gmt[[i]]
    a <- sum(gs %in% sig.genes)
    b <- length(sig.genes) - a
    c <- length(gs) - a
    d <- length(bg) - (a+b+c)
    res <- fisher.test(matrix(c(a, c, b, d), 2, 2), alternative='greater')
    res
    pv[i] <- res$p.value
    or[i] <- res$estimate
  }
  names(pv) <- names(gmt)
  names(or) <- names(gmt)  
  cbind(pval=pv, odd.ratio=or)
}

system.time(res1 <- playbase::gset.fisher(sig.genes, gmt, fdr=1, background=names(fc)) )
system.time(res2 <- gset.fastFET(sig.genes, G, bg=names(fc)))
res1 <- res1[match(names(gmt),rownames(res1)),]
res2 <- res2[match(names(gmt),rownames(res2)),]
head(res1$p.value)
head(res2$p.value)

fres <- list()
qres <- list()
mres <- list()
sig.genes <- list()

z=1
for(z in c(1,2,3)) {
  k <- paste0("z=",z)
  sig <- names(which(abs(fc) > z))
  sig.genes[[k]] <- sig
  res <- playbase::gset.fisher(sig, gmt, fdr=1,
    background=names(fc), nmin=1, min.genes=1, max.genes=99999)
  res <- res[match(names(gmt),rownames(res)),]
  rownames(res) <- names(gmt)
  fres[[k]] <- res

  res <- quasiFET(sig, G)
  res <- res[match(names(gmt),res$pathway),]
  rownames(res) <- names(gmt)  
  qres[[k]] <- res

  res <- manual.ft(sig, gmt, bg=names(fc))
  res <- res[match(names(gmt),rownames(res)),]
  mres[[k]] <- res
}

sapply(sig.genes,length)

S1 <- do.call(cbind, lapply(fres, function(f) f$odd.ratio))
S2 <- do.call(cbind, lapply(qres, function(f) f$tau))
S3 <- do.call(cbind, lapply(mres, function(f) f[,"odd.ratio"]))
rownames(S1) <- colnames(G)
rownames(S2) <- colnames(G)
rownames(S3) <- colnames(G)
ss <- cbind(S1, S2, S3)
ss[is.na(ss)] <- 1
pairs(-log(ss))


P1 <- do.call(cbind, lapply(fres, function(f) f$p.value))
P2 <- do.call(cbind, lapply(qres, function(f) f$pval))
P3 <- do.call(cbind, lapply(mres, function(f) f[,"pval"]))
rownames(P1) <- colnames(G)
rownames(P2) <- colnames(G)
rownames(P3) <- colnames(G)
colnames(P1) <- paste0("gset.fisher\n",colnames(P1))
colnames(P2) <- paste0("quasi.fisher\n",colnames(P2))
colnames(P3) <- paste0("fisher.test\n",colnames(P3))

pp <- cbind( gFT=P1, qFT=P2, FT=P3)
pp[is.na(pp)] <- 1

pairs(pp)
pairs(-log10(pp))
pr <- apply(-log10(pp),2,rank)
pairs(pr)

cx1 <- ss[,1]
plot( pp[,1], pp[,4], cex=cx1)
plot(-log10(pp[,1]), -log10(pp[,4]), cex=0.8*cx1)
plot(ss[,1], -log10(pp[,4]), cex=cx1)

##-----------------------------------------------------
##-----------------------------------------------------
##-----------------------------------------------------
source("../R/fastFET.R")

G <- Matrix::t(playdata::GSETxGENE)
sig.genes <- head(names(sort(Matrix::rowMeans(G),decreasing=TRUE)),200)
gmt <- mat2gmt(G)




system.time(f1 <- playbase::gset.fisher(sig.genes, gmt, fdr=1, background=rownames(G)))
system.time(f2 <- quasiFET(sig.genes, G, bg=rownames(G)))
system.time(f3 <- gset.fastFET(sig.genes, G, bg=rownames(G)))

gs <- colnames(G)
f1 <- f1[match(gs,rownames(f1)),]
f2 <- f2[match(gs,rownames(f2)),]
f3 <- f3[match(gs,rownames(f3)),]

S <- cbind( gset.fisher=log(f1$odd.ratio), fastFET=log(f3$odd.ratio), quasiFET=f2$tau)
pairs(S)

P <- cbind( gset.fisher=f1$p.value, fastFET=f3$p.value, quasiFET=f2$pval)
cx1 <- log(1+Matrix::colSums(G!=0))**2
range(cx1)
pairs(-log10(P), cex=0.1*cx1, pch=20)




source("../R/fastFET.R")
G <- Matrix::t(playdata::GSETxGENE)
bg <- rownames(G)
genes <- head(names(sort(Matrix::rowMeans(G),dec=TRUE)),200)
dim(G)

system.time(pv1 <- gset.fastFET(genes, G, bg=bg))

gsize <- Matrix::colSums(G!=0)
genes <- intersect(genes, rownames(G))
a <- Matrix::colSums(G[genes,]!=0)
b <- length(genes) - a
c <- gsize - a
d <- length(bg) - (a+b+c) 
##pv <- genesetr::fastFET(a,b,c,d)
system.time(pv <- fastFET(a,b,c,d))
names(pv) <- colnames(G)
head(pv)

system.time(pv2 <- corporaFET(a,b,c,d))
names(pv2) <- colnames(G)
head(pv2)


odd.ratio <- (a/b)/(c/d)
qv <- p.adjust(pv, method="fdr")
overlap <- paste0(a,"/",gsize)
data.frame(p.value=pv, q.value=qv, odd.ratio=odd.ratio, overlap=overlap)
