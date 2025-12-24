library(playbase)
library(ultragsea)
source("../R/plaidtest.R")
source("../R/gsetcor.R")

pgx <- playdata::GEIGER_PGX
X <- pgx$X
y <- pgx$samples$activated
G <- Matrix::t(playdata::GSETxGENE)
gg <- intersect(rownames(G),rownames(X))
X <- X[gg,]
G <- G[gg,]
G <- G[,grep("^GO",colnames(G))]
gmt <- mat2gmt(G)
length(gmt)
gsize <- Matrix::colSums(G!=0)

ref="notact"
fc <- rowMeans(X[,y!=ref],na.rm=TRUE) - rowMeans(X[,y==ref],na.rm=TRUE)
  
test1 <- plaid.ttest(X, G, y, ref, nsmooth=10)
test1 <- test1[colnames(G),]

test1a <- plaid.limma(X, G, y, ref, nsmooth=10)
test1a <- test1a[colnames(G),]

test2 <- plaid.cortest(X, G, y, ref)
test2 <- test2[colnames(G),]

test2a <- plaid.dualtest(X, G, y, ref, test1="ttest")
test2a <- test2a[colnames(G),]

test2b <- plaid.dualtest(X, G, y, ref, test1="cortest")
test2b <- test2b[colnames(G),]

cX <- X - rowMeans(X)
gsetR <- as.matrix(gset.cortest(G, cX)$rho)
test3 <- matrixTests::row_t_welch(gsetR[,which(y=="act"),drop=FALSE],
  gsetR[,which(y=="notact"),drop=FALSE] )
test3 <- test3[colnames(G),]

test4 <- gset.cortest(G, fc, compute.p = TRUE,
  use.rank = FALSE, corshrink = 10)
test4 <- sapply(test4, function(x) x[,1])
test4 <- test4[colnames(G),]

test5 <- goat(gmt, fc)
test5 <- test5[match(colnames(G),test5$pathway),]
rownames(test5) <- colnames(G)

P = cbind(
  plaidTT = test1$pval,
  plaidCOR = test2$pval,
  dualTT = test2a$pval,
  dualCOR = test2b$pval,
  gsrhoTT = test3$pvalue,
  gsetCOR = test4[,"p.value"],
  goat = test5[,"pval"]
)
dim(P)

cx1 <- log(1+Matrix::colSums(G!=0))
range(cx1)
pairs(-log10(P), cex=0.2*cx1)

##------------------------------------------------------------------
##------------------------------------------------------------------
##------------------------------------------------------------------
source("~/Playground/playbase/dev/include.R",chdir=TRUE)
library(playbase)

pp <- apply(P, 2, function(p) head(rownames(P)[order(p)],200))
pp <- as.list(data.frame(pp))
length(pp)

aa <- pgx$samples
gsetX <- plaid::plaid(X, G)
gx.heatmap(gsetX[pp[[1]],], col.annot=aa, mar=c(10,18))
gx.heatmap(gsetX[pp[[2]],], col.annot=aa, mar=c(10,18))
gx.heatmap(gsetX[pp[[3]],], col.annot=aa, mar=c(10,18))
gx.heatmap(gsetX[pp[[4]],], col.annot=aa, mar=c(10,18))
gx.heatmap(gsetX[pp[[5]],], col.annot=aa, mar=c(10,18))
gx.heatmap(gsetX[pp[[6]],], col.annot=aa, mar=c(10,18))

head(pp[[1]],10)
head(pp[[2]],10)
head(pp[[3]],20)
head(pp[[4]],20)
head(pp[[5]],20)

par(mfrow=c(6,5), mar=c(2,4,2,1))
for(k in 1:6) {
  for(i in 1:5) {
    m <- pp[[k]][[i]]
    gsea.enplot(fc, gmt[[m]], main=substring(m,1,40))
  }
}

