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
test1 <- test1[,c("mean.x","mean.y","mean.diff","pvalue")]
test1 <- test1[colnames(G),]

test2 <- gset.cor(G, fc, compute.p = TRUE,
  use.rank = FALSE, corshrink = 10)
test2 <- sapply(test2, function(x) x[,1])
test2 <- test2[colnames(G),]

test3 <- goat(gmt, fc)
test3 <- test3[match(colnames(G),test3$pathway),]
rownames(test3) <- colnames(G)

cX <- X - rowMeans(X)
gsetR <- as.matrix(gset.cor(cX, G)$rho)
test4 <- matrixTests::row_t_welch(gsetR[,which(y=="act"),drop=FALSE],
  gsetR[,which(y=="notact"),drop=FALSE] )
test4 <- test4[colnames(G),]

P = cbind(
  plaidTT = test1$pvalue,
  gsrhoTT = test4$pvalue,
  gsetCOR = test2[,"p.value"],
  goat = test3[,"pval"]
)
cx1 <- log(1+Matrix::colSums(G!=0))
range(cx1)
pairs(-log10(P), cex=0.3*cx1)

##------------------------------------------------------------------
##------------------------------------------------------------------
##------------------------------------------------------------------
source("~/Playground/playbase/dev/include.R",chdir=TRUE)

tt1 <- rownames(test1)[order(test1$pvalue)]
tt4 <- rownames(test4)[order(test4$pvalue)]
tt2 <- rownames(test2)[order(test2[,"p.value"])]
tt3 <- test3$pathway[order(test3[,"pval"])]

pp <- list(plaidTT=tt1, gsrhoTT=tt4, gsetCOR=tt2, goat=tt3)
pp <- lapply(pp, head, 200)
aa <- pgx$samples

gsetX <- plaid::plaid(X, G)
gx.heatmap(gsetX[pp[[1]],], col.annot=aa, mar=c(10,18))
gx.heatmap(gsetX[pp[[2]],], col.annot=aa, mar=c(10,18))
gx.heatmap(gsetX[pp[[3]],], col.annot=aa, mar=c(10,18))
gx.heatmap(gsetX[pp[[4]],], col.annot=aa, mar=c(10,18))

head(pp[[1]],10)
head(pp[[2]],10)
head(pp[[3]],20)
head(pp[[4]],20)

par(mfrow=c(4,4), mar=c(2,4,2,1))
for(k in 1:4) {
  for(i in 1:4) {
    m <- pp[[k]][[i]]
    gsea.enplot(fc, gmt[[m]], main=substring(m,1,40))
  }
}

k=2
k=3
k=1
head(test1[pp[[k]],])
par(mfrow=c(4,4), mar=c(2,4,2,1))
for(i in 1:16) {
  m <- pp[[k]][[i]]
  gsea.enplot(fc, gmt[[m]], main=substring(m,1,40) )
}
