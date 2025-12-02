library(playbase)
library(limma)

source("~/Playground/playbase/dev/include.R",chdir=TRUE)
source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/goat.R")

counts  <- playbase::COUNTS
samples <- playbase::SAMPLES
head(samples)
X <- playbase::logCPM(counts)
y <- samples$activated

pgx <- playbase::pgx.load("~/Playground/pgx/tcga-brca_pub.pgx")
X <- pgx$X
head(pgx$samples)
y <- 1*(pgx$samples$ER_STATUS=="Positive")
#y <- 1*(pgx$samples$.cell_cycle=="G1")

## compute logFC
res <- playbase::gx.limma(X, pheno=y, fdr=1, lfc=0)
fc <- res$logFC
names(fc) <- rownames(res)
head(res)

G <- Matrix::t(playdata::GSETxGENE)
G <- G[,grep("HALLMARK",colnames(G))]
G <- G[,grep("^GO_",colnames(G))]
gg <- intersect(rownames(X),rownames(G))
G <- G[gg,]
X <- X[gg,]
fc <- fc[gg]

G <- G[, Matrix::colSums(G!=0)>=10]
gmt <- mat2gmt(G)
G <- G[gg,names(gmt)]
dim(G)

peakRAM::peakRAM(G<-gmt2mat(gmt))

res <- list()
tt <- peakRAM::peakRAM(
  res$fgsea <- fgsea::fgsea(gmt, fc, eps=0),
  res$ultragsea.z <- ultragsea(fc, G, method='ztest'),
  res$ultragsea.t <- ultragsea(fc, G, method='ttest'),
  res$ultragsea.c <- ultragsea(fc, G, method='cor'),
  res$cor <- gset.cor(fc, G, compute.p=TRUE, use.rank=FALSE),
  res$rankcor  <- gset.cor(fc, G, compute.p=TRUE, use.rank=TRUE),
  res$goatPC <- goat(gmt, fc, filter=FALSE, method="goat_precomputed"),
  res$goatBS <- goat(gmt, fc, filter=FALSE, method="goat_bootstrap"),
  res$cameraPR <- limma::cameraPR(fc, gmt, use.ranks=FALSE)
)
tt
tt[,1] <- names(res)
tt

par(mfrow=c(2,2), mar=c(4,10,4,2))
rtime <- tt[,2]; names(rtime)=tt[,1]
barplot(sort(rtime, dec=TRUE), horiz=TRUE,
  xlab="runtime (sec)", las=1)

## align results
gs <- names(gmt)
res$cameraPR <- res$cameraPR[gs,]
ii <- grep("gsea|goat",names(res))
for(i in ii) {
  res[[i]] <- res[[i]][match(gs,res[[i]]$pathway),]
}

res$camera$score <- -log10(res$camera$PValue) * (-1 + 2*(res$camera$Direction=="Up"))

## show test-statistics (e.g. NES,rho,t)
ii <- grep("ultragsea",names(res))
Z <- cbind( fgsea=res$fgsea$NES)
Z <- cbind(Z, do.call(cbind, lapply(res[ii], function(x) x$score)))
Z <- cbind(Z, cor=res$cor$rho[gs,1], rankcor=res$rankcor$rho[gs,1] )
Z <- cbind(Z, camera=res$camera$score)
head(Z)
pairs(Z, cex=3, pch='.')

x11()
P <- cbind( fgsea=res$fgsea$pval)
P <- cbind(P, do.call(cbind, lapply(res[ii], function(x) x$pval)))
P <- cbind(P, cor=res$cor$p.value[gs,1], rankcor=res$rankcor$p.value[gs,1] )
P <- cbind(P, camera=res$camera[gs,"PValue"] )

pairs(-log10(P), cex=0.6, pch=20)
pairs(-log10(P), cex=3, pch='.')

## volcano plots
mm <- intersect(colnames(Z),colnames(P))
par(mfrow=c(3,3))
for(m in mm) {
  plot( Z[,m], -log10(P[,m]), pch=20, cex=0.8, main=m)
}

##---------------------------------------------------------------
##----------------------- zGSEA ---------------------------------
##---------------------------------------------------------------
source("../R/ultragsea.R")

fres <- fgsea::fgsea(gmt, fc)
zres <- ultragsea(fc, G, method='cor')
zres <- ultragsea(fc, G, method='ztest')
head(zres)

gs <- names(gmt)
fres <- fres[match(gs,fres$pathway),]
zres <- zres[match(gs,zres$pathway),]

k=10
gg1 <- fres$leadingEdge[[k]]
gg2 <- zres$leadingEdge[[k]]
gg1
gg2
gg1 %in% gg2
gg2 %in% gg1

