library(playbase)
source("~/Playground/playbase/dev/include.R",chdir=TRUE)
source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")

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
gmt <- mat2gmt(G)
G <- G[gg,names(gmt)]
dim(G)

peakRAM::peakRAM(G<-gmt2mat(gmt))

rc <- playbase::signedRank(fc)

res <- list()
tt <- peakRAM::peakRAM(
  res$fgsea <- fgsea::fgsea(gmt, fc, eps=0),
  res$ultragsea.z <- ultragsea(fc, G, method='ztest'),
  res$ultragsea.t <- ultragsea(fc, G, method='ttest'),
  res$ultragsea.c <- ultragsea(fc, G, method='cor'),
  res$ultragsea.rz <- ultragsea(rc, G, method='ztest'),
  res$ultragsea.rt <- ultragsea(rc, G, method='ttest'),
  res$ultragsea.rc <- ultragsea(rc, G, method='cor'),  
  res$cor <- gset.cor(fc, G, compute.p=TRUE, use.rank=FALSE),
  res$rankcor  <- gset.cor(fc, G, compute.p=TRUE, use.rank=TRUE),
  res$camera <- cameraPR(fc, gmt, use.ranks=FALSE),
  res$camera.rnk <- cameraPR(fc, gmt, use.ranks=TRUE)  
)
tt
tt[,1] <- names(res)
tt

## align results
gs <- names(gmt)
res$fgsea <- res$fgsea[match(gs,res$fgsea$pathway),]
res$camera <- res$camera[gs,]
res$camera.rnk <- res$camera.rnk[gs,]
ii <- grep("gsea",names(res))
for(i in ii) {
  res[[i]] <- res[[i]][match(gs,res[[i]]$pathway),]
}

## show test-statistics (e.g. NES,rho,t)
ii <- grep("ultragsea",names(res))
z2 <- ifelse(res$fgsea$NES < 0, res$fgsea$NES + 1, res$fgsea$NES - 1)
Z <- cbind( fgsea=res$fgsea$NES, fgsea2=z2)
Z <- cbind(Z, do.call(cbind, lapply(res[ii], function(x) x$score)))
Z <- cbind(Z, cor=res$cor$rho[gs,1], rankcor=res$rankcor$rho[gs,1] )
Z <- cbind(Z, camera=1-res$camera[gs,"PValue"], camera.rnk=1-res$camera[gs,"PValue"] ) 
head(Z)
pairs(Z, cex=3, pch='.')

x11()
P <- cbind( fgsea=res$fgsea$pval)
P <- cbind(P, do.call(cbind, lapply(res[ii], function(x) x$pval)))
P <- cbind(P, cor=res$cor$p.value[gs,1], rankcor=res$rankcor$p.value[gs,1] )
P <- cbind(P, camera=res$camera[gs,"PValue"], camera.rnk=res$camera[gs,"PValue"] ) 
##pairs(-log10(P), cex=0.6, pch=20)
pairs(-log10(P), cex=3, pch='.')

par(mfrow=c(3,3))
plot(res$fgsea$NES, -log10(res$fgsea$pval),
  pch=20, cex=0.8, main="fgsea")
ii <- grep("ultragsea",names(res))
for(i in ii) {
  plot(res[[i]]$score, -log10(res[[i]]$pval), 
    pch=20, cex=0.8, main=names(res)[i])
}
plot(res$cor$rho[gs,1], -log10(res$cor$p.value[gs,1]),
  pch=20, cex=0.8, main="gset.cor")
plot(res$rankcor$rho[gs,1], -log10(res$rankcor$p.value[gs,1]),
  pch=20, cex=0.8, main="gset.rankcor")

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

##---------------------------------------------------------------
##----------------------- zGSEA ---------------------------------
##---------------------------------------------------------------
library(limma)

mod <- model.matrix(~y)
m1 <- roast(X, gmt, design=mod, contrast=2)
m2 <- mroast(X, gmt, design=mod, contrast=2)
m3 <- fry(X, gmt, design=mod, contrast=2)

## Default S3 method:
m4 <- cameraPR(fc, gmt, use.ranks = FALSE)

m1 <- m1[gs,]
m2 <- m2[gs,]
m3 <- m3[gs,]
m4 <- m4[gs,]

P <- cbind( fgsea=fres$pval, ultragsea=zres$pval,
  roast=m1$PValue, mroast=m2$PValue, fry=m3$PValue,
  camera=m4$PValue)
pairs(P, pch='.', cex=3)
pairs(-log10(P), pch='.', cex=3)

head(sort(m1$PValue))
head(sort(m2$PValue))
head(sort(m3$PValue))

head(rownames(m1)[order(m1$PValue)])
head(rownames(m2)[order(m2$PValue)])
head(rownames(m3)[order(m3$PValue)])


head(fres[order(fres$pval),])
head(zres[order(zres$pval),])
