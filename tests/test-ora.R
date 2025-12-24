##--------------------------------------------------------
## Comparison of ORA
##--------------------------------------------------------

source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")
source("../R/goat.R")
source("../R/fastFET.R")

library(fgsea)
library(goat)

G <- Matrix::t(playdata::GSETxGENE)
fc <- scale(Matrix::rowMeans(G!=0))[,1]
gmt <- mat2gmt(G)
table(abs(fc)>1)

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

sig.genes <- names(which(abs(fc) > 1))
genes=sig.genes
gs <- gmt[[6]]
bg = names(fc)

bg = names(fc)
system.time(res1 <- playbase::gset.fisher(sig.genes, gmt, fdr=1, background=bg))
system.time(res2 <- playbase::gset.fisher(sig.genes, G, fdr=1, background=bg))
system.time(res3 <- gset.fastFET(sig.genes, G, bg=bg, report.genes=FALSE))
system.time(res4 <- manual.ft(sig.genes, gmt, bg=bg))

res1 <- res1[match(names(gmt),rownames(res1)),]
res2 <- res2[match(names(gmt),rownames(res2)),]
res3 <- res3[match(names(gmt),rownames(res3)),]
res4 <- res4[match(names(gmt),rownames(res4)),]
head(res1$p.value)
head(res2$p.value)

P <- cbind( res1$p.value, res2$p.value, res3$p.value, res4[,"pval"])
pairs(-log10(P), cex=0.4)


##-----------------------------------------------------
##-----------------------------------------------------
##-----------------------------------------------------
