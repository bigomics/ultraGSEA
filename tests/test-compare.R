##--------------------------------------------------------
## Comparison of different methods
##--------------------------------------------------------

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

P <- cbind(
  fgsea = f0[ss,]$pval,
  ultragsea = z1[ss,]$pval,
  goat = g1$pval,
  goat.bs = g2$pval,
  goat.gsea = g3$pval
)
pairs(-log10(P))


##--------------------------------------------------------
##--------------------------------------------------------
##--------------------------------------------------------

gs=names(gmt)

f0 <- fgsea::fgsea(gmt, fc, eps=0, gseaParam=0)
f1 <- fgsea::fgsea(gmt, fc, eps=0, gseaParam=1)
f0 <- f0[match(gs,f0$pathway),]
f1 <- f1[match(gs,f1$pathway),]

z0 <- ultragsea(fc, G, method="ztest", format='as.gsea')
z0 <- z0[match(gs,z0$pathway),]
z1 <- ultragsea(fc, G, method="cor", format='as.gsea')
z1 <- z1[match(gs,z1$pathway),]

g1 <- goat(gmt, fc)
g2 <- goat(gmt, fc, method="goat_bootstrap")
g1 <- g1[match(gs,g1$pathway),]
g2 <- g2[match(gs,g2$pathway),]

c1 <- limma::cameraPR(fc, gmt, use.ranks=FALSE)
c2 <- limma::cameraPR(fc, gmt, use.ranks=TRUE)
c1$score <- -log10(c1$PValue) * (-1 + 2*(c1$Direction=="Up"))
c2$score <- -log10(c2$PValue) * (-1 + 2*(c2$Direction=="Up"))
c1 <- c1[gs,]
c2 <- c2[gs,]

ungap <- function(x) ifelse(x<0, x+1, x-1)
S <- cbind(
  ES0=f0$ES, NES0=f0$NES, 
  ES1=f1$ES, NES1=f1$NES,
  ultragsea.ztest = z0$NES,
  ultragsea.cor = z1$NES, 
  goat = g1$score,
  goat.bs = g2$score,
  cameraPR = c1$score,
  cameraPR.rnk = c2$score
)
gsize <- Matrix::colSums(G!=0)[names(gmt)]
pairs(S, pch=20, cex = 0.2*log10(1+gsize)**2 )

cor.test( ungap(f0$NES), g1$score )
cor.test( ungap(f1$NES), z0$NES )
cor.test( g1$score, z0$NES )

P <- cbind(
  fgsea0 = f0$pval,
  fgsea1 = f1$pval,
  ultragsea.ztest = z0$pval,
  ultragsea.cor = z1$pval, 
  goat = g1$pval,
  goat.bs = g2$pval,
  cameraPR = c1$PValue,
  cameraPR.rnk = c2$PValue  
)
pairs(-log(P), pch=20, cex=1)


##----------------------------------------------
##----------------------------------------------
##----------------------------------------------
