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

u1 <- ultragsea(fc, G, format='simple')
u2 <- ultragsea(fc, G, format='as.gsea')
u3 <- ultragsea(fc, G, format='as.gsea2')

z1 <- fc_ztest(fc, G, zmat=TRUE)
z2 <- fc_ztest( cbind(fc,fc), G, zmat=TRUE)
plot( z1$z_statistic, z2$z_statistic[,1])
plot( z1$p_value, z2$p_value[,1])

F <- cbind(fc,fc,fc,fc,fc)
F <- do.call(cbind, rep(list(fc),1000))
dim(F)
system.time(z1 <- apply(F, 2, function(f) fc_ztest(f, G)))
system.time(z2 <- fc_ztest(F, G))



source("../R/goat.R")
g1 <- goat(gmt, fc, method="goat", filter=FALSE)
head(g1)


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
