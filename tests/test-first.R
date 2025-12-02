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
gmt <- gmt[sapply(gmt,length)>0]
G <- gmt2mat(gmt)

gg <- intersect(rownames(G),names(fc))
G <- G[gg,]
fc <- fc[gg]

length(fc)
length(gmt)
range(sapply(gmt,length))

plot(sort(fc))
fc1 = fc
fc1 = fc - mean(fc)
res1 <- fgsea::fgsea(gmt, fc1, eps=0)
res1 <- res1[match(colnames(G),res1$pathway),]
res2 <- ultragsea(fc1, G, method='cor')
res3 <- ultragsea(fc1, G, method='ztest')
res4 <- ultragsea(fc1, G, method='ttest')
res5 <- ultragsea(fc1, G, method='goat')
res5 <- res5[match(colnames(G),res5$pathway),]

S <- cbind( res1$NES, res2$score, res3$score, res4$score, res5$score)
colnames(S) <- c("fgsea","cor","ztest","ttest","goat")
pairs(S)

P <- cbind( res1$pval, res2$pval, res3$pval, res4$pval, res5$pval)
colnames(P) <- c("fgsea","cor","ztest","ttest","goat")
pairs(-log10(P))

r1 <- gset.cor(fc1, G, compute.p=TRUE)
r2 <- fc_ztest(fc1, G)
r3 <- matrix_onesample_ttest(cbind(fc1), G)
r4 <- matrix_onesample_ztest(cbind(fc1), G)

S <- cbind( cor=r1$rho[,1], fcz=r2$z, tt=r3$t[,1], zt=r4$z[,1])
pairs(S)

P <- cbind( cor=r1$p.value[,1], fcz=r2$p, tt=r3$p[,1], zt=r4$p[,1])
pairs(-log10(P))
