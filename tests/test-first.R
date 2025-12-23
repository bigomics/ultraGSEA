library(playbase)
library(limma)

source("~/Playground/playbase/dev/include.R",chdir=TRUE)
source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/goat.R")
source("../R/utils.R")
source("../R/plaidtest.R")

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
F <- cbind(fc)

length(fc)
length(gmt)
range(fc)
range(sapply(gmt,length))

plot(sort(fc))
fc1 <- scale(fc)[,1]
#fc1 <- fc1 + 1

plot(sort(fc1));abline(h=0)
summary(fc1)
table(sign(fc1))

res1 <- fgsea::fgsea(gmt, fc1, eps=0)
res1 <- res1[match(colnames(G),res1$pathway),]

res2 <- ultragsea(G, fc1, method='cor')
res3a <- ultragsea(G, fc1, method='ztest', alpha=0)
res3b <- ultragsea(G, fc1, method='ztest', alpha=1)
res3c <- ultragsea(G, fc1, method='ztest', alpha=0, center=FALSE)
res3d <- ultragsea(G, fc1, method='ztest', alpha=1, center=FALSE)
res4a <- ultragsea(G, fc1, method='ttest', alpha=0)
res4b <- ultragsea(G, fc1, method='ttest', alpha=1)
res5 <- ultragsea(G, fc1, method='goat')
res6 <- ultragsea(G, fc1, method='camera')


table(sign(res1$NES))
table(sign(res2$score))
table(sign(res3a$score))
table(sign(res3b$score))
table(sign(res3c$score))
table(sign(res3d$score))
table(sign(res4a$score))
table(sign(res4b$score))
table(sign(res5$score))

S <- cbind( res1$NES, res2$score,
  res3a$score, res3b$score, res3c$score, res3d$score,
  res4a$score, res4b$score, res5$score,
  res6$score)
colnames(S) <- c("fgsea","cor",
  "ztest.a0","ztest.a1", "ztestZ.a0","ztestZ.a1",
  "ttest.a0","ttest.a1","goat","cameraPR")
pairs(S, panel=function(x,y) {
  abline(h=0,v=0,lty=2,lwd=0.7)
  points(x,y,cex=0.4,pch=19)
})

P <- cbind( res1$pval, res2$pval,
  res3a$pval, res3b$pval,res3c$pval, res3d$pval,
  res4a$pval, res4b$pval, res5$pval,res6$pval)
colnames(P) <- c("fgsea","cor",
  "ztest.a0","ztest.a1", "ztestZ.a0","ztestZ.a1",
  "ttest.a0","ttest.a1","goat","cameraPR")
pairs(-log10(P), cex=0.6)
cor(-log10(P))

res <- ultragsea(G, fc1, method="cor", alpha=0)

source("../R/ultragsea.R")
res1 <- fgsea::fgsea(gmt, fc1, eps=0)
colnames(res1)

res2 <- fgsea(gmt, fc1, nperm=10)
res3 <- fgsea(gmt, fc1, nperm=20)
res4 <- fgsea(gmt, fc1, nperm=1000)
N <- cbind(res1$NES, res2$NES, res3$NES, res4$NES)
pairs(N)

P <- cbind(res1$pval, res2$pval, res3$pval, res4$pval)
pairs(-log(P))

head(res1)
head(res3)
colnames(res1)
colnames(res3)


G1 <- Matrix::Matrix(G, sparse=TRUE)
