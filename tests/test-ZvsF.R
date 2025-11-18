library(playbase)
source("~/Playground/playbase/dev/include.R",chdir=TRUE)
source("../R/zgsea.R")

counts  <- playbase::COUNTS
samples <- playbase::SAMPLES
head(samples)
X <- logCPM(counts)
y <- samples$activated

## compute logFC
res <- gx.limma(X, pheno=y, fdr=1, lfc=0)
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
head(names(gmt))
length(gmt)

f0 <- fgsea::fgsea(gmt, fc)
f1 <- fgsea::fgseaSimple(gmt, fc, nperm=10)
f2 <- fgsea::fgseaMultilevel(gmt, fc, nPermSimple=100, eps=1)

par(mfrow=c(3,3))
f2 <- f2[match(f1$pathway,f2$pathway),]
head(f1)
head(f2)
plot(f1$NES, f2$NES)

z1 <- zgsea(fc, G)
head(z1)
z1 <- z1[match(f1$pathway,z1$pathway),]
plot(f1$NES, z1$NES)

z2 <- zgsea(sign(fc)*abs(fc)**2, G)
head(z2)
z2 <- z2[match(f1$pathway,z2$pathway),]
plot(f1$NES, z2$NES)

##--------------------------------------------------------
##--------------------------------------------------------
##--------------------------------------------------------

library(Rcpp)
sourceCpp('../R/KS.cpp')

ks.stat <- function(x,y) {
  n <- length(x)
  w <- c(x,y)
  z <- order(w)
  w <- rep(-1,length(w))
  w[which(z <= (n-1))] <- 1
  max(abs(cumsum(w)))/n
}

dim(G)
length(fc)

ks <- rep(NA, ncol(G))
kp <- rep(NA, ncol(G))
names(ks) <- colnames(G)
names(kp) <- colnames(G)
for(j in 1:ncol(G)) {
  ii <- which(G[,j]!=0)
  if(length(ii)==0) next
  ksign <- -1 + 2*( mean(fc[ii]) > mean(fc[-ii]) )
  kt <- ks.test(x=fc[ii], y=fc[-ii])
  ks[j] <- ksign * kt$statistic
  kp[j] <- kt$p.value
}

kss <- rep(NA, ncol(G))
names(kss) <- colnames(G)
for(j in 1:ncol(G)) {
  ii <- which(G[,j]!=0)
  if(length(ii) <= 2) next
  ksign <- -1 + 2*( mean(fc[ii]) > mean(fc[-ii]) )
  kss[j] <- ksign * ks.stat(x=fc[ii], y=fc[-ii])
}

f0 <- fgsea::fgsea(gmt, fc, eps=0)
f0 <- f0[match(names(ks),f0$pathway),]

z0 <- zgsea(fc, G)
z0 <- z0[match(names(ks),z0$pathway),]

K <- cbind(ES=f0$ES, NES=f0$NES, zES=z0$NES, ks=ks, kss=kss)
pairs(K, pch='.', cex=3)

P <- cbind(fgsea=f0$pval, KS=kp, zgsea=z0$pval)
pairs(-log(P), pch='.', cex=3)


mean(abs(f0$NES), na.rm=TRUE)

plot(ks, f0$ES)
plot(ks, f0$NES)


##----------------------------------------------
##----------------------------------------------
##----------------------------------------------

library(Rcpp)
sourceCpp('../R/KS.cpp')
ks.stat <- function(x,y) {
  n <- length(x)
  w <- c(x,y)
  z <- order(w)
  w <- rep(-1,length(w))
  w[which(z <= (n-1))] <- 1
  max(abs(cumsum(w)))/n
}

x = seq(0,1,0.05)
y = seq(0.5,1,0.05)
ks.test(x,y)
ks.test(y,x)
KSxy(x,y)
KSxy(y,x)
ks.stat(x,y)
ks.stat(y,x)


plot(cumsum(w))

##----------------------------------------------
##----------------------------------------------
##----------------------------------------------

library(Rcpp)
sourceCpp('../R/KS.cpp')

set.seed(1942)
mt <- matrix(rnorm(400*30000), nrow=30000)
dim(mt)

results <- rep(0, nrow(mt))
for (i in 2:nrow(mt)) {
  results[i] <- ks.test(x = mt[i - 1, ], y = mt[i, ])$statistic
}

result <- runKS(t(mt))
all.equal(result, results)

result
results
