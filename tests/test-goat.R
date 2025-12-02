##install.packages("goat")
library(goat)

#source("~/Playground/playbase/dev/include.R",chdir=TRUE)
source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")
source("../R/goat.R")

# download an example gene list
output_dir = "/tmp"
datasets = download_goat_manuscript_data(output_dir)
names(datasets)
sapply(datasets, nrow)

# download GO gene sets
genesets_asis = download_genesets_goatrepo(output_dir)
genelist <- datasets[[3]]
genesets_filtered = filter_genesets(genesets_asis, genelist,
  min_overlap=10, max_overlap=1500)

gmt <- genesets_filtered$genes
names(gmt) <- genesets_filtered$name
range(sapply(gmt,length))
fc <- genelist$log2fc
names(fc) <- genelist$gene

source("../R/goat.R")
res <- goat(gmt, fc, filter=FALSE, method="goat")
res <- goat(gmt, fc, filter=FALSE, method="goat_fitfunction")
res <- goat(gmt, fc, filter=FALSE, method="goat_bootstrap")
head(res)

names(datasets)
dim(datasets[[1]])

fc.list <- list()
for(d in names(datasets)) {
  genelist <- datasets[[d]]
  fc <- genelist$log2fc
  names(fc) <- genelist$gene
  fc.list[[d]] <- fc
}

lapply(fc.list, length)
range(sapply(gmt, length))
G <- gmt2mat(gmt)

res.goat <- lapply(fc.list, function(f) goat(gmt,f,filter=TRUE))
res.goat <- lapply(res.goat, function(r) r[match(names(gmt),r$pathway),])
res.goat <- lapply(res.goat, function(d) as.data.frame(d, row.names=names(gmt)))

res.gsea <- lapply(fc.list, function(f) fgsea::fgsea(gmt,f,eps=0))
#res.gsea <- lapply(fc.list, function(f) fgsea::fgsea(gmt,f,eps=0,nPermSimple=10000))
res.gsea <- lapply(res.gsea, function(r) r[match(names(gmt),r$pathway),])
res.gsea <- lapply(res.gsea, function(d) as.data.frame(d, row.names=names(gmt)))

res.ultra <- lapply(fc.list, function(f) ultragsea(f,G,method="cor"))
res.ultraB <- lapply(fc.list, function(f) ultragsea(f,G,method="cor",cor0=3))
res.ultraZ <- lapply(fc.list, function(f) ultragsea(f,G,method="ztest"))

res.ultra <- lapply(res.ultra, function(r) r[match(names(gmt),r$pathway),])
res.ultraB <- lapply(res.ultraB, function(r) r[match(names(gmt),r$pathway),])
res.ultraZ <- lapply(res.ultraZ, function(r) r[match(names(gmt),r$pathway),])

cexx <- lapply(fc.list, function(f) {
  gsize <- sapply(gmt, function(s) sum(s %in% names(f)))
  0.2 + log2(1+gsize / mean(gsize))
})

score.list <- list()
pv.list <- list()
for(i in 1:length(res.gsea)) {
  score.list[[i]] <- cbind( 
    gsea =  res.gsea[[i]]$NES,
    goat =  res.goat[[i]]$score,
    ultra = res.ultra[[i]]$score,
    ultraB = res.ultraB[[i]]$score,
    ultraZ = res.ultraZ[[i]]$score
  )
  pv.list[[i]] <- cbind( 
    gsea =  res.gsea[[i]]$pval,
    goat =  res.goat[[i]]$pval,
    ultra = res.ultra[[i]]$pval,
    ultraB = res.ultraB[[i]]$pval,
    ultraZ = res.ultraZ[[i]]$pval
  )
}

do.call(rbind,lapply(score.list, function(S) colSums(is.na(S))))
do.call(rbind,lapply(pv.list, function(P) colSums(is.na(P))))

ii <- which(is.na(score.list[[1]][,3]))
head(res.ultra[[1]][ii,])
head(res.ultraZ[[1]][ii,])

combi <- list(
  c("gsea","goat"),
  c("gsea","ultra"),
  c("gsea","ultraZ"),
  c("goat","ultra"),
  c("goat","ultraZ"),
  c("ultraZ","ultra"),
  c("ultraZ","ultraB"),  
  c("ultra","ultraB")
)

k=6
pdf("datatest-ultra-vs-fgsea-vs-goat-max1500.pdf",7,7)
for(k in 1:length(combi)) {
  a <- combi[[k]][1]
  b <- combi[[k]][2]
  par(mfrow=c(3,3), mar=c(4,4,2,2))
  i=1
  for(i in 1:length(res.gsea)) {
    x <- score.list[[i]][,a]
    y <- score.list[[i]][,b]
    s <- score.list[[i]][,"gsea"]    
    d <- substring(names(res.gsea)[i],1,20)
    col1 <- 'black'
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
    col1 <- c('black','red')[1+1*(x==0|y==0)]
    col1[is.na(s)|s==0] <- 'green2'  ## missing in gsea    
    plot(x, y, pch=".",cex=2*cexx[[i]], main=d,
      xlab = paste(a,"(score)"),
      ylab = paste(b,"(score)"),
      col = col1)
    jj <- which(col1=="red")
    points(x[jj], y[jj], pch=".",
      cex=2*cexx[[i]][jj], col = 'red')
  }
  par(mfrow=c(3,3), mar=c(4,4,2,1))
  for(i in 1:length(res.gsea)) {
    s <- score.list[[i]][,"gsea"]
    x <- pv.list[[i]][,a]
    y <- pv.list[[i]][,b]
    d <- substring(names(res.gsea)[i],1,20)
    x[is.na(x)] <- 1
    y[is.na(y)] <- 1
    col1 <- c('black','red')[1+1*(x==1|y==1)]
    col1[is.na(s)|s==0] <- 'green2'  ## missing in gsea        
    x <- -log10(x)
    y <- -log10(y)
    plot(x, y,
      pch=".",cex=2*cexx[[i]], main=d,
      xlab = paste(a,"(-log10p)"),
      ylab = paste(b,"(-log10p)"),
      col=col1)
    jj <- which(col1=="red")
    points(x[jj], y[jj], pch=".",
      cex=2*cexx[[i]][jj], col = 'red')    
  }
}
dev.off()

pdf("datatest-scores-max1500.pdf",7,7)
i=1
for(i in 1:length(res.gsea)) {
  Z <- score.list[[i]]
  col1 <- c("black","red")[1 + 1*(rowSums(is.na(Z))>0)]
  col1[is.na(Z[,'gsea'])] <- 'green2'
  Z[is.na(Z)] <- 0
  d <- substring(names(res.gsea)[i],1,99)
  pairs( Z, pch=".",cex=2*cexx[[i]], col=col1, main=d)
}
dev.off()

pdf("datatest-pvalues-max1500.pdf",7,7)
i=1
for(i in 1:length(res.gsea)) {
  Z <- score.list[[i]]  
  P <- pv.list[[i]]
  d <- substring(names(res.gsea)[i],1,99)
  col1 <- c("black","red")[1 + 1*(rowSums(is.na(P))>0)]
  col1[is.na(Z[,'gsea'])] <- 'green2'  
  P[is.na(P)] <- 1
  pairs( -log10(P),
    labels = c("fgsea","goat","ultra","ultraZ"),
    pch=".",cex=2*cexx[[i]], col=col1,
    main=d)
}
dev.off()

##----------------------------------
##----------------------------------
##----------------------------------


ii <- which(is.na(res.gsea[[1]]$NES))
ii <- which(is.na(res.goat[[1]]$score))
R = res.ultra[[1]][ii,]
R = R[order(-R$score),]
head(R)

source("~/Playground/playbase/dev/include.R",chdir=TRUE)

pdf("datatest-gseaNA.pdf",w=12,h=7)
par(mfrow=c(4,5), mar=c(2,3,3,1))
ppw <- head(R$pathway,40)
for(pw in ppw) {
  ##fgsea::plotEnrichment(gmt[[k]], fc.list[[1]])
  gsea.enplot(fc.list[[1]], gmt[[pw]], main=pw)
}
dev.off()

