##install.packages("goat")
library(goat)

#source("~/Playground/playbase/dev/include.R",chdir=TRUE)
source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")
source("../R/goat.R")

# TODO: change the output directory to some folder on your computer
# e.g. "~/data/goat" on unix systems, "C:/data/goat" on Windows,
# or set to getwd() to write output to the current working directory
output_dir = getwd()

# download an example gene list
datasets = download_goat_manuscript_data(output_dir)
names(datasets)

# download GO gene sets
genesets_asis = download_genesets_goatrepo(output_dir)
genesets_filtered = filter_genesets(genesets_asis, genelist)

## genesets: tibble with genesets, must contain columns 'source',
##           'source_version', 'id', 'name', 'genes', 'ngenes',
##           'ngenes_signif'
head(genesets_filtered)
genesets_filtered$genes[[1]]

## genelist: tibble with genes, must contain column 'gene' and 'test'.
##           gene = character column, which are matched against list
##           column 'genes' in genesets tibble. test = boolean column (you
##           can set all to FALSE if not performing Fisher-exact or
##           hypergeometric test downstream)

head(genelist$gene)
head(genelist$test)
head(genelist$effectsize)

# apply GOAT; score_type indicates the input gene list should be ranked by respective
# "effectsize". This is the recommended setting. Alternatively, this can be set to "pvalue"
# as well but will yield fewer genesets (c.f. GOAT manuscript, real-world benchmarks)
result = test_genesets(genesets_filtered, genelist, method = "goat",
  score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05)

# print first 10 rows of the result table
print(result |> select(source, name, ngenes, pvalue_adjust) |> utils::head(n=10))

str(genesets_filtered)

gmt <- genesets_filtered$genes
names(gmt) <- genesets_filtered$name
fc <- genelist$log2fc
names(fc) <- genelist$gene

source("../R/goat.R")
res1 <- goat(gmt, fc, filter=FALSE, method="goat")
res2 <- goat(gmt, fc, filter=FALSE, method="goat_fitfunction")
res3 <- goat(gmt, fc, filter=FALSE, method="goat_bootstrap")
res2 <- res2[match(res2$pathway, res1$pathway),]
res3 <- res3[match(res3$pathway, res1$pathway),]

S <- cbind( res1$score, res2$score, res3$score)
pairs(S)

P <- cbind( res1$pval, res2$pval, res3$pval)
pairs(-log10(P), pch='.', cex=3)

##----------------------------------------------------
##----------------------------------------------------
##----------------------------------------------------

goat::download_genesets_goatrepo(".")
genesets_asis = download_genesets_goatrepo()
genesets_filtered = filter_genesets(genesets_asis, genelist)
gmt <- genesets_filtered$genes
names(gmt) <- genesets_filtered$name

load("goat_manuscript_datasets.rda",verbose=1)
str(goat_manuscript_datasets)
names(goat_manuscript_datasets)
dim(goat_manuscript_datasets[[1]])

fc.list <- list()
for(d in names(goat_manuscript_datasets)) {
  genelist <- goat_manuscript_datasets[[d]]
  fc <- genelist$log2fc
  names(fc) <- genelist$gene
  fc.list[[d]] <- fc
}

lapply(fc.list, length)
G <- gmt2mat(gmt)
res.goat <- lapply(fc.list, function(f) goat(gmt,f,filter=TRUE))
res.gsea <- lapply(fc.list, function(f) fgsea::fgsea(gmt,f,eps=0))

res.goat <- lapply(res.goat, function(r) r[match(names(gmt),r$pathway),])
res.gsea <- lapply(res.gsea, function(r) r[match(names(gmt),r$pathway),])
res.goat <- lapply(res.goat, function(d) as.data.frame(d, row.names=names(gmt)))
res.gsea <- lapply(res.gsea, function(d) as.data.frame(d, row.names=names(gmt)))

res.ultra <- lapply(fc.list, function(f) ultragsea(f,G,method="cor"))
res.ultraZ <- lapply(fc.list, function(f) ultragsea(f,G,method="ztest"))
res.ultra <- lapply(res.ultra, function(r) r[match(names(gmt),r$pathway),])
res.ultraZ <- lapply(res.ultraZ, function(r) r[match(names(gmt),r$pathway),])

i=1
for(i in 1:length(res.gsea)) {
  ii <- which(is.na(res.gsea[[i]]$NES))
  if(length(ii)) {
    res.gsea[[i]]$NES[ii] <- 0
    res.gsea[[i]]$pval[ii] <- 1
  }
}


cexx <- lapply(fc.list, function(f) {
  gsize <- sapply(gmt, function(s) sum(s %in% names(f)))
  0.2 + log2(1+gsize / mean(gsize))
})

combn(1:4, 2)

pdf("datatest-ultra-vs-fgsea.pdf",7,7)
par(mfrow=c(3,3), mar=c(4,4,2,2))
i=1
for(i in 1:length(res.gsea)) {
  x <- res.gsea[[i]]$NES
  #x <- res.goat[[i]]$score
  y <- res.ultra[[i]]$score

  d <- substring(names(res.gsea)[i],1,20)
  plot(x, y, pch=".",cex=2*cexx[[i]], main=d,
    xlab="fgsea", ylab="ultra (cor test)")
}
par(mfrow=c(3,3), mar=c(4,4,2,2))
for(i in 1:length(res.gsea)) {
  x <- res.gsea[[i]]$NES
  #x <- res.goat[[i]]$score
  y <- res.ultraZ[[i]]$score
  d <- substring(names(res.gsea)[i],1,20)
  plot(x, y, pch=".",cex=2*cexx[[i]], main=d,
    xlab="fgsea", ylab="ultraZ (z-test)")
}
par(mfrow=c(3,3), mar=c(4,4,2,1))
for(i in 1:length(res.gsea)) {
  x <- res.gsea[[i]]$pval
  #x <- res.goat[[i]]$pval
  y <- res.ultra[[i]]$pval
  d <- substring(names(res.gsea)[i],1,20)
  plot( -log10(x), -log10(y),
    pch=".",cex=2*cexx[[i]], main=d,
    xlab="fgsea", ylab="ultra (cor test)")
}
par(mfrow=c(3,3), mar=c(4,4,2,1))
for(i in 1:length(res.gsea)) {
  x <- res.gsea[[i]]$pval
  #x <- res.goat[[i]]$pval
  y <- res.ultraZ[[i]]$pval
  d <- substring(names(res.gsea)[i],1,20)
  plot( -log10(x), -log10(y),
    pch=".",cex=2*cexx[[i]], main=d,
    xlab="fgsea", ylab="ultraZ (z-test)")
}
dev.off()


pdf("goat-datatest-scores.pdf",7,7)
i=1
for(i in 1:length(res.gsea)) {
  x <- res.gsea[[i]]$NES
  y <- res.goat[[i]]$score
  z <- res.ultra[[i]]$score
  d <- substring(names(res.gsea)[i],1,99)
  pairs( cbind(x,y,z),
    labels = c("fgsea","goat","ultra"),
    pch=".",cex=2*cexx[[i]], main=d)
}
dev.off()

pdf("goat-datatest-pvalues.pdf",7,7)
i=1
for(i in 1:length(res.gsea)) {
  x <- res.gsea[[i]]$pval
  y <- res.goat[[i]]$pval
  z <- res.ultra[[i]]$pval
  d <- substring(names(res.gsea)[i],1,99)
  P <- cbind(x,y,z)
  pairs( -log10(P),
    labels = c("fgsea","goat","ultra"),
    pch=".",cex=2*cexx[[i]], main=d)
}
dev.off()

R = res.goat[[1]]
R = R[order(-R$score),]
head(R)

k = R$pathway[2]
res.gsea[[1]][k,]
res.goat[[1]][k,]
res.ultra[[1]][k,]


fgsea::plotEnrichment(gmt[[k]], fc.list[[1]])
source("~/Playground/playbase/dev/include.R",chdir=TRUE)
gsea.enplot(fc.list[[1]], gmt[[k]])

fc <- fc.list[[1]]
gg <- intersect(names(fc),rownames(G))
cor( G[gg,k], fc[gg])
cor( G[gg,k], rank(fc[gg]))
cor( G[gg,k], signedRank(fc[gg]))

R[k,]
head(R)
