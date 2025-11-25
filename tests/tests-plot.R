##install.packages("tidydr")

library(playbase)
source("~/Playground/playbase/dev/include.R",chdir=TRUE)

source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")

pgx <- playdata::GEIGER_PGX
counts  <- pgx$counts
samples <- pgx$samples
head(samples)
y <- samples$activated
fc <- pgx.getMetaMatrix(pgx)$fc[,1]
head(fc)

G <- Matrix::t(playdata::GSETxGENE)
G <- G[,grep("GO",colnames(G))]
gg <- intersect(names(fc),rownames(G))
G <- G[gg,]
fc <- fc[gg]
gmt <- mat2gmt(G)
head(names(gmt))
length(gmt)
system.time(G <- gmt2mat(gmt))
object.size(gmt)/1e6
object.size(G)/1e6

res1 <- fgsea::fgsea(gmt, fc, eps=0)
res2 <- ultragsea(fc, G, method='ztest', format='as.gsea')
res3 <- ultragsea(fc, G, method='cor', format='as.gsea')

gs <- names(gmt)
res1 <- res1[match(gs,res1$pathway),]
res2 <- res2[match(gs,res2$pathway),,]
res3 <- res3[match(gs,res3$pathway),,]
head(res1)
head(res2)
head(res3)

Z <- cbind(res1$NES, res2$NES, res3$NES)
pairs(Z, pch='.', cex=3)
P <- cbind(res1$pval, res2$pval, res3$pval)
pairs(P, pch='.', cex=3)
pairs(-log10(P), pch='.', cex=3)

##BiocManager::install(c("clusterProfiler","enrichplot"))
library(fgsea)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

contrast=1
pgx.makeGseaResult <- function(pgx, contrast, sig=0.05, filter=NULL) {
  gx.meta <- pgx.getMetaFoldChangeMatrix(pgx)
  fc <- gx.meta$fc[,contrast]
  pv <- gx.meta$pv[,contrast]
  qv <- gx.meta$qv[,contrast]  
  names(fc) <- names(pv) <- rownames(gx.meta$fc)
  sig.genes <- names(which(qv <= sig))
  sig.genes <- intersect(sig.genes, rownames(pgx$GMT))
  gs.meta <- pgx$gset.meta$meta[[contrast]]
  if(!is.null(filter)) {
    gs.meta <- gs.meta[grep(filter,rownames(gs.meta)),]
    if(nrow(gs.meta)==0) {
      warning("filtering caused empty geneset table")
      return(NULL)
    }
  }
  gs <- intersect(rownames(gs.meta), colnames(pgx$GMT))
  G <- pgx$GMT[,gs]
  gs.size <- Matrix::colSums(G!=0)
  le.genes <- mat2gmt(G[sig.genes,gs])
  gmt <- mat2gmt(G[,gs])  
  res <- data.table::data.table(
    pathway = rownames(gs.meta),
    pval = gs.meta$meta.p,
    padj = gs.meta$meta.q,
    log2err = NA,
    ES = gs.meta$meta.fx,
    NES = gs.meta$meta.fx,
    size = gs.size,
    leadingEdge = le.genes
  )
  res <- makeGseaResult(res, geneList=fc, geneSets=gmt)
  res
}

res=res3
res <- res[order(-res$NES),]
res <- res[order(res$pval),] 
head(res)

gs <- res$pathway[1]
playbase::gsea.enplot(fc, gmt[[gs]])

detach("package:clusterProfiler", unload=TRUE)
detach("package:enrichplot", unload=TRUE)
library(enrichplot)

## Here we create a gseaResults object from fc+gmt
obj <- makeGseaResult(res, geneList=fc, geneSets=gmt)
obj <- pgx.makeGseaResult(pgx, contrast=1, filter="GO_")
head(obj@result)
pgx.getContrasts(pgx)

obj@result <- obj@result[order(-abs(obj@result$NES)),]
head(obj@result)

cex=1.4
#cex=1
#barplot(obj, showCategory=20) + ggtitle("barplot for GSEA")
dotplot(obj, showCategory=40, x="NES", label_format=100, font.size=12*cex) +
  ggtitle("dotplot for GSEA")
dotplot(obj, showCategory=40, x="GeneRatio", label_format=100, font.size=12*cex) +
  ggtitle("dotplot for GSEA")

## enrichment running score plot
cex=1.2
#cex=0.9
gseaplot2(obj, geneSetID = 1, base_size=20*cex)
gseaplot2(obj, geneSetID = 1:5, base_size=20*cex)

## enrichmap
bp <- pairwise_termsim(obj, showCategory=200)
gsname <- gsub(".*:|\\(.*|GO","",bp@result$Description)
gsname <- gsub("[-_]"," ",gsname)
gsname <- gsub(".*:|\\(.*","",gsname)
names(gsname) <- bp@result$Description

rownames(bp@termsim) <- colnames(bp@termsim) <- gsname[colnames(bp@termsim)]
bp@result$Description <- gsname[bp@result$Description]
names(bp@geneSets) <- gsname[names(bp@geneSets)]
emapplot(bp, showCategory=40, color="NES", nWords=4,
  group=TRUE, group_style="ggforce", label_group_style="ggforce",
  label_format = 10)


devtools::load_all("~/src/enrichplot")

## Bipartite graph
par(mar=c(8,8,8,8))
cnetplot(obj, foldChange=fc, fc_threshold=0, size_item=1.5) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
x11()
cnetplot(obj, foldChange=fc, fc_threshold=1.3, size_item=1.5) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


## Upset plot (hate this...)
upsetplot(obj, n=20) + theme_dose(font.size=22)

## Ridge plot
ridgeplot(obj, showCategory=40, label_format = 120) + theme_minimal(base_size=16)

## Heatplot
p1 <- heatplot(obj, foldChange=fc, showCategory=20, showTop=NULL, label_format=120) + tt
p2 <- heatplot(obj, foldChange=fc, showCategory=20, showTop=80, label_format=120) + tt
plot_list(p1, p2, ncol=1, tag_levels='A', tag_size = 24*cex)


## Cluster tree
obj2 <- pairwise_termsim(obj)
gsname <- gsub(".*:|\\(.*","",colnames(obj2@termsim))
gsname <- gsub(".*:|\\(.*","",gsname)
rownames(obj2@termsim) <- colnames(obj2@termsim) <- toupper(gsname)

devtools::load_all("~/src/enrichplot")
gg <- treeplot(obj2,  
  showCategory = 40, fontsize_tiplab = 5*cex,
  label_format = 20, fontsize_cladelab=8*cex,
  cladelab_offset = 20, tiplab_offset=1,
  color = "NES", size_var="setSize"
)
gg




