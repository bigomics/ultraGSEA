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
X <- playbase::logCPM(counts)
y <- samples$activated

## compute logFC
res <- playbase::gx.limma(X, pheno=y, fdr=1, lfc=0)
fc <- res$logFC
names(fc) <- rownames(res)
head(res)

G <- Matrix::t(playdata::GSETxGENE)
#G <- G[,grep("HALLMARK",colnames(G))]
G <- G[,grep("GO",colnames(G))]
gg <- intersect(rownames(X),rownames(G))
G <- G[gg,]
X <- X[gg,]
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

res1 <- res1[order(res1$pval),]
res2 <- res2[order(res2$pval),]
res3 <- res3[order(res3$pval),]
head(res1)

library(fgsea)
library(dplyr)
library(ggplot2)
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)

contrast=1
pgx.createGseaResult <- function(pgx, contrast, sig=0.05, filter=NULL) {
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
  res <- new.gseaResult(res, geneList=fc, geneSets=gmt)
  res
}

res=res3
res <- res[order(-res$NES),]
res <- res[order(res$pval),] 
head(res)

gs <- res$pathway[1]
playbase::gsea.enplot( fc, gmt[[gs]])

## Here we create a gseaResults object from fc+gmt
obj <- new.gseaResult(res, geneList=fc, geneSets=gmt)
obj <- pgx.createGseaResult(pgx, contrast=1, filter="GO_")
head(obj@result)
pgx.getContrasts(pgx)

obj@result <- obj@result[order(obj@result$pvalue),]
obj@result <- obj@result[order(-abs(obj@result$NES)),]
head(obj@result)

#barplot(obj, showCategory=20) + ggtitle("barplot for GSEA")
dotplot(obj, showCategory=40, x="NES", label_format=100, font.size=12) +
  ggtitle("dotplot for GSEA")
dotplot(obj, showCategory=40, x="GeneRatio", label_format=100, font.size=12) +
  ggtitle("dotplot for GSEA")

## enrichment running score plot
gseaplot2(obj, geneSetID = 1, base_size=20)
gseaplot2(obj, geneSetID = 1:5, base_size=20)

## enrichmap 
bp <- pairwise_termsim(obj, showCategory=200)
emapplot(bp)
emapplot(bp,
  showCategory=40,
  color = "NES",
  group_category = TRUE,
  group_legend = FALSE,
  cex_label_group = 1,
  nWords = 4,
  label_format = 20
)

## Bipartite graph
cnetplot(obj, foldChange=fc)

## Upset plot (hate this...)
upsetplot(obj, n=20) + theme_dose(font.size=22)

## Ridge plot
ridgeplot(obj, showCategory=40, label_format = 120) + theme_minimal(base_size=16)

## Heatplot
p1 <- heatplot(obj, foldChange=fc, showCategory=20, label_format=120) +
  theme_minimal(base_size=18)
p2 <- heatplot(obj, showCategory=20, label_format=120) +
  theme_minimal(base_size=18)
plot_list(p1, p2, ncol=1, tag_levels = 'A')

## Cluster tree
obj2 <- pairwise_termsim(obj)
gsname <- gsub(".*:|\\(.*","",colnames(obj2@termsim))
gsname <- gsub(".*:|\\(.*","",gsname)
rownames(obj2@termsim) <- colnames(obj2@termsim) <- toupper(gsname)
gg <- treeplot(obj2,
  cluster.params = list(method = "ward.D", n=3, color=NULL, label_words_n=4,
    label_format=100),
  showCategory = 40, fontsize = 4,
  label_format=20, size=6, color="NES")
gg
##gg + ggtree::geom_tiplab(size=7)


