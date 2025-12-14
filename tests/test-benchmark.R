##BiocManager::install("GSEABenchmarkeR")

library(GSEABenchmarkeR)
geo <- loadEData("geo2kegg")
geo <- maPreproc(geo)
length(geo)

geo <- runDE(geo, de.method="limma", padj.method="BH")
head(rowData(geo[[1]], use.names=TRUE))

library(EnrichmentBrowser)
kegg.gmt <- getGenesets(org="hsa", db="kegg")
go.gmt <- getGenesets(org="hsa", db="go")
gmt <- kegg.gmt
gmt <- c(kegg.gmt, go.gmt)
length(gmt)
matG <- ultragsea::gmt2mat(gmt)
dim(matG)

##--------------------------------------------
## Benchmarking
##--------------------------------------------
source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")
source("../R/goat.R")
source("../R/fastFET.R")

length(geo)
se <- geo[[1]]
table(se$GROUP)
rd <- rowData(se, use.names=TRUE)
head(rd)

run.fgsea <- function(se,gs) {
  rd <- rowData(se, use.names=TRUE)
  fc <- rd$FC
  names(fc) <- rownames(rd)
  res <- fgsea::fgsea(gs, fc)
  pv <- res$pval
  names(pv) <- res$pathway
  pv
}

run.goat <- function(se,gs) {
  rd <- rowData(se, use.names=TRUE)
  fc <- rd$FC
  names(fc) <- rownames(rd)
  res <- ultragsea::goat(gs, fc)
  pv <- res$pval
  names(pv) <- res$pathway
  pv
}

run.cameraPR <- function(se,gs) {
  rd <- rowData(se, use.names=TRUE)
  fc <- rd$FC
  names(fc) <- rownames(rd)
  res <- limma::cameraPR(fc,gs)
  pv <- res$PValue
  names(pv) <- rownames(res)
  pv
}

run.ultraC <- function(se,gs) {
  rd <- rowData(se, use.names=TRUE)
  fc <- rd$FC
  names(fc) <- rownames(rd)
  res <- ultragsea::ultragsea(matG, fc, method='cor', alpha=0)
  pv <- res$pval
  names(pv) <- rownames(res)
  pv
}

run.ultraZ <- function(se,gs) {
  rd <- rowData(se, use.names=TRUE)
  fc <- rd$FC
  names(fc) <- rownames(rd)
  res <- ultragsea::ultragsea(matG, fc, method='ztest', alpha=0)
  pv <- res$pval
  names(pv) <- rownames(res)
  pv
}

run.gsetcor <- function(se, gs) {
  rd <- rowData(se, use.names=TRUE)
  fc <- rd$FC
  names(fc) <- rownames(rd)
  res <- ultragsea::gset.cor(matG, fc,
    compute.p=TRUE, corshrink=3)
  pv <- res$p.value[,1]
  pv
}

run.fastFET <- function(se,gs) {
  rd <- rowData(se, use.names=TRUE)
  head(rd)
  sig <- rownames(rd)[ rd$PVAL < 0.05 ]
  if(length(sig)<20) sig <- head(rownames(rd)[order(rd$PVAL)],100)
  #  res <- playbase::gset.fisher(sig, gs, background=rownames(rd),
  #    fdr=1, nmin=1, min.genes=1, max.genes=999999)
  res <- gset.fastFET(sig, matG, bg=rownames(rd), method=2)
  res <- res[match(names(gs),rownames(res)),]
  pv <- res$p.value
  names(pv) <- rownames(res)
  pv
}

##-----------------------------------------------
## run methods
##-----------------------------------------------

par.methods = list(ora="ora", camera="camera")
par.methods = list(
  #ora = "ora",
  fGSEA = run.fgsea,
  #camera = "camera",
  cameraPR = run.cameraPR,
  fastFET = run.fastFET,
  GOAT = run.goat,
  #gsetcor = run.gsetcor,
  ultraZ = run.ultraZ,
  ultraC = run.ultraC  
)
methods <- names(par.methods)
methods

res.dir <- tempdir()
res <- runEA(
  geo,
  methods = par.methods,
  gs = gmt,
  save2file = TRUE,
  out.dir = res.dir
)

##-----------------------------------------------
## gather results
##-----------------------------------------------

png("figures/benchmark-runtime.png",w=650, h=450, pointsize=20)
par(mar=c(5,4.5,2,1))
## plot runtime
ea.rtimes <- readResults(
  res.dir, names(geo), 
  methods = methods,
  type = "runtime")
ea.rtimes
bpPlot(ea.rtimes, what="runtime")
title("Runtime",cex.main=1.3)
dev.off()

## Number of significant genesets
ea.ranks <- readResults(
  res.dir, names(geo), 
  methods = methods,
  type = "ranking"
)
i=j=1
for(i in 1:length(ea.ranks)) {
  for(j in 1:length(ea.ranks[[i]])) {
    pv <- ea.ranks[[i]][[j]]$PVAL
    ii <- which(is.na(pv))
    if(length(ii)) ea.ranks[[i]][[j]]$PVAL[ii] <- 1
  }
}
lengths(ea.ranks)
sig.sets <- evalNrSigSets(ea.ranks, alpha=0.05, padj="fdr")
sig.sets <- evalNrSigSets(ea.ranks, alpha=0.05, padj="none")
bpPlot(sig.sets, what="sig.sets")

## MALA relevance ranking
data.dir <- system.file("extdata", package="GSEABenchmarkeR")
mala.kegg.file <- file.path(data.dir, "malacards", "KEGG.rds")
mala.kegg <- readRDS(mala.kegg.file)
sapply(mala.kegg, nrow)
d2d.file <- file.path(data.dir, "malacards", "GseId2Disease.txt")
d2d.map <- readDataId2diseaseCodeMap(d2d.file)
all.kegg.res <- evalRelevance(ea.ranks, mala.kegg, d2d.map[names(geo)])
bpPlot(all.kegg.res, what="rel.sets")

## relevance vs runtime
rt <- log10(do.call(cbind, ea.rtimes))
rt <- rt - rowMeans(rt)
avg.rt <- colMeans(rt)
avg.rel <- colMeans(all.kegg.res)
avg.rt <- 1/(10**avg.rt)

png("figures/benchmark-relevance-vs-runtime.png",w=700, h=500, pointsize=20)
par(mar=c(5,4.5,2,1))
plot( avg.rt, avg.rel, pch=19,
  xlab="relative speedup (larger is faster)",
  ylab="avg. relevance (larger is better)",
  las=1, log='x',
  xlim = range(avg.rt) * c(0.8,1.2),
  ylim = range(avg.rel) + c(-1.2,1)*2
)
abline(v=1,lty=3)
text( avg.rt, avg.rel, labels=names(avg.rt), cex=1, pos=2)
title("relevance vs. runtime",cex.main=1.3)
dev.off()

