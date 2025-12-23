##===========================================================
##
## Benchmarking using GSEABenchmarkeR
##
##===========================================================

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
range(sapply(gmt,length))
gmt <- gmt[sapply(gmt,length)>=10 & sapply(gmt,length) < 500]
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
source("../R/plaidtest.R")

length(geo)
head(names(geo))
se <- geo[[1]]
table(se$GROUP)
rd <- rowData(se, use.names=TRUE)
head(rd)
gs <- gmt
X <- assay(se)
y <- se$GROUP
ref=0
G=matG

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
  res <- goat(gs, fc, verbose=0)
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
  res <- ultragsea::gset.cortest(matG, fc,
    compute.p=TRUE, corshrink=3)
  pv <- res$p.value[,1]
  pv
}

run.fastFET <- function(se,gs) {
  rd <- rowData(se, use.names=TRUE)
  sig <- rownames(rd)[ rd$PVAL < 0.05 ]
  if(length(sig)<20) sig <- head(rownames(rd)[order(rd$PVAL)],100)
  res <- gset.fastFET(sig, matG, bg=rownames(rd))
  res <- res[match(names(gs),rownames(res)),]
  pv <- res$p.value
  names(pv) <- rownames(res)
  pv
}

run.plaid <- function(se,gs) {
  rd <- rowData(se, use.names=TRUE)
  fc <- rd$FC
  names(fc) <- rownames(rd)
  y <- se$GROUP
  X <- assay(se)
  res <- plaid.ttest(X, matG, y, ref=0)
  res <- res[match(names(gs),rownames(res)),]
  pv <- res$pval
  names(pv) <- rownames(res)
  pv
}

run.dual <- function(se,gs) {
  y <- se$GROUP
  X <- assay(se)
  res <- plaid.dualtest(X, matG, y, ref=0, test1="ttest")
  res <- res[match(names(gs),rownames(res)),]
  pv <- res$pval
  names(pv) <- rownames(res)
  pv
}

run.dual2 <- function(se,gs) {
  y <- se$GROUP
  X <- assay(se)
  res <- plaid.dualtest(X, matG, y, ref=0, test1="cortest")
  res <- res[match(names(gs),rownames(res)),]
  pv <- res$pval
  names(pv) <- rownames(res)
  pv
}

run.dual3 <- function(se,gs) {
  y <- se$GROUP
  X <- assay(se)
  res <- plaid.dualtest(X, matG, y, ref=0, pbalance=FALSE)
  res <- res[match(names(gs),rownames(res)),]
  pv <- res$pval
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
  ultraC = run.ultraC,
  plaid = run.plaid,
  dual = run.dual,
  dual2 = run.dual2,
  dual3 = run.dual3
)
methods <- names(par.methods)
methods

## RunEA 
res.dir <- tempdir()
res <- runEA(
  geo,
  methods = par.methods,
  gs = gmt,
  save2file = TRUE,
  out.dir = res.dir
)

## Run manually (better timing)
length(geo)
head(names(geo))
rtimes <- list()
k=names(geo)[1]
for(k in names(geo)) {
  cat("computing ",k,"...\n")
  se <- geo[[k]]
  tt <- peakRAM::peakRAM(
    run.fgsea(se,gs),
    run.cameraPR(se,gs),
    run.fastFET(se,gs),
    run.goat(se,gs),
    run.ultraZ(se,gs),
    run.ultraC(se,gs),
    run.plaid(se,gs),      
    run.dual(se,gs),
    run.dual2(se,gs),
    run.dual3(se,gs)          
  )
  tt2 <- tt[,2]
  names(tt2) <- methods
  rtimes[[k]] <- tt2
}

RT <- do.call(cbind, rtimes)
ea.rtimes <- apply(RT, 1, function(x)
  array(x,dimnames=list(colnames(RT))), simplify=FALSE)
save(ea.rtimes, file="eartimes.rda")


##-----------------------------------------------
## gather results
##-----------------------------------------------

png("figures/benchmark-runtime.png",w=700, h=550, pointsize=20)
par(mar=c(5.8,4.5,2,1))
## plot runtime
## ea.rtimes <- readResults(
##   res.dir, names(geo), 
##   methods = methods,
##   type = "runtime")
load(file="eartimes.rda",verbose=1)
#bpPlot(ea.rtimes, what="runtime")
sel <- order(sapply(lapply(ea.rtimes,log),mean))
boxplot(ea.rtimes[sel], log='y', col=2:10,
  ylab="runtime (sec)", las=2, srt=45)
title("Runtime",cex.main=1.3)
dev.off()

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

png("figures/benchmark-relevance-vs-runtime.png",w=650, h=600, pointsize=20)
par(mar=c(5,4.5,2,1))
plot( avg.rt, avg.rel, pch=19, cex=1.2,
  xlab="relative speedup (larger is faster)",
  ylab="avg. relevance (larger is better)",
  las=1, log='x',
  xlim = range(avg.rt) * c(0.8,1.2),
  ylim = range(avg.rel) + c(-1.2,1)*2
)
abline(v=1,lty=3)
text( avg.rt, avg.rel, labels=names(avg.rt), cex=1.05, pos=3)
title("relevance vs. runtime",cex.main=1.3)
dev.off()

