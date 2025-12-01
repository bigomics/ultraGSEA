# Comparing ultragsea with fGSEA

## Introduction

ultragsea is a novel, ultrafast and memory optimized gene set enrichment
scoring algorithm much like fGSEA. ultragsea typically is 10-100x faster
than fGSEA but demonstrates highly correlated enrichment scores and
highly similar p-values.

## Preparing the sparse gene set matrix

``` r
library(fgsea)
library(goat)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ultragsea)
data(examplePathways)
data(exampleRanks)
fc = exampleRanks
length(fc)
#> [1] 12000
gmt = examplePathways
range(sapply(gmt,length))
#> [1]    1 2366
gmt <- lapply(gmt, function(s) intersect(s,names(fc)))
gmt <- gmt[sapply(gmt,length)>10]
G <- gmt2mat(gmt)
```

## Running the methods

For benchmarking we run some of the currently fastest gene set
enrichment algorithms: fGSEA (fast GSEA), cameraPR (pre-ranked Camera
from the limma package), goat (Gene set Ordinal Association Test). In
particular we are interested as how these compare to fGSEA in runtime,
score and p-value, as fGSEA is perhaps currently the most widely used
algorithm for pre-ranked enrichment testing.

``` r
if(0) {
  ## run this if you want to augment the number of gene sets.
  gmt <- rep(gmt,150)
  length(gmt)
  names(gmt) <- make.unique(names(gmt))
  G <- do.call(cbind, rep(list(G),150))
  colnames(G) <- make.unique(colnames(G))
  dim(G)
}
```

``` r
library(peakRAM)
gs <- names(gmt)
tt <- peakRAM::peakRAM(
  res.fgsea <- fgsea::fgsea(gmt, fc, eps=0),
  res.camera <- limma::cameraPR(fc, gmt, use.ranks=FALSE)[gs,],
  res.ultragsea.z <- ultragsea(fc, G, method='ztest')[gs,],
  res.ultragsea.c <- ultragsea(fc, G, method='cor')[gs,],
  res.cor <- gset.cor(fc, G, compute.p=TRUE, use.rank=FALSE),
  res.ztest <- fc_ztest(fc, G),
  res.goat <- goat(gmt, fc)
)
#> filtering genesets...
res.fgsea <- res.fgsea[match(gs,res.fgsea$pathway),]
res.goat <- res.goat[match(gs,res.goat$pathway),]
tt[,1] <- gsub("res[.]|<-.*","",tt[,1])
kableExtra::kable(tt)
```

| Function_Call | Elapsed_Time_sec | Total_RAM_Used_MiB | Peak_RAM_Used_MiB |
|:--------------|-----------------:|-------------------:|------------------:|
| fgsea         |            4.320 |                4.0 |              43.9 |
| camera        |            0.320 |                0.8 |              42.6 |
| ultragsea.z   |            0.111 |                1.9 |              21.2 |
| ultragsea.c   |            0.040 |                0.9 |               8.9 |
| cor           |            0.005 |                0.0 |               2.7 |
| ztest         |            0.009 |                0.0 |               6.5 |
| goat          |            0.318 |                5.0 |              35.5 |

``` r
rt <- tt[,2]
names(rt) <- tt[,1]
par(mar=c(5,8,2,2))
barplot(sort(rt,decreasing=TRUE), horiz=TRUE, las=1,
  xlab="runtime (sec)")
```

![](compare-methods_files/figure-html/unnamed-chunk-5-1.png)

## Comparing the scores

We can compare the scores between the methods. We see that all scores
are very correlated to each other. CameraPR does not return a score
value, so for the score we computed `-log(p)*sign` as its score.

``` r
res.camera$score <- -log(res.camera$PValue)*(-1+2*(res.camera$Direction=="Up"))
Z <- cbind(
  fgsea = res.fgsea$NES,
  cameraPR = res.camera$score,
  ultragsea.ztest = res.ultragsea.z$score,
  ultragsea.c = res.ultragsea.c$score,
  goat = res.goat$score
)
pairs(Z, pch='.',cex=4)
```

![](compare-methods_files/figure-html/unnamed-chunk-6-1.png)

Notice the gap and the particular S-curve that fgsea shows. Also
cameraPR shows some curve compared to the others. We can compare them
better by comparing the ranks of the scores:

``` r
pairs(apply(Z,2,rank), pch='.',cex=4)
```

![](compare-methods_files/figure-html/unnamed-chunk-7-1.png)

## Comparing the p-values

Next, we can compare the p-values between the methods. Also here, we see
that all p-values are highly correlated to each other.

``` r
P <- cbind(
  fgsea = res.fgsea$pval,
  cameraPR = res.camera$PValue,
  ultragsea.ztest = res.ultragsea.z$pval,
  ultragsea.c = res.ultragsea.c$pval,
  goat = res.goat$pval
)
mlp <- -log10(P)
pairs(mlp, pch='.',cex=4)
```

![](compare-methods_files/figure-html/unnamed-chunk-8-1.png)
