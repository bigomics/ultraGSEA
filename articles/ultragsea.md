# Getting Started with ultragsea

## Introduction

ultragsea is a novel, ultrafast and memory optimized gene set enrichment
scoring algorithm not unlike fGSEA. ultragsea demonstrates highly
correlated gene set scores and p-values and typically is 100 times
faster than fGSEA.

## Example

### Preparing gensets

For this vignette, we retrieve the Hallmark gene sets using the
`msigdbr` package. Of course you can use your own gene sets or download
from somewhere else.

``` r
library(ultragsea)
library(msigdbr)
gs <- msigdbr::msigdbr(collection = "H")
gmt <- tapply(gs$gene_symbol,gs$gs_name,list)
length(gmt)
#> [1] 50
```

### Performing the statistical test

Just for this example we use some random values for `fc` but in general
this is a fold change from some statistical comparison of your
experiment.

``` r
G <- gmt2mat(gmt)
fc <- rnorm(nrow(G))
names(fc) <- rownames(G)
res <- ultragsea(G, fc, format='as.gsea', method='cor')
res <- ultragsea(G, fc, format='as.gsea', method='ztest')
head(res)
#>                               pathway      pval      padj log2err          ES
#>                                <char>     <num>     <num>   <num>       <num>
#> 1:           HALLMARK_APICAL_JUNCTION 0.8346855 0.9696221      NA  0.08236827
#> 2:           HALLMARK_HEME_METABOLISM 0.8726599 0.9696221      NA  0.06325971
#> 3:     HALLMARK_INFLAMMATORY_RESPONSE 0.7118918 0.9696221      NA  0.14576208
#> 4:         HALLMARK_KRAS_SIGNALING_UP 0.9757536 0.9757536      NA -0.01199556
#> 5: HALLMARK_OXIDATIVE_PHOSPHORYLATION 0.8677728 0.9696221      NA -0.06570964
#> 6:              HALLMARK_ADIPOGENESIS 0.2364614 0.7389418      NA  0.46725361
#>            NES  size                             leadingEdge
#>          <num> <int>                                  <list>
#> 1:  0.19461526   200 CRAT,HRAS,INSIG1,PIK3R3,ACTN2,ADAM9,...
#> 2:  0.14946660   200    CA2,BNIP3L,CDC27,H1-0,BCAM,BLVRA,...
#> 3:  0.34439870   200    IL4R,CCL2,HBEGF,INHBA,ITGB3,EREG,...
#> 4: -0.02834246   200        BMP2,CCND2,LIF,CFB,CSF2,PLAT,...
#> 5: -0.15525515   200 ACAA1,SDHC,ECHS1,GPI,GRPEL1,ATP6V1F,...
#> 6:  1.10400136   200  ALDOA,ECH1,GPX4,CRAT,DHCR7,GADD45A,...
```

For the fastest speed, especially when computing enrichment for multiple
foldchanges, it is better to call the low-level functions for geneset
correlation or geneset z-test directly. They do not return the nicely
formatted output tables as
[`ultragsea()`](https://bigomics.github.io/ultragsea/reference/ultragsea.md)
but they are much faster as they are optimized for matrix inputs `F`.

``` r
F <- cbind(fc,fc,fc,fc,fc)
res1 <- gset.cor(G, F, compute.p=TRUE)
res2 <- gset.ztest(G, F)
```

### Plotting

You can proceed by plotting the gene set using e.g. the plotting
functions of the fgsea package.

``` r
#library(fgsea)
#k <- res$pathway[which.min(res$pval)]
#fgsea::plotEnrichment(gmt[[k]],fc,ticksSize = 0.2)
```

## Session info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] msigdbr_25.1.1   ultragsea_0.1.12 BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-4        babelgene_22.9      jsonlite_2.0.0     
#>  [4] dplyr_1.1.4         compiler_4.5.2      BiocManager_1.30.27
#>  [7] tidyselect_1.2.1    slam_0.1-55         parallel_4.5.2     
#> [10] assertthat_0.2.1    jquerylib_0.1.4     systemfonts_1.3.1  
#> [13] textshaping_1.0.4   yaml_2.3.12         fastmap_1.2.0      
#> [16] lattice_0.22-7      R6_2.6.1            generics_0.1.4     
#> [19] qlcMatrix_0.9.9     curl_7.0.0          knitr_1.50         
#> [22] tibble_3.3.0        bookdown_0.46       desc_1.4.3         
#> [25] bslib_0.9.0         pillar_1.11.1       rlang_1.1.6        
#> [28] cachem_1.1.0        xfun_0.55           fs_1.6.6           
#> [31] sass_0.4.10         cli_3.6.5           pkgdown_2.2.0      
#> [34] withr_3.0.2         magrittr_2.0.4      digest_0.6.39      
#> [37] grid_4.5.2          docopt_0.7.2        lifecycle_1.0.4    
#> [40] vctrs_0.6.5         sparsesvd_0.2-3     data.table_1.17.8  
#> [43] evaluate_1.0.5      glue_1.8.0          ragg_1.5.0         
#> [46] rmarkdown_2.30      matrixStats_1.5.0   tools_4.5.2        
#> [49] pkgconfig_2.0.3     htmltools_0.5.9
```
