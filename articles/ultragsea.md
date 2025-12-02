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
res <- ultragsea(fc, G, format='as.gsea', method='cor')
res <- ultragsea(fc, G, format='as.gsea', method='ztest')
head(res)
#>                               pathway      pval      padj log2err          ES
#>                                <char>     <num>     <num>   <num>       <num>
#> 1:           HALLMARK_APICAL_JUNCTION 0.8340206 0.9703000      NA  0.08495955
#> 2:           HALLMARK_HEME_METABOLISM 0.8732700 0.9703000      NA  0.06467059
#> 3:     HALLMARK_INFLAMMATORY_RESPONSE 0.7109902 0.9703000      NA  0.15022712
#> 4:         HALLMARK_KRAS_SIGNALING_UP 0.9758927 0.9758927      NA -0.01225189
#> 5: HALLMARK_OXIDATIVE_PHOSPHORYLATION 0.8634175 0.9703000      NA -0.06974639
#> 6:              HALLMARK_ADIPOGENESIS 0.2248032 0.7025100      NA  0.49214802
#>            NES  size                             leadingEdge
#>          <num> <int>                                  <list>
#> 1:  0.19721074   200     EGFR,CRAT,THY1,HRAS,INSIG1,MMP2,...
#> 2:  0.15011538   200       BTG2,CA2,GCLC,SLC2A1,BLVRB,PC,...
#> 3:  0.34871186   200        IL6,BTG2,CSF1,IL1B,IL4R,CCL2,...
#> 4: -0.02843947   200 PLAUR,BMP2,CCND2,CXCL10,LIF,TNFAIP3,...
#> 5: -0.16189748   200     IDH1,RETSAT,ACAA1,DLD,IDH2,SDHC,...
#> 6:  1.14238929   200  ALDOA,ECH1,GPX4,CRAT,DHCR7,GADD45A,...
```

For the fastest speed, especially when computing enrichment for multiple
foldchanges, it is better to call the low-level functions for geneset
correlation or geneset z-test directly. They do not return the nicely
formatted output tables as
[`ultragsea()`](https://bigomics.github.io/ultragsea/reference/ultragsea.md)
but they are much faster as they are optimized for matrix inputs `F`.

``` r
F <- cbind(fc,fc,fc,fc,fc)
res1 <- gset.cor(F, G, compute.p=TRUE)
res2 <- fc_ztest(F, G)
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
#>  [1] Matrix_1.7-4             babelgene_22.9           jsonlite_2.0.0          
#>  [4] dplyr_1.1.4              compiler_4.5.2           BiocManager_1.30.27     
#>  [7] Rcpp_1.1.0               tidyselect_1.2.1         slam_0.1-55             
#> [10] parallel_4.5.2           assertthat_0.2.1         jquerylib_0.1.4         
#> [13] systemfonts_1.3.1        textshaping_1.0.4        yaml_2.3.11             
#> [16] fastmap_1.2.0            lattice_0.22-7           R6_2.6.1                
#> [19] generics_0.1.4           qlcMatrix_0.9.9          curl_7.0.0              
#> [22] knitr_1.50               tibble_3.3.0             bookdown_0.45           
#> [25] MatrixGenerics_1.22.0    desc_1.4.3               bslib_0.9.0             
#> [28] pillar_1.11.1            rlang_1.1.6              cachem_1.1.0            
#> [31] xfun_0.54                fs_1.6.6                 sass_0.4.10             
#> [34] cli_3.6.5                pkgdown_2.2.0            withr_3.0.2             
#> [37] magrittr_2.0.4           digest_0.6.39            grid_4.5.2              
#> [40] sparseMatrixStats_1.22.0 docopt_0.7.2             lifecycle_1.0.4         
#> [43] vctrs_0.6.5              sparsesvd_0.2-3          data.table_1.17.8       
#> [46] evaluate_1.0.5           glue_1.8.0               ragg_1.5.0              
#> [49] rmarkdown_2.30           matrixStats_1.5.0        tools_4.5.2             
#> [52] pkgconfig_2.0.3          htmltools_0.5.8.1
```
