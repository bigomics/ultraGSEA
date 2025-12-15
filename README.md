# UltraGSEA: Ultrafast preranked gene set enrichment analysis (WIP) <a href='https://bigomics.github.io/ultragsea'><img src='man/figures/logo.png' align="right" height="138"/></a>

[![codecov](https://codecov.io/github/bigomics/ultragsea/graph/badge.svg?token=66J6W41C0G)](https://codecov.io/github/bigomics/ultragsea)

(Warning: ***work in progress***! do not use for production purposes)

[ultragsea](https://bigomics.github.io/ultragsea) is an ultrafast
method to compute gene set enrichment on a preranked list of genes not unlike fGSEA. ultragsea can be used as replacement of fGSEA. Although ultragsea uses a different statistical
test (namely z-test or correlation), its scores are highly correlated and its p-values are closely similar to those from GSEA's weighted
Kolmogorov-Smirnov test. ultragsea computes its scores using fast sparse computation and parametric p-values and typically is
100-1000x faster than fGSEA.


## Installation

You can install ultragsea from from GitHub:

```r
remotes::install_github("bigomics/ultragsea")
```

## Usage

For detailed usage examples and tutorials, please see our vignettes:

- [Getting Started with ultragsea](https://bigomics.github.io/ultragsea/articles/ultragsea.html)
- [Comparing ultragsea with other methods](https://bigomics.github.io/ultragsea/articles/compare-methods.html)
- [Benchmarking results](https://bigomics.github.io/ultragsea/articles/benchmark.html)

ultragsea is the main gene set scoring algorithm in OmicsPlayground,
our Bioinformatics platform at BigOmics Analytics. In OmicsPlayground,
you can perform ultragsea without coding needs.

## Example

```{r}
library(ultragsea)
gs <- msigdbr::msigdbr(collection = "H")
gmt <- tapply(gs$gene_symbol,gs$gs_name,list)
G <- gmt2mat(gmt)
fc <- rnorm(nrow(G))
names(fc) <- rownames(G)
res <- ultragsea(fc, G, format='as.gsea', method='ztest')
head(res)
```

## References

For more technical details please refer to our papers. Please cite us when you use
ultragsea as part of your research.

- ultragsea: Ultrafast preranked gene set enrichment scoring. 
- Akhmedov M., et al., Omics Playground: a comprehensive self-service platform for visualization, analytics and exploration of Big Omics Data, NAR Genomics and Bioinformatics, 2020, [lqz019](https://doi.org/10.1093/nargab/lqz019).

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics: help@bigomics.ch

If you like ultragsea, please recommend us to your friends, buy us [coffee](https://buymeacoffee.com/bigomics)
and brag about ultragsea on your social media. 
