# UltraGSEA: Ultrafast preranked gene set enrichment analysis (WIP)

[![codecov](https://codecov.io/github/bigomics/zgsea/graph/badge.svg?token=66J6W41C0G)](https://codecov.io/github/bigomics/ultragsea)

(Warning: ***work in progress***! do not use for production purposes)

[ultragsea](https://bigomics.github.io/ultragsea) is an ultrafast method
to compute gene set enrichment on a preranked list of genes like GSEA
and fGSEA. Although ultragsea uses a different statistical test (namely
z-test and correlation), its scores are highly correlated and its
p-values are highly similar to those from GSEA’s weighted
Kolmogorov-Smirnov test. ultragsea computes its scores using fast sparse
computation and parametric p-values and typically is more than 100x
faster than fGSEA.

#### Key features

- Ultra-fast preranked gene set enrichment scoring
- Works with regular matrices, sparse matrices, and Bioconductor data
  structures
- Automatically detects and handles Bioconductor objects
  ([`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html),
  [`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html),
  [`BiocSet`](https://bioconductor.org/packages/release/bioc/html/BiocSet.html))
- Built-in differential enrichment testing

#### Warning

ultragsea is fast. Ludicrously fast. Please fasten your seatbelts before
usage.

## Installation

You can install ultragsea from from GitHub:

``` r
remotes::install_github("bigomics/ultragsea")
```

## Usage

For detailed usage examples and tutorials, please see our vignettes:

- [Getting Started with
  ultragsea](https://bigomics.github.io/ultragsea/articles/01_ultragsea-vignette.html)
- [Comparing ultragsea with
  fGSEA](https://bigomics.github.io/ultragsea/articles/02_compare-vignette.html)

ultragsea is the main gene set scoring algorithm in OmicsPlayground, our
Bioinformatics platform at BigOmics Analytics. In OmicsPlayground, you
can perform ultragsea without coding needs.

## Example

`{r} library("ultragsea") gs <- msigdbr::msigdbr(collection = "H") gmt <- tapply(gs$gene_symbol,gs$gs_name,list) G <- gmt2mat(gmt) fc <- rnorm(nrow(G)) names(fc) <- rownames(G) res <- ultragsea(fc, G, format='as.gsea', method='ztest') head(res)`

## References

For more technical details please refer to our papers. Please cite us
when you use ultragsea as part of your research.

- ultragsea: Ultrafast preranked gene set enrichment scoring.
- Akhmedov M., et al., Omics Playground: a comprehensive self-service
  platform for visualization, analytics and exploration of Big Omics
  Data, NAR Genomics and Bioinformatics, 2020,
  [lqz019](https://doi.org/10.1093/nargab/lqz019).

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics: <help@bigomics.ch>

If you like ultragsea, please recommend us to your friends, buy us
[coffee](https://buymeacoffee.com/bigomics) and brag about ultragsea on
your social media.
