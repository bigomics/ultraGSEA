# ultraGSEA: Ultrafast preranked gene set enrichment analysis (WIP)

[![codecov](https://codecov.io/github/bigomics/zgsea/graph/badge.svg?token=66J6W41C0G)](https://codecov.io/github/bigomics/ultraGSEA)

(***work in progress***! do not use for production purposes)

[ultraGSEA](https://bigomics.github.io/ultraGSEA) is an ultrafast
method to compute gene set enrichment on a preranked list of genes
like GSEA and fGSEA. Although UltraGSEA uses a different statistical
test (namely z-test and correlation), its scores are highly correlated
and its p-values are highly similar to those from GSEA's weighted
Kolmogorov-Smirnov test. UltraGSEA computes its scores using fast
sparse computation and parametric p-values and typically is more than
100x faster than fGSEA.

#### Key features

- Ultra-fast preranked gene set enrichment scoring
- Works with regular matrices, sparse matrices, and Bioconductor data structures
- Automatically detects and handles Bioconductor objects 
([`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html), 
[`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html), [`BiocSet`](https://bioconductor.org/packages/release/bioc/html/BiocSet.html))
- Built-in differential enrichment testing

#### Warning

UltraGSEA is fast. Ludicrously fast. Please fasten your seatbelts
before usage.

## Installation

You can install ultraGSEA from from GitHub:

```r
remotes::install_github("bigomics/ultraGSEA")
```

## Usage

For detailed usage examples and tutorials, please see our vignettes:

- [Getting Started with ultraGSEA](https://bigomics.github.io/ultraGSEA/articles/01_ultraGSEA-vignette.html)
- [Comparing ultraGSEA with fGSEA](https://bigomics.github.io/ultraGSEA/articles/02_compare-vignette.html)

UltraGSEA is the main gene set scoring algorithm in OmicsPlayground,
our Bioinformatics platform at BigOmics Analytics. In OmicsPlayground,
you can perform ultraGSEA without coding needs.

## Example

```{r}
library("ultraGSEA")
gs <- msigdbr::msigdbr(collection = "H")
gmt <- tapply(gs$gene_symbol,gs$gs_name,list)
G <- gmt2mat(gmt)
fc <- rnorm(nrow(G))
names(fc) <- rownames(G)
res <- ultraGSEA(fc, G, format='as.gsea', method='ztest')
head(res)
```

## References

For more technical details please refer to our papers. Please cite us when you use
ultraGSEA as part of your research. 

- UltraGSEA: Ultrafast preranked gene set enrichment scoring. 
- Akhmedov M., et al., Omics Playground: a comprehensive self-service platform for visualization, analytics and exploration of Big Omics Data, NAR Genomics and Bioinformatics, 2020, [lqz019](https://doi.org/10.1093/nargab/lqz019).

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics: help@bigomics.ch

If you like ultraGSEA, please recommend us to your friends, buy us [coffee](https://buymeacoffee.com/bigomics)
and brag about ultraGSEA on your social media. 
