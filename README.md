# ultraGSEA: Ultrafast preranked gene set enrichment analysis  <img src='man/figures/logo.png' align="right" height="150"/>

[![codecov](https://codecov.io/github/bigomics/zgsea/graph/badge.svg?token=66J6W41C0G)](https://codecov.io/github/bigomics/ultraGSEA)

[ultraGSEA](https://bigomics.github.io/ultraGSEA) is an ultrafast
method to compute gene set enrichment on a preranked list of genes not
unlike GSEA and its faster cousin fGSEA. UltraGSEA computes the
enrichment scores using fast sparse computation and parametric
p-values and typically is more than 100 times faster than fGSEA.

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
