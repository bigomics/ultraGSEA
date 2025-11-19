# Ultra-fast GSEA using z-score or geneset correlation. Results of these methods highly correlate with GSEA/fGSEA but are much faster.

Ultra-fast GSEA using z-score or geneset correlation. Results of these
methods highly correlate with GSEA/fGSEA but are much faster.

## Usage

``` r
ultragsea(
  fc,
  G,
  alpha = 0.5,
  minLE = 1,
  method = c("ztest", "ttest", "cor", "rankcor")[1],
  format = c("simple", "as.gsea")[1]
)
```
