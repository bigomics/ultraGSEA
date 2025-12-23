# Ultra-fast GSEA using z-score or geneset correlation. Results of these methods highly correlate with GSEA/fGSEA but are much faster.

Ultra-fast GSEA using z-score or geneset correlation. Results of these
methods highly correlate with GSEA/fGSEA but are much faster.

## Usage

``` r
ultragsea(
  G,
  fc,
  alpha = 0,
  minLE = 1,
  corshrink = 3,
  minsize = 1L,
  maxsize = 9999L,
  center = TRUE,
  method = c("cor", "ztest", "ttest", "goat", "camera")[1],
  format = c("simple", "as.gsea")[1]
)
```
