# Calculate fast Fisher exact test.

Calculate fast Fisher exact test.

## Usage

``` r
gset.fastFET(
  genes,
  G,
  bg,
  report.genes = FALSE,
  method = c("phyper", "genesetr")[1]
)
```

## Arguments

- genes:

  Vector of significant genes

- G:

  Sparse matrix containing gene sets
