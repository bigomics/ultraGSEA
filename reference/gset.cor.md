# Calculate gene set correlation

Compute correlation between a foldchange vector/matrix and gene sets

## Usage

``` r
gset.cor(gset, FC, compute.p = FALSE, use.rank = FALSE, corshrink = 0)
```

## Arguments

- gset:

  Numeric matrix of gene sets, with genes as row/column names

- FC:

  Numeric vector or matrix of (log)foldchanges, with genes as row names

- compute.p:

  Logical indicating whether to compute p-values

- use.rank:

  Logical indicating whether to rank transform FC before correlation

## Value

Named list with components:

- rho - Matrix of correlation coefficients between FC and gset

- p.value - Matrix of p-values for correlation (if compute.p = TRUE)

- q.value - Matrix of FDR adjusted p-values (if compute.p = TRUE)

## Details

This function calculates sparse rank correlation between FC and each
column of gset using
[`qlcMatrix::corSparse()`](https://rdrr.io/pkg/qlcMatrix/man/corSparse.html).
It handles missing values in FC by computing column-wise correlations.

P-values are computed from statistical distribution

## Examples

``` r
if (FALSE) { # \dontrun{
library(playbase)
ranks <- sample(1:10000, 1000, replace = TRUE)
names(ranks) <- replicate(1000, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
genesets <- matrix(rnorm(1000 * 20), ncol = 20)
rownames(genesets) <- names(ranks)

gset.rankcor(ranks, genesets, compute.p = TRUE)
} # }
```
