# Compute z-score matrix for gene sets

Compute z-score matrix for gene sets

## Usage

``` r
fc_zmat(fc, G, alpha = 0.5)
```

## Arguments

- fc:

  Numeric vector of (log)foldchanges, with genes as names

- G:

  Numeric matrix of gene sets, with genes as row names

- alpha:

  Numeric value between 0 and 1 for Bayesian shrinkage

## Value

Z-score matrix
