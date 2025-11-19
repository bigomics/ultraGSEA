# Fast one sample z-test for matrix object F (e.g. foldchanges) and grouping matrix G (e.g. gene sets).

Fast one sample z-test for matrix object F (e.g. foldchanges) and
grouping matrix G (e.g. gene sets).

## Usage

``` r
fc_ztest(fc, G, zmat = FALSE, alpha = 0.5)
```

## Arguments

- fc:

  Numeric vector of (log)foldchanges, with genes as names

- G:

  Numeric matrix of gene sets, with genes as row names

- zmat:

  Logical indicating whether to compute z-score matrix

- alpha:

  Numeric value between 0 and 1 for Bayesian shrinkage

## Value

List with components:

- z_statistic - Vector of z-test statistics

- p_value - Vector of p-values

- zmat - Z-score matrix (if zmat = TRUE)
