# Fast one sample z-test for matrix object F (e.g. foldchanges) and grouping matrix G (e.g. gene sets).

Fast one sample z-test for matrix object F (e.g. foldchanges) and
grouping matrix G (e.g. gene sets).

## Usage

``` r
gset.ztest(F, G, alpha = 0.5, center = TRUE, pdist = "norm")
```

## Arguments

- G:

  Numeric matrix of gene sets, with genes as row names

- alpha:

  Numeric value between 0 and 1 for Bayesian shrinkage

- center:

  Logical value to center population mean or assume zero

- fc:

  Numeric vector of (log)foldchanges, with genes as names

## Value

List with components:

- z_statistic - Vector of z-test statistics

- p_value - Vector of p-values
