# Fast one sample t-test for matrix object F (e.g. foldchanges) and grouping matrix G (e.g. gene sets).

Fast one sample t-test for matrix object F (e.g. foldchanges) and
grouping matrix G (e.g. gene sets).

## Usage

``` r
matrix_onesample_ttest(F, G)
```

## Arguments

- F:

  Numeric vector or matrix of (log)foldchanges

- G:

  Numeric matrix of gene sets, with genes as row names

## Value

List with components:

- mean - Mean values

- t - t-test statistics

- p - p-values from t-test
