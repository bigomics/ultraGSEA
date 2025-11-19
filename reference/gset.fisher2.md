# Perform Fisher's exact test on gene sets

This function performs Fisher's exact test on two sets of genes,
`genes.up` and `genes.dn`, within a given set of gene sets (`genesets`).
It returns a data frame containing the results of the test, including
the sign of the fold change (positive or negative) and relevant
statistics such as p-values, q-values, and odds ratios.

## Usage

``` r
gset.fisher2(
  genes.up,
  genes.dn,
  genesets,
  background = NULL,
  fdr = 0.05,
  mc = TRUE,
  sort.by = "zratio",
  nmin = 3,
  verbose = 1,
  min.genes = 15,
  max.genes = 500,
  method = "fast.fisher",
  check.background = TRUE,
  common.genes = TRUE
)
```

## Arguments

- genes.up:

  A character vector containing the names of genes in the "up" set.

- genes.dn:

  A character vector containing the names of genes in the "down" set.

- genesets:

  A list of gene sets, where each element is a character vector
  representing a gene set.

- background:

  A character vector containing the names of genes in the background
  set. Defaults to `NULL`, which means all genes are considered.

- fdr:

  The false discovery rate (FDR) threshold for multiple testing
  adjustment. Defaults to 0.05.

- mc:

  A logical value indicating whether to perform multiple testing
  adjustment using Monte Carlo simulation. Defaults to `TRUE`.

- sort.by:

  The statistic used to sort the results. Defaults to "zratio".

- nmin:

  The minimum number of genes required in a gene set for the test to be
  performed.

- verbose:

  A numeric value indicating the level of verbosity. Defaults to 1.

- min.genes:

  The minimum number of genes in a gene set to be considered. Defaults
  to 15.

- max.genes:

  The maximum number of genes in a gene set to be considered. Defaults
  to 500.

- method:

  The method used for computing p-values. Defaults to "fast.fisher".

- check.background:

  A logical value indicating whether to check the presence of genes in
  the background set. Defaults to `TRUE`.

- common.genes:

  A logical value indicating whether to use only genes common to both
  the input gene sets and the background set. Defaults to `TRUE`.

## Value

A data frame containing the results of the Fisher's exact test. The data
frame includes columns such as "sign" (fold change sign), "odd.ratio"
(odds ratio), "p.value" (p-value), "q.value" (adjusted p-value), and
others.
