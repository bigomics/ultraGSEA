# Perform Fisher's exact test on gene sets

This function performs Fisher's exact test on a set of genes within a
given set of gene sets. It returns a data frame containing the results
of the test, including p-values, q-values, odds ratios, and gene set
overlaps.

## Usage

``` r
gset.fisher(
  genes,
  genesets,
  background = NULL,
  fdr = 0.05,
  mc = TRUE,
  sort.by = "zratio",
  nmin = 3,
  min.genes = 15,
  max.genes = 500,
  method = "fast.fisher",
  check.background = TRUE,
  common.genes = TRUE,
  no.pass = NA,
  verbose = 1
)
```

## Arguments

- genes:

  A character vector containing the names of genes.

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
  the input gene set and the background set. Defaults to `TRUE`.

- verbose:

  A numeric value indicating the level of verbosity. Defaults to 1.

## Value

A data frame containing the results of the Fisher's exact test. The data
frame includes columns such as "p.value" (p-value), "q.value" (adjusted
p-value), "odd.ratio" (odds ratio), "overlap" (gene set overlap), and
optionally "genes" (common genes).
