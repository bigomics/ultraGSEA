# Read GMT File

Read data from a GMT file (Gene Matrix Transposed). The GMT format is
commonly used to store gene sets or gene annotations.

## Usage

``` r
read.gmt(gmt.file, dir = NULL, add.source = FALSE, nrows = -1)
```

## Arguments

- gmt.file:

  Path to GMT file.

- dir:

  (Optional) The directory where the GMT file is located.

- add.source:

  (optional) Include the source information in the gene sets' names.

- nrows:

  (optional) Number of rows to read from the GMT file.

## Value

A list of gene sets: each gene set is represented as a character vector
of gene names.

## Examples

``` r
# \donttest{
# Read GMT file (requires file to exist)
gmt_file <- system.file("extdata", "hallmarks.gmt", package = "plaid")
if (file.exists(gmt_file)) {
  gmt <- read.gmt(gmt_file)
  print(names(gmt))
  print(head(gmt[[1]]))
  
  # Read with source information
  gmt_with_source <- read.gmt(gmt_file, add.source = TRUE)
  print(head(names(gmt_with_source)))
}
# }
```
