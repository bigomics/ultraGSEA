# Matrix version for combining p-values using fisher or stouffer method. Much faster than doing metap::sumlog() and metap::sumz()

Matrix version for combining p-values using fisher or stouffer method.
Much faster than doing metap::sumlog() and metap::sumz()

## Usage

``` r
matrix_metap(plist, method = c("stouffer", "fisher", "maxp")[1])
```
