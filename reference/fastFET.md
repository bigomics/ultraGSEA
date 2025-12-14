# Fast version of Fisher Exact Test from 'genesetr' R package.

Original code from https://github.com/MaayanLab/genesetr

## Usage

``` r
fastFET(a, b, c, d, alternative = "greater")
```

## Arguments

- a:

  A vector of values where a_i corresponds to t_i

- b:

  A vector of values where b_i corresponds to t_i

- c:

  A vector of values where c_i corresponds to t_i

- d:

  A vector of values where d_i corresponds to t_i

## Value

A vector of p-values P=p_1, p_2,...,p_n-1,p_n where p_i corresponds to
the FET result of t_i.

## Details

Quickly compute FET p-values for n 2x2 contingency tables T =
{t_1,t_2,...,t_n-1,t_n}

## Examples
