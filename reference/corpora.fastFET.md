# Wrapper superfast version of Fisher Exact Test from 'corpora' R package. This is the fastest implementation currently available. Uses phyper inside.

           setAn ¬setA
       setB  a     b | a+b
      ¬setB  c     d | c+d
         ------------|-----
            a+c   b+d| a+b+c+d

## Usage

``` r
corpora.fastFET(
  a,
  b,
  c,
  d,
  alternative = c("two.sided", "less", "greater")[3],
  log.p = FALSE
)
```
