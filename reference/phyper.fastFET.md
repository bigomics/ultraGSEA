# Wrapper superfast version of Fisher Exact Test. This is the fastest implementation currently available. Uses phyper inside. Original code from 'corpora' R package

           setAn ¬setA
       setB  a     b | a+b
      ¬setB  c     d | c+d
         ------------|-----
            a+c   b+d| a+b+c+d

## Usage

``` r
phyper.fastFET(
  a,
  b,
  c,
  d,
  alternative = c("two.sided", "less", "greater")[3],
  log.p = FALSE
)
```

## Details

Example: Given three contingency tables: setA ¬setA setA ¬setA setB
¬setB setB 3 100 setC 5 6 setD 20 45 ¬setB 123 500 ¬setC 10 100 ¬setD 60
1000

a=c(3, 5, 20) b=c(100, 6, 45) c=c(123, 10, 60) d=c(500, 100, 1000) pvals
= fastFET(a,b,c,d)
