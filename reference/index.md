# Package index

## All functions

- [`cor_pvalue()`](https://bigomics.github.io/ultragsea/reference/cor_pvalue.md)
  : Calculate p-value from correlation rho
- [`cor_sparse_matrix()`](https://bigomics.github.io/ultragsea/reference/cor_sparse_matrix.md)
  : Calculate sparse correlation matrix handling missing values
- [`corpora.fastFET()`](https://bigomics.github.io/ultragsea/reference/corpora.fastFET.md)
  : Wrapper superfast version of Fisher Exact Test from 'corpora' R
  package. This is the fastest implementation currently available. Uses
  phyper inside.
- [`fgsea()`](https://bigomics.github.io/ultragsea/reference/fgsea.md) :
  Fast replacement of fgsea with p/q values computed with ultragsea.cor
- [`genesetr.fastFET()`](https://bigomics.github.io/ultragsea/reference/genesetr.fastFET.md)
  : Fast version of Fisher Exact Test from 'genesetr' R package. Quickly
  compute FET p-values for n 2x2 contingency tables T =
  {t_1,t_2,...,t_n-1,t_n}
- [`gmt2mat()`](https://bigomics.github.io/ultragsea/reference/gmt2mat.md)
  : Convert GMT to Binary Matrix
- [`goat()`](https://bigomics.github.io/ultragsea/reference/goat.md) :
  Wrapper around GOAT
- [`gset.cor()`](https://bigomics.github.io/ultragsea/reference/gset.cor.md)
  : Calculate gene set correlation
- [`gset.fastFET()`](https://bigomics.github.io/ultragsea/reference/gset.fastFET.md)
  : Calculate fast Fisher exact test.
- [`gset.ztest()`](https://bigomics.github.io/ultragsea/reference/gset.ztest.md)
  : Fast one sample z-test for matrix object F (e.g. foldchanges) and
  grouping matrix G (e.g. gene sets).
- [`makeGseaResult()`](https://bigomics.github.io/ultragsea/reference/makeGseaResult.md)
  : Utility to convert fgsea/ultragsea style results to 'gseaResult'
  object compatible with clusterProfiler and enrichplot functions. Using
  this function you can use all plotting functions from enrichplot.
- [`mat2gmt()`](https://bigomics.github.io/ultragsea/reference/mat2gmt.md)
  : Convert Binary Matrix to GMT
- [`matrix_metap()`](https://bigomics.github.io/ultragsea/reference/matrix_metap.md)
  : Matrix version for combining p-values using fisher or stouffer
  method. Much faster than doing metap::sumlog() and metap::sumz()
- [`read.gmt()`](https://bigomics.github.io/ultragsea/reference/read.gmt.md)
  : Read GMT File
- [`ultragsea()`](https://bigomics.github.io/ultragsea/reference/ultragsea.md)
  : Ultra-fast GSEA using z-score or geneset correlation. Results of
  these methods highly correlate with GSEA/fGSEA but are much faster.
- [`write.gmt()`](https://bigomics.github.io/ultragsea/reference/write.gmt.md)
  : Write GMT File
