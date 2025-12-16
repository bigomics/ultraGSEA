#' Calculate fast Fisher exact test.  
#'
#' @param genes Vector of significant genes
#' @param G Sparse matrix containing gene sets
#'
#' @export
gset.fastFET <- function(genes, G, bg, method=2) {
  if(is.null(bg)) bg <- unique(c(genes,rownames(G)))
  length.bg <- NULL
  if(length(bg)==1 && is.integer(bg)) {
    length.bg <- as.integer(bg)
  } else if(length(bg)>1) {
    genes <- intersect(genes, bg)
    G <- G[intersect(bg,rownames(G)),,drop=FALSE]
    length.bg <- length(bg)
  }

  if(is.null(length.bg)) stop("error: invalid background. bg:", head(bg))
  if(length(genes)==0) stop("error: zero genes length")
  if(nrow(G)==0) stop("error: empty gene set matrix G")
  
  genes <- intersect(genes, rownames(G))  
  gsize <- Matrix::colSums(G!=0)
  genes <- intersect(genes, rownames(G))
  a <- Matrix::colSums(G[genes,]!=0)
  b <- length(genes) - a
  c <- gsize - a
  d <- length.bg - (a+b+c) 
  if(method==1) {
    pv <- genesetr.fastFET(a,b,c,d)
  } else if(method==2) {
    pv <- corpora.fastFET(a,b,c,d)
  } else {
    stop("[gset.fastFET] unknown method")
  }

  names(pv) <- colnames(G)
  odd.ratio <- (a/b)/(c/d)
  qv <- p.adjust(pv, method="fdr")
  overlap <- paste0(a,"/",gsize)
  data.frame(p.value=pv, q.value=qv, odd.ratio=odd.ratio, overlap=overlap)
}

#' Fast version of Fisher Exact Test from 'genesetr' R
#' package. Quickly compute FET p-values for n 2x2 contingency tables
#' T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}
#' 
#' Original code from https://github.com/MaayanLab/genesetr
#'
#'
#' @param a A vector of values where a_{i} corresponds to t_{i}
#' @param b A vector of values where b_{i} corresponds to t_{i}
#' @param c A vector of values where c_{i} corresponds to t_{i}
#' @param d A vector of values where d_{i} corresponds to t_{i}
#' @return A vector of p-values P={p_{1}, p_{2},...,p_{n-1},p_{n}} where p_i corresponds to the FET result of t_i.
#' @examples
#' T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}, with each t_i of the form:
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#'
#' Example: Given three contingency tables:
#'      setA  ¬setA          setA  ¬setA         setB  ¬setB
#' setB  3     100        setC  5     6        setD  20     45
#'¬setB  123   500       ¬setC  10    100     ¬setD  60     1000
#'
#' a=c(3, 5, 20)
#' b=c(100, 6, 45)
#' c=c(123, 10, 60)
#' d=c(500, 100, 1000)
#' pvals = fastFET(a,b,c,d)
#' 
genesetr.fastFET = function(a, b, c, d, alternative = "greater") {

  checkContingTableVals = function(a,b,c,d){
    l_a = length(a)
    l_b = length(b)
    l_c = length(c)
    l_d = length(d)
    
    if(any(!rep(l_a,3)==c(l_b,l_c,l_d))){
      stop("Vectors a,b,c and d must all be of equal length.")
    }
    if(!is.numeric(c(a, b, c, d))){
      stop("All elements of a,b,c, and d must be numeric.")
    }
    if(any(is.na(c(a,b,c,d)))){
      warning("NA elements will produce NA result.")
    }
  }

  checkContingTableVals(a,b,c,d)

  l_a = length(a)
  #if a,b,c,d are length 1 use FET from standard R distribution
  if(l_a == 1) return(fisher.test(matrix(c(a,c,b,d),ncol = 2),
    alternative = alternative)$p.value)

  #populate factorial look-up vector
  f = numeric(length = max(rowSums(cbind(a,b,c,d))))
  f[1] = log(1)
  for(i in 2:length(f)){
    f[i] = f[i-1] + log(i)
  }

  pval = cbind(numeric(length = l_a),numeric(length = l_a))
  coeff = 1
  min = cbind(pmin(b,c),rep(NA,l_a))
  sides = 1

  if(alternative == "less"){
    coeff = -1
    min = cbind(pmin(a,d),rep(NA,l_a))
  }

  if(alternative == "two.sided"){
    coeff = c(-1,1)
    min = cbind(pmin(a,d),pmin(b,c))
    sides = 2
  }
for(j in 1:sides){
  for(i in 1:l_a){
    temp_a = a[i]:(a[i]+coeff[j]*min[i,j])
    temp_b = b[i]:(b[i]-coeff[j]*min[i,j])
    temp_c = c[i]:(c[i]-coeff[j]*min[i,j])
    temp_d = d[i]:(d[i]+coeff[j]*min[i,j])

    #replace 0s with 1 since R indexing starts at 1
    temp_a[temp_a==0] = 1
    temp_b[temp_b==0] = 1
    temp_c[temp_c==0] = 1
    temp_d[temp_d==0] = 1

    invariant = (f[a[i]+b[i]] + f[c[i]+d[i]] + f[a[i]+c[i]] + f[b[i]+d[i]])-
      f[a[i]+b[i]+c[i]+d[i]]#n
    pval[i,j] = sum(exp(invariant -
        (f[temp_a] + f[temp_b] + f[temp_c] + f[temp_d])))
  }}
  #computed as in fisher.exact in R stats package
  #unclear whether this is appropriate -Ali 5/21/2018
  if(alternative == "two.sided"){
    pval[,1] = pmin(pval[,1],pval[,2])
  }
  return(pval[,1])
}


#' Wrapper superfast version of Fisher Exact Test from 'corpora' R
#' package. This is the fastest implementation currently
#' available. Uses phyper inside.
#' 
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#' 
corpora.fastFET <- function(a, b, c, d, alternative = c("two.sided", "less", 
    "greater")[3], log.p = FALSE) {
  ## this is really-really-really fast...
  pv <- rep(NA, length(a))
  ii <- 1:length(a)
  ii <- which((a + b) > 0)
  d1 <- d + 1 * (d == 0) ## hack to avoid crash...
  c1 <- c + 1 * (c == 0) ## hack to avoid crash...

  k1 <- a[ii]
  n1 <- (a + c1)[ii]
  k2 <- b[ii]
  n2 <- (b + d1)[ii]
  
  .match.len <-  function (vars, len = NULL, adjust = FALSE, check.numeric = TRUE, 
                           envir = parent.frame()) {
    vecs <- setNames(lapply(vars, get, envir = envir), vars)
    ok <- sapply(vecs, is.numeric)
    if (check.numeric && any(!ok)) 
      stop("argument(s) ", paste(vars[!ok], collapse = ", "), 
        " must be numeric vector(s)")
    if (is.null(len)) 
      len <- max(sapply(vecs, length))
    for (v in vars) {
      if (length(vecs[[v]]) == 1) {
        if (adjust) 
          assign(v, rep(vecs[[v]], len), envir = envir)
      }
      else if (length(vecs[[v]]) != len) {
        stop(sprintf("argument %s should be of length %d or a scalar (%s must have same length)", 
          v, len, paste(vars, collapse = ", ")))
      }
    }
    invisible(len)
  }
  
  alternative <- match.arg(alternative)
  l <- .match.len(c("k1", "n1", "k2", "n2"), adjust = TRUE)
  if (any(k1 < 0) || any(k1 > n1) || any(n1 <= 0)) 
    stop("k1 and n1 must be integers with 0 <= k1 <= n1")
  if (any(k2 < 0) || any(k2 > n2) || any(n2 <= 0)) 
    stop("k2 and n2 must be integers with 0 <= k2 <= n2")
  if (any(k1 + k2 <= 0)) 
    stop("either k1 or k2 must be non-zero")
  k <- k1 + k2
  if (alternative == "two.sided") {
    if (log.p) {
      pval <- pmin(phyper(k1 - 1, n1, n2, k, lower.tail = FALSE, 
        log.p = TRUE), phyper(k1, n1, n2, k, lower.tail = TRUE, 
          log.p = TRUE)) + log(2)
      pval <- pmin(pval, 0)
    }
    else {
      pval <- 2 * pmin(phyper(k1 - 1, n1, n2, k, lower.tail = FALSE), 
        phyper(k1, n1, n2, k, lower.tail = TRUE))
      pval <- pmax(0, pmin(1, pval))
    }
  }
  else if (alternative == "greater") {
    pval <- phyper(k1 - 1, n1, n2, k, lower.tail = FALSE, 
      log.p = log.p)
  }
  else if (alternative == "less") {
    pval <- phyper(k1, n1, n2, k, lower.tail = TRUE, log.p = log.p)
  }
  pval
}


