#' Calculate fast Fisher exact test.  
#'
#' @param genes Vector of significant genes
#' @param G Sparse matrix containing gene sets
#'
#' @export
gset.fastFET <- function(genes, G, bg, report.genes=FALSE,
                         method=c('phyper','genesetr')[1]) {
  bgnull <- FALSE
  if(is.null(bg)) {
    message("[gset.fisher] note: it is recommended to specify background")
    bg <- unique(c(genes,rownames(G)))
    bgnull <- TRUE
  }
  if(length(bg)>1 && !bgnull) {
    genes <- intersect(genes, bg)
    G <- G[intersect(bg,rownames(G)),,drop=FALSE]
  }  
  length.bg <- NULL
  if(length(bg)==1 && is.integer(bg)) {
    length.bg <- as.integer(bg)
  } else if(length(bg)>1) {
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
  if(method=='phyper') {
    pv <- phyper.fastFET(a,b,c,d,alternative="greater")
  } else if(method=='genesetr') {
    pv <- genesetr::fastFET(a,b,c,d,alternative="greater")
  } else {
    stop("invalid method")
  }
  names(pv) <- colnames(G)
  odd.ratio <- (a/b)/(c/d)
  qv <- p.adjust(pv, method="fdr")
  overlap <- paste0(a,"/",gsize)

  gsgenes <- NULL
  if(report.genes) {
    G1 <- G[sort(genes),]
    idx <- Matrix::which(G1 != 0, arr.ind=TRUE)
    idx[,1] <- rownames(G1)[idx[,1]]
    gsgenes <- tapply( idx[,1], idx[,2], c)
    gsgenes <- sapply(gsgenes, function(s) paste(s,collapse="|"))    
    names(gsgenes) <- colnames(G1)[as.integer(names(gsgenes))]    
    gsgenes <- gsgenes[match(colnames(G),names(gsgenes))]
  }
  
  df <- data.frame(p.value=pv, q.value=qv, odd.ratio=odd.ratio, overlap=overlap)
  if(report.genes) df <- cbind(df, genes=gsgenes)
  return(df)
}


#' Wrapper superfast version of Fisher Exact Test. This is the fastest
#' implementation currently available. Uses phyper inside. Original
#' code from 'corpora' R package
#' 
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
phyper.fastFET <- function(a, b, c, d, alternative = c("two.sided", "less", 
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


