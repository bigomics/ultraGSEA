##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

## ========================================================================
## =============== single-sample based enrichment tests ===================
## ========================================================================

plaid.ttest <- function(X, G, y, ref, nsmooth=3) {
  if(length(unique(y[!is.na(y)]))>2) message("[plaid.test] warning: more than 2 classes")
  if(length(unique(y[!is.na(y)]))<2) stop("[plaid.test] warning: less than 2 classes")
  sel <- which(!is.na(y))
  y1 <- y[sel]
  X1 <- X[,sel,drop=FALSE]
  gsetX <- plaid::plaid(X1, G, nsmooth=nsmooth)
  res <- matrixTests::row_t_welch(gsetX[,which(y1!=ref),drop=FALSE],
    gsetX[,which(y1==ref),drop=FALSE] )
  res
}

plaid.dualtest <- function(X, G, y, ref=0) {
  if(!ref %in% y) stop("need to specify ref")
  t1 <- plaid.ttest(X=X, G=G, y=y, ref=ref)
  fc <- rowMeans(X[,which(y==1)],na.rm=TRUE) -
    rowMeans(X[,which(y==0)],na.rm=TRUE)
  t2 <- gset.cor(gset=G, fc, compute.p = TRUE,
    use.rank = FALSE, corshrink = 3)
  t2 <- do.call(cbind, lapply(t2, function(x) x[,1]))
  t1 <- t1[colnames(G),]
  t2 <- t2[colnames(G),]    
  metap <- matrix_metap(cbind(t1$pvalue, t2[,"p.value"]))
  df <- data.frame(
    pathway = colnames(G),
    pval = metap,
    padj = p.adjust(metap),
    pval.TT = t1$pvalue,
    pval.COR = t2[,"p.value"],
    rho = t2[,"rho"],
    mean.x = t1$mean.x,
    mean.y = t1$mean.y,
    mean.diff = t1$mean.diff    
  )
  return(df)
}

#' Matrix version for combining p-values using fisher or stouffer
#' method. Much faster than doing metap::sumlog() and metap::sumz()
#'
#' @export
matrix_metap <- function(plist, method = c("stouffer","fisher","maxp")[1]) {
  if (inherits(plist, "matrix")) {
    plist <- as.list(data.frame(plist))
  }
  if (method %in% c("fisher", "sumlog")) {
    chisq <- (-2) * Reduce("+", lapply(plist, log))
    df <- 2 * length(plist)
    pv <- pchisq(chisq, df, lower.tail = FALSE)
  } else if (method %in% c("stouffer", "sumz")) {
    np <- length(plist)
    zz <- lapply(plist, qnorm, lower.tail = FALSE)
    zz <- Reduce("+", zz) / sqrt(np)
    pv <- pnorm(zz, lower.tail = FALSE)
  } else if(method %in% c("maxp","pmax","maximump")) {
    pv <- Reduce(pmax, plist)    
  } else {
    stop("Invalid method: ", method)
  }
  dimnames(pv) <- dimnames(plist[[1]])
  return(pv)
}
