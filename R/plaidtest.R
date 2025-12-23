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
  G1 <- gsetX[,which(y1!=ref),drop=FALSE]
  G0 <- gsetX[,which(y1==ref),drop=FALSE]
  res <- matrixTests::row_t_welch(G1, G0)
  df <- data.frame(
    pathway = rownames(res),
    logFC = res[,"mean.diff"],
    pval = res[,"pvalue"],
    padj = p.adjust(res[,"pvalue"], method="BH")
  )
  rownames(df) <- rownames(res)
  df
}

plaid.cortest <- function(X, G, y, ref, nsmooth=3) {
  if(length(unique(y[!is.na(y)]))>2) message("[plaid.test] warning: more than 2 classes")
  if(length(unique(y[!is.na(y)]))<2) stop("[plaid.test] warning: less than 2 classes")
  sel <- which(!is.na(y))
  y1 <- 1*(y[sel]!=ref)
  X1 <- X[,sel,drop=FALSE]
  gsetX <- plaid::plaid(X1, G, nsmooth=nsmooth)
  rho <- cor(Matrix::t(gsetX), y1)[,1]
  pv <- cor_pvalue(rho, ncol(gsetX))
  df <- data.frame(
    pathway = rownames(gsetX),
    rho = rho,
    pval = pv,
    padj = p.adjust(pv, method="BH")
  )
  rownames(df) <- rownames(gsetX)  
  df
}

plaid.dualtest <- function(X, G, y, ref=0, test1="ttest", pbalance=TRUE) {
  if(!ref %in% y) stop("need to specify ref")
  if(test1 == "ttest") {
    t1 <- plaid.ttest(X=X, G=G, y=y, ref=ref)
    t1 <- t1[colnames(G),]
    stat1 <- t1$logFC    
  } else if(test1 == "cortest") {
    t1 <- plaid.cortest(X=X, G=G, y=y, ref=ref)
    t1 <- t1[colnames(G),]
    stat1 <- t1$rho
  } else {
    stop("Stop invalid method: test1 = ",test1)
  }
  fc <- rowMeans(X[,which(y!=ref)],na.rm=TRUE) -
    rowMeans(X[,which(y==ref)],na.rm=TRUE)
  t2 <- gset.cortest(gset=G, fc, compute.p = TRUE,
    use.rank = FALSE, corshrink = 3)
  t2 <- do.call(cbind, lapply(t2, function(x) x[,1]))
  t2 <- t2[colnames(G),]
  pp <- cbind(t1$pval, t2[,"p.value"])
  if(pbalance) {
    ## balance p-values by scaling
    q99 <- apply(-log10(pp),2,quantile,probs=0.99,na.rm=TRUE)
    q99 <- q99 / mean(q99,na.rm=TRUE)
    pp <- 10^(-sweep(-log10(pp), 2, q99, '/'))
  }
  metap <- matrix_metap(pp)
  df <- data.frame(
    pathway = colnames(G),
    pval = metap,
    padj = p.adjust(metap),
    stat1 = stat1,
    pval1 = t1$pval,
    stat2 = t2[,"rho"],
    pval2 = t2[,"p.value"]
  )
  rownames(df) <- colnames(G)  
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
