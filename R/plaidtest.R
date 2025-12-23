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

plaid.dualtest <- function(X, G, y, ref) {
  test1 <- plaid.ttest(X, G, y, ref)
  test2 <- gset.cor(G, FC, compute.p = FALSE,
    use.rank = FALSE, corshrink = 0)
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
