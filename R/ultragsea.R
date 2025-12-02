
##---------------------------------------------------------------
##-------------------- ultragsea --------------------------------
##---------------------------------------------------------------

#' Ultra-fast GSEA using z-score or geneset correlation. Results of
#' these methods highly correlate with GSEA/fGSEA but are much faster.
#'
#' @export
ultragsea <- function(fc, G, alpha=0.5, minLE=1, cshrink=3,
                      minsize = 1L, maxsize = 9999L,
                      method=c("cor","ztest","ttest","goat")[1],
                      format=c("simple","as.gsea","as.gsea2")[1]) {

  gg <- intersect(names(fc), rownames(G))
  fc <- fc[gg]
  G <- G[gg,]
  size <- Matrix::colSums(G!=0)
  if(method == "goat") minsize <- max(minsize,10L)
  sel <- which(size >= minsize & size <= maxsize)
  G <- G[, sel, drop=FALSE]
  size <- size[sel]
  
  addLE <- (format=="as.gsea")
  p_value <- NULL
  q_value <- NULL
  stat_value <- NULL
  zmat <- NULL
  gmt <- NULL
  if(method == "ztest") {
    zres <- fc_ztest(fc, G, zmat=addLE, alpha=alpha)
    p_value <- zres$p
    stat_value <- zres$z
    zmat <- zres$zmat
  } else if(method == "ttest") {
    res <- matrix_onesample_ttest(cbind(fc), G)
    p_value <- res$p[,1]    
    stat_value <- res$t[,1]    
  } else if(method == "cor") {
    res <- gset.cor(fc, G, compute.p=TRUE, use.rank=FALSE,
      cshrink = cshrink) 
    p_value <- res$p.value[,1]
    q_value <- res$q.value[,1]
    stat_value <- res$rho[,1]
  } else if(method=="goat" && require("goat")) {
    gmt <- mat2gmt(G)
    res <- goat(pathways=gmt, stats=fc, filter=FALSE, method="goat") 
    res <- res[match(colnames(G),res$pathway),]
    p_value <- res$pval
    q_value <- res$padj
    stat_value <- res$score
  } else {
    stop("unknown method: ", method)
  }

  ## Add leading edge list. The leading edge list is computed as those
  ## genes in the gene set that have fold-change larger than z score
  ## z>minLE (default minLE=1).
  leading_edge <- NULL
  if(addLE) {
    if(is.null(zmat)) zmat <- fc_zmat(fc, G, alpha=alpha)
    zsign <- sign(stat_value)
    leading_idx <- Matrix::which( abs(zmat) >= minLE, arr.ind=TRUE)  
    ii <- which(sign(zmat[leading_idx]) == zsign[leading_idx[,1]])
    leading_idx <- leading_idx[ii,]
    leading_edge <- tapply( leading_idx[,2], leading_idx[,1], list)
    leading_edge <- lapply( leading_edge, function(i) colnames(zmat)[i])
    names(leading_edge) <- rownames(zmat)[as.integer(names(leading_edge))]
    leading_edge <- leading_edge[match(rownames(zmat),names(leading_edge))]
  }
  
  #size <- Matrix::colSums(G!=0)
  if(is.null(q_value)) {
    q_value <- p.adjust(p_value, method="fdr")
  }

  if(format %in% c("as.gsea","as.gsea2")) {
    if(format == "as.gsea") {
      # normalize like GSEA
      stat_value0 <-  stat_value / max(abs(stat_value),na.rm=TRUE) 
      stat_value1 <-  stat_value / sd(stat_value,na.rm=TRUE)
    }
    if(format == "as.gsea2") {
      # run fgseaSimple only for few iterations only for NES
      if(is.null(gmt)) gmt <- mat2gmt(G)  
      suppressWarnings(res <- fgsea::fgseaSimple(gmt, fc, nperm=10))
      res <- res[match(colnames(G),res$pathway),]  
      stat_value0 <-  res$ES
      stat_value1 <-  res$NES
      leading_edge <- res$leadingEdge                
    }
    df <- data.table::data.table(
      pathway = colnames(G),
      pval = p_value,
      padj = q_value,
      log2err = as.numeric(NA),
      ES = stat_value0,
      NES = stat_value1,
      size = size,
      leadingEdge = leading_edge        
    )
  } else {
    ## simple format
    df <- data.frame(
      pathway = colnames(G),
      pval = p_value,
      padj = q_value,
      score = stat_value,
      size = size
    )
    rownames(df) <- names(p_value)
  }
  df
}

#' Fast one sample z-test for matrix object F (e.g. foldchanges) and
#' grouping matrix G (e.g. gene sets).
#'
#' @param fc Numeric vector of (log)foldchanges, with genes as names
#' @param G Numeric matrix of gene sets, with genes as row names
#' @param zmat Logical indicating whether to compute z-score matrix
#' @param alpha Numeric value between 0 and 1 for Bayesian shrinkage
#'
#' @return List with components:
#' \itemize{
#'  \item z_statistic - Vector of z-test statistics
#'  \item p_value - Vector of p-values
#'  \item zmat - Z-score matrix (if zmat = TRUE)
#' }
#' 
#' @export
fc_ztest <- function(F, G, zmat=FALSE, alpha=0.5) {
  if(NCOL(F)==1) F <- cbind(F)
  gg <- intersect(rownames(G),rownames(F))
  sample_size <- Matrix::colSums(G[gg,]!=0)
  sample_size <-sample_size + 1e-20 ## avoid div-by-zero
  sample_mean <- (Matrix::t(G[gg,]!=0) %*% F[gg,,drop=FALSE]) / sample_size
  population_mean <- Matrix::colMeans(F, na.rm=TRUE)
  population_var <- matrixStats::colVars(F, na.rm=TRUE)
  sample_var <- apply(F, 2, function(fc) {
    gfc <- (G[gg,]!=0) * fc[gg]
    sparseMatrixStats::colVars(gfc) * nrow(G) / sample_size
  })
  alpha <- pmin(pmax(alpha,0), 0.999) ## limit
  estim_sd <- Matrix::t(alpha*sample_var) + (1-alpha)*population_var
  estim_sd <- Matrix::t(sqrt(estim_sd)) / sqrt(sample_size)
  delta_mean <- Matrix::t(Matrix::t(sample_mean) -  population_mean)
  z_statistic <- delta_mean / estim_sd
  z_statistic <- as.matrix(z_statistic)
  p_value <- 2 * pnorm(abs(z_statistic), lower.tail = FALSE)  
  if(zmat) {
    zmat <- lapply(1:ncol(F), function(i) {
      gfc <- (G[gg,]!=0) * F[gg,i]
      (Matrix::t(gfc) / estim_sd[,i])
    })
  } else {
    zmat = NULL
  }
  if(NCOL(F)==1) {
    ## collapse
    z_statistic <- z_statistic[,1]
    p_value <- p_value[,1]
    zmat <- zmat[[1]]
  }
  list(
    z = z_statistic,
    p = p_value,
    zmat = zmat
  )
}

#' Compute z-score matrix for gene sets
#'
#' @param fc Numeric vector of (log)foldchanges, with genes as names
#' @param G Numeric matrix of gene sets, with genes as row names
#' @param alpha Numeric value between 0 and 1 for Bayesian shrinkage
#'
#' @return Z-score matrix
#' @keywords internal
fc_zmat <- function(fc, G, alpha=0.5) {
  gg <- intersect(rownames(G),names(fc))
  sample_size <- Matrix::colSums(G[gg,]!=0)
  sample_size <- pmax(sample_size, 1) ## avoid div-by-zero
  population_var <- var(fc, na.rm=TRUE)
  gfc <- (G[gg,]!=0) * fc[gg]
  sample_var <- sparseMatrixStats::colVars(gfc) * nrow(G) / sample_size
  alpha <- pmin(pmax(alpha,0), 0.999) ## limit
  estim_sd <- sqrt( alpha*sample_var + (1-alpha)*population_var )
  zmat <- (Matrix::t(gfc) / estim_sd)
  zmat
}

#' Fast one sample t-test for matrix object F (e.g. foldchanges) and
#' grouping matrix G (e.g. gene sets).
#'
#' @param F Numeric vector or matrix of (log)foldchanges
#' @param G Numeric matrix of gene sets, with genes as row names
#'
#' @return List with components:
#' \itemize{
#'  \item mean - Mean values
#'  \item t - t-test statistics
#'  \item p - p-values from t-test
#' }
#' @keywords internal
matrix_onesample_ttest <- function(F, G) {
  sumG <- Matrix::colSums(G != 0)
  sum_sq <- Matrix::crossprod(G != 0, F^2)
  meanx <- Matrix::crossprod(G != 0, F) / (1e-8 + sumG)
  sdx <- sqrt((sum_sq - meanx^2 * sumG) / (sumG - 1))
  t_stats <- meanx / (1e-8 + sdx) * sqrt(sumG)
  p_stats <- apply(abs(t_stats), 2, function(tv) {
    2 * pt(tv, df = pmax(sumG - 1, 1), lower.tail = FALSE)
  })
  list(mean = as.matrix(meanx), t = as.matrix(t_stats), p = p_stats)
}

#' Fast one sample z-test for matrix object F (e.g. foldchanges) and
#' grouping matrix G (e.g. gene sets).
#'
#' @param F Numeric vector or matrix of (log)foldchanges
#' @param G Numeric matrix of gene sets, with genes as row names
#'
#' @return List with components:
#' \itemize{
#'  \item mean - Mean values
#'  \item z - z-test statistics
#'  \item p - p-values from z-test
#' }
#' @keywords internal
matrix_onesample_ztest <- function(F, G) {
  sumG <- Matrix::colSums(G != 0)
  sum_sq <- Matrix::crossprod(G != 0, F^2)
  meanx <- Matrix::crossprod(G != 0, F) / (1e-8 + sumG)
  sdx <- sqrt((sum_sq - meanx^2 * sumG) / sumG)
  z_stats <- meanx / (1e-8 + sdx) * sqrt(sumG)
  p_stats <- apply(abs(z_stats), 2, function(zv) {
    2 * pnorm(zv, lower.tail = FALSE)
  })
  list(mean = as.matrix(meanx), z = as.matrix(z_stats), p = p_stats)
}
