##---------------------------------------------------------------
##-------------------- ultragsea --------------------------------
##---------------------------------------------------------------

#' Ultra-fast GSEA using z-score or geneset correlation. Results of
#' these methods highly correlate with GSEA/fGSEA but are much faster.
#'
#' @export
ultragsea <- function(G, fc, alpha=0.5, minLE=1, corshrink=3,
                      minsize = 1L, maxsize = 9999L, center=TRUE,
                      method=c("cor","ztest","ttest","goat","camera")[1],
                      format=c("simple","as.gsea")[1]) {

  gg <- intersect(names(fc), rownames(G))
  fc <- fc[gg]
  G <- G[gg,]
  size <- Matrix::colSums(G!=0)
  if(method == "goat") minsize <- max(minsize,10L)
  sel <- which(size >= minsize & size <= maxsize)
  G <- G[, sel, drop=FALSE]
  size <- size[sel]

  gmt <- NULL
  if(method %in% c("goat","camera") && is.null(gmt)) {
    gmt <- mat2gmt(G)
  }
  
  addLE <- (format=="as.gsea")
  p_value <- NULL
  q_value <- NULL
  stat_value <- NULL
  if(method == "ztest") {
    zres <- gset.ztest(G, fc, alpha=alpha, center=center)
    p_value <- zres$p
    stat_value <- zres$stat
  } else if(method == "ttest") {
    zres <- gset.ztest(G, fc, alpha=alpha, pdist="norm", center=center)
    p_value <- zres$p
    stat_value <- zres$stat
  } else if(method == "cor") {
    res <- gset.cor(G, fc, compute.p=TRUE, use.rank=FALSE,
      corshrink = corshrink) 
    p_value <- res$p.value[,1]
    q_value <- res$q.value[,1]
    stat_value <- res$rho[,1]
  } else if(method=="goat" && require("goat")) {
    res <- goat(pathways=gmt, stats=fc, filter=FALSE, method="goat") 
    res <- res[match(colnames(G),res$pathway),]
    p_value <- res$pval
    q_value <- res$padj
    stat_value <- res$score
  } else if(method=="camera" && require("limma")) {
    res <- limma::cameraPR(fc, gmt)
    res <- res[match(colnames(G),rownames(res)),]
    p_value <- res$PValue
    q_value <- res$FDR
    stat_value <- -log10(res$PValue)*(1 - 2*(res$Direction=="Down"))
  } else {
    stop("unsupported method: ", method)
  }

  ## Add leading edge list. The leading edge list is computed as those
  ## genes in the gene set that have fold-change larger than z score
  ## z>minLE (default minLE=1).
  leading_edge <- NULL
  if(addLE) {
    zsign <- sign(stat_value) 
    zmat <- sweep((G!=0), 2, zsign, FUN='*')
    zmat <- (zmat * fc > minLE) ## signed            
    leading_idx <- Matrix::which(zmat == 1, arr.ind=TRUE)  
    leading_edge <- tapply( leading_idx[,1], leading_idx[,2], list)
    leading_edge <- lapply( leading_edge, function(i) rownames(zmat)[i])
    names(leading_edge) <- colnames(zmat)[as.integer(names(leading_edge))]
    leading_edge <- leading_edge[match(colnames(zmat),names(leading_edge))]
    names(leading_edge) <- colnames(zmat)
  }
  
  #size <- Matrix::colSums(G!=0)
  if(is.null(q_value)) {
    q_value <- p.adjust(p_value, method="fdr")
  }

  if(format == "as.gsea") {
    # normalize like GSEA
    stat_value0 <-  stat_value / max(abs(stat_value),na.rm=TRUE) 
    stat_value1 <-  stat_value / sd(stat_value,na.rm=TRUE)
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

#' Fast replacement of fgsea with p/q values computed with ultragsea.cor
#' 
#' @export
fgsea <- function(pathways, stats, minSize = 1, maxSize = length(stats)-1,
                  nperm=20) {
  G <- gmt2mat(pathways)
  gg <- intersect(names(stats),rownames(G))
  G <- G[gg,match(names(pathways),colnames(G))]
  stats <- stats[gg]
  gsize <- Matrix::colSums(G)
  G <- G[, gsize >= minSize & gsize <= maxSize]
  res <- ultragsea(G, stats, method="ztest", alpha=0, format="as.gsea2")
  # run fgseaSimple only for few iterations only for NES
  suppressWarnings(fres <- fgsea::fgseaSimple(pathways, stats, nperm=nperm))
  res <- res[match(res$pathway,fres$pathway),]  
  fres$pval <-  fres$pval
  fres$padj <-  fres$padj  
  fres
}


#' Fast one sample z-test for matrix object F (e.g. foldchanges) and
#' grouping matrix G (e.g. gene sets).
#'
#' @param fc Numeric vector of (log)foldchanges, with genes as names
#' @param G Numeric matrix of gene sets, with genes as row names
#' @param alpha Numeric value between 0 and 1 for Bayesian shrinkage
#' @param center Logical value to center population mean or assume zero
#'
#' @return List with components:
#' \itemize{
#'  \item z_statistic - Vector of z-test statistics
#'  \item p_value - Vector of p-values
#' }
#' 
#' @export
gset.ztest <- function(G, F, alpha=0.5, center=TRUE, pdist="norm") {
  if(NCOL(F)==1) F <- cbind(F)
  gg <- intersect(rownames(G),rownames(F))
  G <- G[gg,,drop=FALSE]
  F <- F[gg,,drop=FALSE]  
  gset_size <- Matrix::colSums(G!=0)
  gset_size <-gset_size + 1e-20 ## avoid div-by-zero
  gset_mean <- (Matrix::t(G!=0) %*% F) / gset_size
  if(center) {
    population_mean <- Matrix::colMeans(F, na.rm=TRUE)
    population_var <- matrixStats::colVars(F, na.rm=TRUE)
  } else {
    population_mean <- rep(0, ncol(F))    
    population_var <- matrixStats::colVars(F, center=rep(0,ncol(F)), na.rm=TRUE)
  }
  if(alpha == 0) {
    ##estim_var <- population_var
    estim_var <- matrix(population_var, ncol=ncol(F), nrow=ncol(G), byrow=TRUE)
  } else {
    alpha <- pmin(pmax(alpha,0), 0.999) ## limit
    gset_sumsq <- Matrix::t(G!=0) %*% F**2 / gset_size
    gset_var <- (gset_sumsq - gset_mean**2) * ( gset_size / (gset_size - 1))
    estim_var <- sweep( alpha*gset_var, 2, (1-alpha)*population_var, '+')
  }
  estim_sd <- sqrt(estim_var) / sqrt(gset_size)
  delta_mean <- sweep(gset_mean, 2, population_mean, FUN='-')
  statistic <- as.matrix( delta_mean / estim_sd )
  if(pdist=="t") {
    df <- gset_size - 1
    p_value <- 2 * pt(abs(statistic), df=df, lower.tail = FALSE)
  } else {
    p_value <- 2 * pnorm(abs(statistic), lower.tail = FALSE)
  } 
  if(NCOL(F)==1) {
    ## collapse
    statistic <- statistic[,1]
    p_value <- p_value[,1]
  }
  list(
    stat = statistic,
    p = p_value
  )
}


