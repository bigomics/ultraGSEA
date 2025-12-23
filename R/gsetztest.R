
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
gset.ztest <- function(G, F, alpha=0, center=TRUE, pdist="norm") {
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


