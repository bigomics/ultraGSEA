
#' Calculate gene set correlation
#'
#' Compute correlation between a foldchange vector/matrix and gene sets
#'
#' @param FC Numeric vector or matrix of (log)foldchanges, with genes as row names
#' @param gset Numeric matrix of gene sets, with genes as row/column names
#' @param compute.p Logical indicating whether to compute p-values
#' @param use.rank Logical indicating whether to rank transform FC before correlation
#'
#' @return Named list with components:
#' \itemize{
#'  \item rho - Matrix of correlation coefficients between FC and gset
#'  \item p.value - Matrix of p-values for correlation (if compute.p = TRUE)
#'  \item q.value - Matrix of FDR adjusted p-values (if compute.p = TRUE)
#' }
#'
#' @details This function calculates sparse rank correlation between FC and each
#' column of gset using \code{qlcMatrix::corSparse()}. It handles missing values in
#' FC by computing column-wise correlations.
#'
#' P-values are computed from statistical distribution
#'
#' @examples
#' \dontrun{
#' library(playbase)
#' ranks <- sample(1:10000, 1000, replace = TRUE)
#' names(ranks) <- replicate(1000, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
#' genesets <- matrix(rnorm(1000 * 20), ncol = 20)
#' rownames(genesets) <- names(ranks)
#'
#' gset.rankcor(ranks, genesets, compute.p = TRUE)
#' }
#' @export
gset.cor <- function(gset, FC, compute.p = FALSE, use.rank = FALSE,
                     corshrink = 0) {
  if (ncol(gset) == 0 || NCOL(FC) == 0) {
    if (ncol(gset) == 0) stop("gset has zero columns")
    if (NCOL(FC) == 0) stop("FC has zero columns")
  }

  #  if (!any(class(gset) %in% c("Matrix", "dgCMatrix", "lgCMatrix", "matrix", "array"))) {
  #    stop("gset must be a matrix")
  #  }
  if (!inherits(gset, "Matrix")) stop("gset must be a matrix")

  is.vec <- (NCOL(FC) == 1 && !inherits(FC, c("matrix", "Matrix")))
  if (is.vec && is.null(names(FC))) stop("rank vector must be named")
  if (!is.vec && is.null(rownames(FC))) stop("rank matrix must have rownames")
  if (is.vec) FC <- matrix(FC, ncol = 1, dimnames = list(names(FC), "FC"))
  n1 <- sum(rownames(FC) %in% colnames(gset), na.rm = TRUE)
  n2 <- sum(rownames(FC) %in% rownames(gset), na.rm = TRUE)
  if (n1 > n2) gset <- Matrix::t(gset)

  # align matrices
  gg <- intersect(rownames(gset), rownames(FC))
  FC1 <- FC[gg, , drop = FALSE]
  gset <- gset[gg, , drop = FALSE]

  if (use.rank) {
    if (inherits(FC1, "dgCMatrix")) {
      ## for sparse dgCMatrix
      ## FC1 <- apply(FC1, 2, base::rank, na.last = "keep", ties.method="random")
      FC1 <- sparseMatrixStats::colRanks(FC1, na.last = "keep", ties.method = "random", preserveShape = TRUE)
    } else {
      FC1 <- matrixStats::colRanks(FC1, na.last = "keep", ties.method = "random", preserveShape = TRUE)
    }
  }

  ## two cases: (1) in case no missing values, just use corSparse on
  ## whole matrix. (2) in case the FC matrix has missing values, we
  ## must proceed 1-column at time and do reduced corSparse on
  ## intersection of genes.
  rho1 <- cor_sparse_matrix(gset, FC1)

  rownames(rho1) <- colnames(gset)
  colnames(rho1) <- colnames(FC1)
  rho1[is.nan(rho1)] <- NA ## ??

  if(corshrink > 0) {
    ## Shrinkage to penalize small gene sets.
    gsize <- Matrix::colSums(gset!=0)
    rho1 <- rho1 / (1 + (corshrink / gsize)) 
  }
  
  ## compute p-value
  .cor.pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))
  if (compute.p) {
    pv <- apply(rho1, 2, function(x) .cor.pvalue(x, n = nrow(FC1)))
    pv[is.nan(pv)] <- NA ## ??
    qv <- apply(pv, 2, stats::p.adjust, method = "fdr")
    df <- list(rho = rho1, p.value = pv, q.value = qv)
  } else {
    df <- list(rho = rho1, p.value = NA, q.value = NA)
  }
  df
}

#' Calculate sparse correlation matrix handling missing values
#'
#' @param G Sparse matrix containing gene sets
#' @param mat Matrix of values
#' @return Correlation matrix between G and mat
#' @details If mat has no missing values, calculates correlation directly using corSparse.
#' Otherwise computes column-wise correlations only using non-missing values.
#'
cor_sparse_matrix <- function(G, mat) {
  if (sum(is.na(mat)) == 0) {
    cor_matrix <- qlcMatrix::corSparse(G, mat)
  } else {
    message("matrix has missing values: computing column-wise reduced cor")
    corSparse.vec <- function(X, y) {
      jj <- which(!is.na(y))
      qlcMatrix::corSparse(X[jj, , drop = FALSE], cbind(y[jj]))
    }
    cor_matrix <- lapply(1:ncol(mat), function(i) corSparse.vec(G, mat[, i]))
    cor_matrix <- do.call(cbind, cor_matrix)
  }
  return(cor_matrix)
}
