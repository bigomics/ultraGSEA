# Unit tests for gset.cortesttest and cor_sparse_matrix functions

test_that("gset.cortest basic functionality with vector input", {
  set.seed(100)
  FC <- rnorm(50)
  names(FC) <- paste0("g", 1:50)

  gset <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(gset) <- names(FC)
  colnames(gset) <- paste0("gs", 1:5)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result <- gset.cortest(FC, gset)

  expect_type(result, "list")
  expect_named(result, c("rho", "p.value", "q.value"))
  expect_equal(nrow(result$rho), 5)
  expect_equal(ncol(result$rho), 1)
})

test_that("gset.cortest with matrix input", {
  set.seed(200)
  FC <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
  rownames(FC) <- paste0("g", 1:50)
  colnames(FC) <- paste0("sample", 1:3)

  gset <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(gset) <- rownames(FC)
  colnames(gset) <- paste0("gs", 1:5)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result <- gset.cortest(FC, gset)

  expect_type(result, "list")
  expect_equal(nrow(result$rho), 5)
  expect_equal(ncol(result$rho), 3)
})

test_that("gset.cortest with compute.p=TRUE", {
  set.seed(300)
  FC <- rnorm(40)
  names(FC) <- paste0("g", 1:40)

  gset <- matrix(rbinom(40 * 4, 1, 0.3), nrow = 40, ncol = 4)
  rownames(gset) <- names(FC)
  colnames(gset) <- paste0("gs", 1:4)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result <- gset.cortest(FC, gset, compute.p = TRUE)

  expect_true(is.matrix(result$p.value))
  expect_true(is.matrix(result$q.value))
  expect_equal(nrow(result$p.value), 4)
  expect_true(all(result$p.value >= 0 & result$p.value <= 1, na.rm = TRUE))
})

test_that("gset.cortest with compute.p=FALSE", {
  set.seed(400)
  FC <- rnorm(30)
  names(FC) <- paste0("g", 1:30)

  gset <- matrix(rbinom(30 * 3, 1, 0.4), nrow = 30, ncol = 3)
  rownames(gset) <- names(FC)
  colnames(gset) <- paste0("gs", 1:3)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result <- gset.cortest(FC, gset, compute.p = FALSE)

  expect_true(is.matrix(result$rho))
  expect_equal(result$p.value, NA)
  expect_equal(result$q.value, NA)
})

test_that("gset.cortest with use.rank=TRUE", {
  set.seed(500)
  FC <- rnorm(45)
  names(FC) <- paste0("g", 1:45)

  gset <- matrix(rbinom(45 * 5, 1, 0.3), nrow = 45, ncol = 5)
  rownames(gset) <- names(FC)
  colnames(gset) <- paste0("gs", 1:5)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result_cor <- gset.cortest(FC, gset, use.rank = FALSE)
  result_rankcor <- gset.cortest(FC, gset, use.rank = TRUE)

  # Results should be different
  expect_false(isTRUE(all.equal(result_cor$rho, result_rankcor$rho)))
})

test_that("gset.cortest with transposed gene set", {
  set.seed(600)
  FC <- rnorm(50)
  names(FC) <- paste0("g", 1:50)

  gset <- matrix(rbinom(5 * 50, 1, 0.3), nrow = 5, ncol = 50)
  colnames(gset) <- names(FC)
  rownames(gset) <- paste0("gs", 1:5)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result <- gset.cortest(FC, gset)

  # Should auto-transpose to correct orientation
  expect_equal(nrow(result$rho), 5)
})

test_that("gset.cortest returns valid correlations", {
  set.seed(700)
  FC <- rnorm(40)
  names(FC) <- paste0("g", 1:40)

  gset <- matrix(rbinom(40 * 4, 1, 0.3), nrow = 40, ncol = 4)
  rownames(gset) <- names(FC)
  colnames(gset) <- paste0("gs", 1:4)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result <- gset.cortest(FC, gset)

  # Correlations should be between -1 and 1 (or NA)
  expect_true(all(result$rho >= -1 & result$rho <= 1, na.rm = TRUE))
})

test_that("gset.cortest throws error for missing names", {
  FC <- rnorm(30)
  gset <- matrix(rbinom(30 * 3, 1, 0.3), nrow = 30, ncol = 3)
  rownames(gset) <- paste0("g", 1:30)
  colnames(gset) <- paste0("gs", 1:3)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  expect_error(gset.cortest(FC, gset), "must be named")
})

test_that("gset.cortest throws error for zero columns", {
  set.seed(800)
  FC <- rnorm(30)
  names(FC) <- paste0("g", 1:30)

  gset <- matrix(0, nrow = 30, ncol = 0)
  rownames(gset) <- names(FC)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  expect_error(gset.cortest(FC, gset), "zero columns")
})

test_that("gset.cortest with gene intersection", {
  set.seed(900)
  FC <- rnorm(100)
  names(FC) <- paste0("g", 1:100)

  gset <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(gset) <- paste0("g", 51:100)  # Partial overlap
  colnames(gset) <- paste0("gs", 1:5)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result <- gset.cortest(FC, gset)

  expect_equal(nrow(result$rho), 5)
})

test_that("cor_sparse_matrix without missing values", {
  set.seed(1000)
  G <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  G <- Matrix::Matrix(G, sparse = TRUE)

  mat <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)

  result <- cor_sparse_matrix(G, mat)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 3)
})

test_that("cor_sparse_matrix with missing values", {
  set.seed(1100)
  G <- matrix(rbinom(40 * 4, 1, 0.3), nrow = 40, ncol = 4)
  G <- Matrix::Matrix(G, sparse = TRUE)

  mat <- matrix(rnorm(40 * 2), nrow = 40, ncol = 2)
  mat[sample(1:40, 5), 1] <- NA  # Add some missing values

  result <- cor_sparse_matrix(G, mat)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 2)
})

test_that("gset.cortest FDR adjustment", {
  set.seed(1200)
  FC <- rnorm(50)
  names(FC) <- paste0("g", 1:50)

  gset <- matrix(rbinom(50 * 10, 1, 0.3), nrow = 50, ncol = 10)
  rownames(gset) <- names(FC)
  colnames(gset) <- paste0("gs", 1:10)
  gset <- Matrix::Matrix(gset, sparse = TRUE)

  result <- gset.cortest(FC, gset, compute.p = TRUE)

  # FDR adjustment should produce q-values >= p-values
  expect_true(all(result$q.value >= result$p.value - 1e-10, na.rm = TRUE))
})
