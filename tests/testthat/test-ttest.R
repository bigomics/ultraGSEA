# Unit tests for matrix_onesample_ttest function

test_that("matrix_onesample_ttest basic functionality", {
  set.seed(100)
  F <- matrix(rnorm(50), nrow = 50, ncol = 1)
  rownames(F) <- paste0("g", 1:50)

  G <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(G) <- rownames(F)
  colnames(G) <- paste0("gs", 1:5)

  result <- matrix_onesample_ttest(F, G)

  expect_type(result, "list")
  expect_named(result, c("mean", "t", "p"))
  expect_equal(nrow(result$mean), 5)
  expect_equal(nrow(result$t), 5)
  expect_equal(nrow(result$p), 5)
})

test_that("matrix_onesample_ttest returns correct dimensions", {
  set.seed(200)
  F <- matrix(rnorm(60 * 3), nrow = 60, ncol = 3)
  rownames(F) <- paste0("g", 1:60)

  G <- matrix(rbinom(60 * 4, 1, 0.35), nrow = 60, ncol = 4)
  rownames(G) <- rownames(F)
  colnames(G) <- paste0("gs", 1:4)

  result <- matrix_onesample_ttest(F, G)

  expect_equal(nrow(result$mean), 4)
  expect_equal(ncol(result$mean), 3)
  expect_equal(nrow(result$t), 4)
  expect_equal(ncol(result$t), 3)
  expect_equal(nrow(result$p), 4)
  expect_equal(ncol(result$p), 3)
})

test_that("matrix_onesample_ttest with single sample", {
  set.seed(300)
  F <- matrix(rnorm(40), nrow = 40, ncol = 1)
  rownames(F) <- paste0("g", 1:40)

  G <- matrix(rbinom(40 * 3, 1, 0.3), nrow = 40, ncol = 3)
  rownames(G) <- rownames(F)
  colnames(G) <- paste0("gs", 1:3)

  result <- matrix_onesample_ttest(F, G)

  expect_equal(ncol(result$mean), 1)
  expect_equal(ncol(result$t), 1)
  expect_equal(ncol(result$p), 1)
})

test_that("matrix_onesample_ttest with multiple samples", {
  set.seed(400)
  F <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(F) <- paste0("g", 1:50)

  G <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(G) <- rownames(F)
  colnames(G) <- paste0("gs", 1:5)

  result <- matrix_onesample_ttest(F, G)

  expect_equal(ncol(result$mean), 10)
  expect_equal(ncol(result$t), 10)
  expect_equal(ncol(result$p), 10)
})

test_that("matrix_onesample_ttest p-values are valid", {
  set.seed(500)
  F <- matrix(rnorm(45 * 2), nrow = 45, ncol = 2)
  rownames(F) <- paste0("g", 1:45)

  G <- matrix(rbinom(45 * 4, 1, 0.3), nrow = 45, ncol = 4)
  rownames(G) <- rownames(F)
  colnames(G) <- paste0("gs", 1:4)

  result <- matrix_onesample_ttest(F, G)

  # P-values should be between 0 and 1
  expect_true(all(result$p >= 0 & result$p <= 1, na.rm = TRUE))
  expect_false(any(is.infinite(result$p)))
})

test_that("matrix_onesample_ttest with sparse matrix", {
  skip_if_not_installed("Matrix")

  set.seed(600)
  F <- matrix(rnorm(35), nrow = 35, ncol = 1)
  rownames(F) <- paste0("g", 1:35)

  G_dense <- matrix(rbinom(35 * 4, 1, 0.3), nrow = 35, ncol = 4)
  rownames(G_dense) <- rownames(F)
  colnames(G_dense) <- paste0("gs", 1:4)
  G_sparse <- Matrix::Matrix(G_dense, sparse = TRUE)

  result_sparse <- matrix_onesample_ttest(F, G_sparse)
  result_dense <- matrix_onesample_ttest(F, G_dense)

  # Results should be very similar
  expect_equal(result_sparse$t, result_dense$t, tolerance = 1e-10)
})

test_that("matrix_onesample_ttest with large foldchanges", {
  set.seed(700)
  F <- matrix(c(rep(10, 25), rep(-10, 25)), ncol = 1)
  rownames(F) <- paste0("g", 1:50)

  G <- matrix(c(rep(1, 25), rep(0, 25)), nrow = 50, ncol = 1)
  rownames(G) <- rownames(F)
  colnames(G) <- "gene_set"

  result <- matrix_onesample_ttest(F, G)

  # Should have significant t-statistic
  expect_true(abs(result$t[1]) > 1)
  expect_true(result$p[1] < 0.05)
})

test_that("matrix_onesample_ttest with null effect", {
  set.seed(800)
  F <- matrix(rnorm(50), nrow = 50, ncol = 1)
  rownames(F) <- paste0("g", 1:50)

  # Gene set with genes randomly selected (no signal)
  G <- matrix(rbinom(50 * 3, 1, 0.5), nrow = 50, ncol = 3)
  rownames(G) <- rownames(F)
  colnames(G) <- paste0("gs", 1:3)

  result <- matrix_onesample_ttest(F, G)

  # P-values should be relatively high (on average)
  expect_true(mean(result$p) > 0.1)
})

test_that("matrix_onesample_ttest handles single gene gene sets", {
  set.seed(900)
  F <- matrix(rnorm(30), nrow = 30, ncol = 1)
  rownames(F) <- paste0("g", 1:30)

  # Gene sets with single genes
  G <- diag(30)[1:5, ]  # Select first 5 genes
  rownames(G) <- rownames(F)[1:5]
  colnames(G) <- paste0("gs", 1:5)

  result <- matrix_onesample_ttest(F, G)

  expect_equal(nrow(result$t), 5)
})

test_that("matrix_onesample_ttest mean computation", {
  set.seed(1000)
  # Create known test data
  F <- matrix(c(1, 2, 3, 0, 0), nrow = 5, ncol = 1)
  rownames(F) <- paste0("g", 1:5)

  G <- matrix(c(1, 1, 1, 0, 0), nrow = 5, ncol = 1)
  rownames(G) <- rownames(F)
  colnames(G) <- "gs1"

  result <- matrix_onesample_ttest(F, G)

  # Mean should be (1+2+3)/3 = 2
  expect_equal(result$mean[1], 2, tolerance = 1e-7)
})

test_that("matrix_onesample_ttest with negative values", {
  set.seed(1100)
  F <- matrix(c(rep(-5, 20), rep(5, 20)), nrow = 40, ncol = 1)
  rownames(F) <- paste0("g", 1:40)

  G <- matrix(c(rep(1, 20), rep(0, 20)), nrow = 40, ncol = 1)
  rownames(G) <- rownames(F)
  colnames(G) <- "negative_genes"

  result <- matrix_onesample_ttest(F, G)

  # t-statistic should be negative (mean < 0)
  expect_true(result$t[1, 1] < 0)
})
