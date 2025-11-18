# Unit tests for fc_ztest and fc_zmat functions

test_that("fc_ztest basic functionality", {
  set.seed(100)
  fc <- rnorm(50)
  names(fc) <- paste0("g", 1:50)

  G <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:5)

  result <- fc_ztest(fc, G)

  expect_type(result, "list")
  expect_named(result, c("z_statistic", "p_value", "zmat"))
  expect_equal(length(result$z_statistic), 5)
  expect_equal(length(result$p_value), 5)
  expect_null(result$zmat)
})

test_that("fc_ztest with zmat=TRUE returns z-score matrix", {
  set.seed(200)
  fc <- rnorm(40)
  names(fc) <- paste0("g", 1:40)

  G <- matrix(rbinom(40 * 4, 1, 0.35), nrow = 40, ncol = 4)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:4)

  result <- fc_ztest(fc, G, zmat = TRUE)

  expect_s4_class(result$zmat, "Matrix")
  expect_equal(nrow(result$zmat), 4)
  expect_equal(ncol(result$zmat), 40)
})

test_that("fc_ztest p-values are valid", {
  set.seed(300)
  fc <- rnorm(60)
  names(fc) <- paste0("g", 1:60)

  G <- matrix(rbinom(60 * 6, 1, 0.3), nrow = 60, ncol = 6)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:6)

  result <- fc_ztest(fc, G)

  # P-values should be between 0 and 1
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
  expect_true(all(!is.na(result$p_value)))
  expect_true(all(!is.nan(result$p_value)))
})

test_that("fc_ztest with alpha parameter", {
  set.seed(400)
  fc <- rnorm(30)
  names(fc) <- paste0("g", 1:30)

  G <- matrix(rbinom(30 * 3, 1, 0.4), nrow = 30, ncol = 3)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:3)

  result_low <- fc_ztest(fc, G, alpha = 0.1)
  result_high <- fc_ztest(fc, G, alpha = 0.9)

  # Different alpha should give different results
  expect_false(isTRUE(all.equal(result_low$z_statistic,
                                 result_high$z_statistic)))
})

test_that("fc_ztest handles gene intersection", {
  set.seed(500)
  fc <- rnorm(100)
  names(fc) <- paste0("g", 1:100)

  G <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(G) <- paste0("g", 51:100)  # Only genes 51-100
  colnames(G) <- paste0("gs", 1:5)

  result <- fc_ztest(fc, G)

  # Should still return results for all 5 gene sets
  expect_equal(length(result$z_statistic), 5)
})

test_that("fc_zmat returns correct matrix dimensions", {
  set.seed(600)
  fc <- rnorm(45)
  names(fc) <- paste0("g", 1:45)

  G <- matrix(rbinom(45 * 7, 1, 0.3), nrow = 45, ncol = 7)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:7)

  zmat <- fc_zmat(fc, G)

  expect_s4_class(zmat, "Matrix")
  expect_equal(nrow(zmat), 7)
  expect_equal(ncol(zmat), 45)
})

test_that("fc_zmat with different alpha values", {
  set.seed(700)
  fc <- rnorm(35)
  names(fc) <- paste0("g", 1:35)

  G <- matrix(rbinom(35 * 3, 1, 0.4), nrow = 35, ncol = 3)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:3)

  zmat_low <- fc_zmat(fc, G, alpha = 0.1)
  zmat_high <- fc_zmat(fc, G, alpha = 0.9)

  # Different alpha should give different z-scores
  if (length(zmat_low@x) > 0 && length(zmat_high@x) > 0) {
    expect_false(isTRUE(all.equal(zmat_low@x, zmat_high@x)))
  }
})

test_that("fc_zmat values are reasonable", {
  set.seed(800)
  fc <- rnorm(50, mean = 5, sd = 2)
  names(fc) <- paste0("g", 1:50)

  G <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:5)

  zmat <- fc_zmat(fc, G)

  # Z-scores should be numeric and finite (allowing 0s in sparse matrix)
  non_zero_vals <- zmat@x
  expect_true(all(is.finite(non_zero_vals)))
})

test_that("fc_ztest z-statistics and p-values are consistent", {
  set.seed(900)
  fc <- rnorm(40)
  names(fc) <- paste0("g", 1:40)

  G <- matrix(rbinom(40 * 4, 1, 0.3), nrow = 40, ncol = 4)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:4)

  result <- fc_ztest(fc, G)

  # For a z-test, p-value = 2 * P(Z > |z|)
  # Manually compute and verify
  expected_p <- 2 * pnorm(abs(result$z_statistic), lower.tail = FALSE)

  expect_equal(result$p_value, expected_p, tolerance = 1e-10)
})

test_that("fc_ztest with single gene set", {
  set.seed(1000)
  fc <- rnorm(30)
  names(fc) <- paste0("g", 1:30)

  G <- matrix(rbinom(30 * 1, 1, 0.3), nrow = 30, ncol = 1)
  rownames(G) <- names(fc)
  colnames(G) <- "single_geneset"

  result <- fc_ztest(fc, G)

  expect_equal(length(result$z_statistic), 1)
  expect_equal(length(result$p_value), 1)
})

test_that("fc_ztest with large dataset", {
  skip_on_ci()  # Skip on CI to save time

  set.seed(1100)
  fc <- rnorm(5000)
  names(fc) <- paste0("g", 1:5000)

  G <- matrix(rbinom(5000 * 100, 1, 0.2), nrow = 5000, ncol = 100)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:100)

  result <- fc_ztest(fc, G)

  expect_equal(length(result$z_statistic), 100)
  expect_true(all(!is.na(result$p_value)))
})
