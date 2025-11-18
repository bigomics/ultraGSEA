# Unit tests for ultraGSEA main function

test_that("ultraGSEA basic functionality works", {
  # Create test data
  set.seed(123)
  fc <- rnorm(100, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:100)

  # Create gene set matrix (binary, sparse)
  G <- matrix(rbinom(100 * 10, 1, 0.3), nrow = 100, ncol = 10)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("geneset_", 1:10)

  # Test basic call
  result <- ultraGSEA(fc, G)

  # Check output structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 5)
  expect_named(result, c("pathway", "pval", "padj", "score", "size"))
})

test_that("ultraGSEA returns correct column types", {
  set.seed(456)
  fc <- rnorm(50)
  names(fc) <- paste0("g", 1:50)
  G <- matrix(rbinom(50 * 5, 1, 0.4), nrow = 50, ncol = 5)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:5)

  result <- ultraGSEA(fc, G)

  expect_is(result$pathway, "character")
  expect_is(result$pval, "numeric")
  expect_is(result$padj, "numeric")
  expect_is(result$score, "numeric")
  expect_is(result$size, "numeric")
})

test_that("ultraGSEA with format='as.gsea' returns data.table", {
  set.seed(789)
  fc <- rnorm(30)
  names(fc) <- paste0("g", 1:30)
  G <- matrix(rbinom(30 * 3, 1, 0.5), nrow = 30, ncol = 3)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:3)

  result <- ultraGSEA(fc, G, format = "as.gsea")

  expect_s3_class(result, "data.table")
  expect_equal(ncol(result), 8)
  expect_named(result, c("pathway", "pval", "padj", "log2err", "ES",
                         "NES", "size", "leadingEdge"))
})

test_that("ultraGSEA with different methods", {
  set.seed(101)
  fc <- rnorm(40)
  names(fc) <- paste0("g", 1:40)
  G <- matrix(rbinom(40 * 4, 1, 0.35), nrow = 40, ncol = 4)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:4)

  # Test ztest method
  result_ztest <- ultraGSEA(fc, G, method = "ztest")
  expect_s3_class(result_ztest, "data.frame")
  expect_equal(nrow(result_ztest), 4)

  # Test ttest method
  result_ttest <- ultraGSEA(fc, G, method = "ttest")
  expect_s3_class(result_ttest, "data.frame")
  expect_equal(nrow(result_ttest), 4)

  # Test cor method
  result_cor <- ultraGSEA(fc, G, method = "cor")
  expect_s3_class(result_cor, "data.frame")
  expect_equal(nrow(result_cor), 4)

  # Test rankcor method
  result_rankcor <- ultraGSEA(fc, G, method = "rankcor")
  expect_s3_class(result_rankcor, "data.frame")
  expect_equal(nrow(result_rankcor), 4)
})

test_that("ultraGSEA with alpha parameter", {
  set.seed(202)
  fc <- rnorm(35)
  names(fc) <- paste0("g", 1:35)
  G <- matrix(rbinom(35 * 3, 1, 0.3), nrow = 35, ncol = 3)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:3)

  # Different alpha values
  result_low <- ultraGSEA(fc, G, alpha = 0.1)
  result_high <- ultraGSEA(fc, G, alpha = 0.9)

  # Results should be different
  expect_false(isTRUE(all.equal(result_low$score, result_high$score)))
})

test_that("ultraGSEA with minLE parameter", {
  set.seed(303)
  fc <- rnorm(25)
  names(fc) <- paste0("g", 1:25)
  G <- matrix(rbinom(25 * 2, 1, 0.5), nrow = 25, ncol = 2)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:2)

  result <- ultraGSEA(fc, G, minLE = 2, format = "as.gsea")

  expect_is(result$leadingEdge, "list")
})

test_that("ultraGSEA handles gene intersection", {
  set.seed(404)
  fc <- rnorm(100)
  names(fc) <- paste0("g", 1:100)

  # Gene set with subset of genes
  G <- matrix(rbinom(50 * 5, 1, 0.4), nrow = 50, ncol = 5)
  rownames(G) <- paste0("g", 51:100)  # Partial overlap
  colnames(G) <- paste0("gs", 1:5)

  result <- ultraGSEA(fc, G)

  expect_equal(nrow(result), 5)
  # Check that size reflects only overlapping genes
  expect_true(all(result$size <= 50))
})

test_that("ultraGSEA produces valid p-values", {
  set.seed(505)
  fc <- rnorm(60)
  names(fc) <- paste0("g", 1:60)
  G <- matrix(rbinom(60 * 6, 1, 0.3), nrow = 60, ncol = 6)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:6)

  result <- ultraGSEA(fc, G)

  # Check p-values are between 0 and 1
  expect_true(all(result$pval >= 0 & result$pval <= 1))
  expect_true(all(result$padj >= 0 & result$padj <= 1))

  # Check FDR adjustment (padj >= pval)
  expect_true(all(result$padj >= result$pval - 1e-10))
})

test_that("ultraGSEA throws error for invalid method", {
  fc <- rnorm(20)
  names(fc) <- paste0("g", 1:20)
  G <- matrix(rbinom(20 * 2, 1, 0.5), nrow = 20, ncol = 2)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:2)

  expect_error(ultraGSEA(fc, G, method = "invalid"),
               "unknown method")
})

test_that("ultraGSEA with sparse matrix input", {
  skip_if_not_installed("Matrix")

  set.seed(606)
  fc <- rnorm(50)
  names(fc) <- paste0("g", 1:50)

  # Create sparse matrix
  G_dense <- matrix(rbinom(50 * 5, 1, 0.3), nrow = 50, ncol = 5)
  rownames(G_dense) <- names(fc)
  colnames(G_dense) <- paste0("gs", 1:5)
  G_sparse <- Matrix::Matrix(G_dense, sparse = TRUE)

  result_sparse <- ultraGSEA(fc, G_sparse)
  result_dense <- ultraGSEA(fc, G_dense)

  # Results should be identical
  expect_equal(as.numeric(result_sparse$score),
               as.numeric(result_dense$score), tolerance = 1e-10)
  expect_equal(as.numeric(result_sparse$pval),
               as.numeric(result_dense$pval), tolerance = 1e-10)
})
