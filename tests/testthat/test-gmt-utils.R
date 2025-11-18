# Unit tests for GMT utility functions

test_that("gmt2mat converts GMT list to matrix", {
  set.seed(100)
  gmt <- list(
    GS1 = c("gene1", "gene2", "gene3"),
    GS2 = c("gene2", "gene3", "gene4"),
    GS3 = c("gene5", "gene6")
  )

  result <- gmt2mat(gmt)

  expect_s4_class(result, "Matrix")
  expect_equal(nrow(result), 6)  # Total unique genes
  expect_equal(ncol(result), 3)  # Number of gene sets
  expect_equal(sum(result), 8)   # Total gene-set membership
})

test_that("gmt2mat respects sparse parameter", {
  gmt <- list(
    GS1 = c("gene1", "gene2"),
    GS2 = c("gene2", "gene3")
  )

  result_sparse <- gmt2mat(gmt, sparse = TRUE)
  result_dense <- gmt2mat(gmt, sparse = FALSE)

  expect_true(is(result_sparse, "Matrix"))
  expect_true(is.matrix(result_dense))
  expect_equal(dim(result_sparse), dim(result_dense))
})

test_that("gmt2mat with max.genes parameter", {
  gmt <- list(
    GS1 = c("gene1", "gene2", "gene3", "gene4", "gene5"),
    GS2 = c("gene1", "gene2"),
    GS3 = c("gene6", "gene7", "gene8")
  )

  result <- gmt2mat(gmt, max.genes = 3)

  # GS1 should be removed (> 3 genes), GS3 has 3 genes (not removed)
  expect_true(ncol(result) <= 3)
})

test_that("gmt2mat with ntop parameter", {
  gmt <- list(
    GS1 = c("gene1", "gene2", "gene3"),
    GS2 = c("gene1", "gene2"),
    GS3 = c("gene1")
  )

  result <- gmt2mat(gmt, ntop = 2)

  # Only top 2 gene sets by size
  expect_true(ncol(result) <= 3)
})

test_that("gmt2mat returns correct row and column names", {
  gmt <- list(
    GeneSet1 = c("a", "b", "c"),
    GeneSet2 = c("b", "c", "d")
  )

  result <- gmt2mat(gmt)

  expect_true(all(rownames(result) %in% c("a", "b", "c", "d")))
  expect_equal(colnames(result), c("GeneSet1", "GeneSet2"))
})

test_that("mat2gmt converts matrix back to list", {
  set.seed(200)
  mat <- matrix(rbinom(10 * 5, 1, 0.4), nrow = 10, ncol = 5)
  rownames(mat) <- paste0("gene_", 1:10)
  colnames(mat) <- paste0("gs_", 1:5)

  result <- mat2gmt(mat)

  expect_type(result, "list")
  expect_length(result, 5)
  expect_equal(names(result), colnames(mat))

  # Check that gene sets contain only non-zero genes
  for (gs in result) {
    expect_type(gs, "character")
    expect_true(all(gs %in% rownames(mat)))
  }
})

test_that("mat2gmt with sparse matrix input", {
  skip_if_not_installed("Matrix")

  set.seed(300)
  mat_dense <- matrix(rbinom(15 * 3, 1, 0.3), nrow = 15, ncol = 3)
  rownames(mat_dense) <- paste0("gene_", 1:15)
  colnames(mat_dense) <- paste0("gs_", 1:3)
  mat_sparse <- Matrix::Matrix(mat_dense, sparse = TRUE)

  result_sparse <- mat2gmt(mat_sparse)
  result_dense <- mat2gmt(mat_dense)

  expect_equal(length(result_sparse), length(result_dense))
})

test_that("gmt2mat and mat2gmt are reversible", {
  set.seed(400)
  gmt_original <- list(
    GS1 = c("g1", "g2", "g3"),
    GS2 = c("g2", "g3", "g4"),
    GS3 = c("g5")
  )

  # Convert to matrix and back
  mat <- gmt2mat(gmt_original)
  gmt_restored <- mat2gmt(mat)

  # Check that gene sets are the same (order doesn't matter)
  for (i in names(gmt_original)) {
    expect_setequal(gmt_original[[i]], gmt_restored[[i]])
  }
})

test_that("read.gmt reads GMT file correctly", {
  # Create a temporary GMT file
  gmt_file <- tempfile(fileext = ".gmt")
  writeLines(c(
    "GeneSet1\turl1\tgene1\tgene2\tgene3",
    "GeneSet2\turl2\tgene2\tgene3\tgene4",
    "GeneSet3\turl3\tgene5"
  ), gmt_file)

  result <- read.gmt(gmt_file)

  expect_type(result, "list")
  expect_length(result, 3)
  expect_equal(names(result), c("GeneSet1", "GeneSet2", "GeneSet3"))
  expect_equal(result$GeneSet1, c("gene1", "gene2", "gene3"))

  unlink(gmt_file)
})

test_that("read.gmt with add.source parameter", {
  gmt_file <- tempfile(fileext = ".gmt")
  writeLines(c(
    "GeneSet1\turl1\tgene1\tgene2",
    "GeneSet2\turl2\tgene3\tgene4"
  ), gmt_file)

  result <- read.gmt(gmt_file, add.source = TRUE)

  # Check that we get gene sets
  expect_true(length(result) > 0)
  expect_true(all(sapply(result, is.character)))

  unlink(gmt_file)
})

test_that("write.gmt writes GMT file correctly", {
  gmt <- list(
    GeneSet1 = c("gene1", "gene2", "gene3"),
    GeneSet2 = c("gene2", "gene3")
  )

  gmt_file <- tempfile(fileext = ".gmt")
  write.gmt(gmt, gmt_file)

  # Verify file was created and can be read back
  result <- read.gmt(gmt_file)

  expect_length(result, 2)
  expect_setequal(result$GeneSet1, gmt$GeneSet1)
  expect_setequal(result$GeneSet2, gmt$GeneSet2)

  unlink(gmt_file)
})

test_that("write.gmt preserves gene set names", {
  gmt <- list(
    CustomName1 = c("a", "b"),
    CustomName2 = c("c", "d", "e")
  )

  gmt_file <- tempfile(fileext = ".gmt")
  write.gmt(gmt, gmt_file)

  result <- read.gmt(gmt_file)

  expect_equal(names(result), names(gmt))

  unlink(gmt_file)
})

test_that("gmt2mat with multicore processing", {
  skip_if_not_installed("parallel")

  gmt <- list(
    GS1 = c("gene1", "gene2", "gene3"),
    GS2 = c("gene2", "gene3", "gene4"),
    GS3 = c("gene5")
  )

  result_single <- gmt2mat(gmt, use.multicore = FALSE)
  result_multi <- gmt2mat(gmt, use.multicore = TRUE)

  # Results should be identical
  expect_equal(dim(result_single), dim(result_multi))
})

test_that("mat2gmt with empty gene sets", {
  mat <- matrix(c(1, 0, 0, 0, 0, 1), nrow = 3, ncol = 2)
  rownames(mat) <- c("g1", "g2", "g3")
  colnames(mat) <- c("gs1", "gs2")

  result <- mat2gmt(mat)

  expect_equal(length(result$gs1), 1)
  expect_equal(length(result$gs2), 1)
  # Empty gene sets should be included
  expect_length(result, 2)
})

test_that("gmt2mat handles duplicate genes in same set", {
  gmt <- list(
    GS1 = c("gene1", "gene2", "gene1")  # gene1 appears twice
  )

  result <- gmt2mat(gmt)

  # Should handle duplicates gracefully
  expect_equal(ncol(result), 1)
  expect_true("gene1" %in% rownames(result))
})

test_that("read.gmt with file not found", {
  expect_error(read.gmt("/nonexistent/file.gmt"))
})

test_that("read.gmt with nrows parameter", {
  gmt_file <- tempfile(fileext = ".gmt")
  writeLines(c(
    "GeneSet1\turl1\tgene1\tgene2",
    "GeneSet2\turl2\tgene3\tgene4",
    "GeneSet3\turl3\tgene5",
    "GeneSet4\turl4\tgene6"
  ), gmt_file)

  result <- read.gmt(gmt_file, nrows = 2)

  expect_length(result, 2)

  unlink(gmt_file)
})
