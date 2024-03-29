library(MultiAssayExperiment)
library(SummarizedExperiment)

mae1_file <- tempfile()
se1_file <- tempfile()
se2_file <- tempfile()

test_that("a MAE object written to H5MU can be read", {
  se <- SummarizedExperiment(matrix(rnorm(100), ncol = 5))
  rownames(colData(se)) <- paste("obs-", 1:ncol(se), sep = "")
  rownames(se) <- paste("var-", 1:nrow(se), sep = "")
  mae <- MultiAssayExperiment(list(x = se))

  writeH5MU(mae, mae1_file)

  expect_error(mae_ <- readH5MU(mae1_file), NA)
  # backed
  expect_error(mae_b <- readH5MU(mae1_file, backed = TRUE), NA)

  # Modality class is recovered
  expect_equal(class(mae_[[1]])[1], class(se)[1])
  expect_equal(class(mae_b[[1]])[1], class(se)[1])

  # MAE class is also recovered
  expect_equal(class(mae_)[1], class(mae)[1])
  expect_equal(class(mae_b)[1], class(mae)[1])
})

test_that("a SE object written to H5AD can be read", {
  se <- SummarizedExperiment(matrix(rnorm(100), ncol = 5))
  rownames(colData(se)) <- paste("obs-", 1:ncol(se), sep = "")
  rownames(se) <- paste("var-", 1:nrow(se), sep = "")
  metadata(se) <- list("origin" = "SummarizedExperiment")
  col1 <- 1:ncol(se)
  colData(se)$col1 <- col1

  writeH5AD(se, se1_file)

  expect_error(se_ <- readH5AD(se1_file), NA)
  # backed
  expect_error(se_b <- readH5AD(se1_file, backed = TRUE), NA)
  expect_true(inherits(assay(se_b), "DelayedArray"))

  # Simple colData is recovered
  expect_equal(colData(se_)$col1, col1)
  expect_equal(colData(se_b)$col1, col1)

  # Simple metadata is recovered
  expect_equal(metadata(se_)$origin, "SummarizedExperiment")
  expect_equal(metadata(se_b)$origin, "SummarizedExperiment")

  # Object class is recovered
  expect_equal(class(se_)[1], class(se)[1])
  expect_equal(class(se_b)[1], class(se)[1])
})

test_that("a SE object with a sparse matrix written to H5AD can be read", {
  x <- matrix(rnorm(100), ncol = 5)
  mask <- matrix(rbinom(100, 1, .5), ncol = 5)
  x[as.logical(mask)] <- 0
  for (mxfmt in c("dgCMatrix", "dgRMatrix")) {
    x_sp <- as(x, mxfmt)

    se <- SummarizedExperiment(x_sp)
    rownames(colData(se)) <- paste("obs-", 1:ncol(se), sep = "")
    rownames(se) <- paste("var-", 1:nrow(se), sep = "")

    writeH5AD(se, se2_file)

    expect_error(se_ <- readH5AD(se2_file), NA)
    expect_true(inherits(assay(se_), "dgCMatrix"))

    expect_error(se_b <- readH5AD(se2_file, backed = TRUE), NA)
    expect_true(inherits(assay(se_b), "DelayedArray"))
  }
})

test_that("Categoricals columns with NA are loaded correctly", {
  values <- c(0, -1)
  categories <- c('a')
  res <- convert_categoricals(values, categories)
  expect_equal(as.character(res[1]), 'a')
  expect_true(is.na(res[2]))
})

test_that("Extra levels in categoricals are ignored", {
  values <- c(0)
  categories <- c('a', 'b')
  res <- convert_categoricals(values, categories)
  expect_equal(levels(res), c('a'))
})

test_that("Extra levels in categoricals are ignored when NA are present", {
  values <- c(0, -1)
  categories <- c('a', 'b')
  res <- convert_categoricals(values, categories)
  expect_equal(levels(res), c('a'))
})