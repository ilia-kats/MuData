library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

mae1_file <- tempfile()
se1_file <- tempfile()
sce1_file <- tempfile()

# writeH5MU

test_that("a MAE object with a single matrix can be written to H5MU", {
  se <- SummarizedExperiment(matrix(rnorm(100), ncol = 5))
  rownames(colData(se)) <- paste("obs-", 1:ncol(se), sep = "")
  mae <- MultiAssayExperiment(list(x = se))

  expect_error(writeH5MU(mae, mae1_file), NA)
})

# writeH5AD

test_that("a SE object with a single matrix can be written to H5AD", {
  se <- SummarizedExperiment(matrix(rnorm(100), ncol = 5))
  rownames(colData(se)) <- paste("obs-", 1:ncol(se), sep = "")
  rownames(se) <- paste("var-", 1:nrow(se), sep = "")

  expect_error(writeH5AD(se, se1_file), NA)
})

test_that("a SCE object with a single matrix can be written to H5AD", {
  sce <- SingleCellExperiment(matrix(rnorm(100), ncol = 5))
  colnames(sce) <- paste("obs-", 1:ncol(sce), sep = "")
  rownames(sce) <- paste("var-", 1:nrow(sce), sep = "")

  expect_error(writeH5AD(sce, sce1_file), NA)
})

test_that("a SE object with a sparse matrix can be written to H5AD", {
  x <- matrix(rnorm(100), ncol = 5)
  mask <- matrix(rbinom(100, 1, .5), ncol = 5)
  x[as.logical(mask)] <- 0
  for (mxfmt in c("dgCMatrix", "dgRMatrix")) {
    x_sp <- as(x, mxfmt)

    se <- SummarizedExperiment(x_sp)
    rownames(colData(se)) <- paste("obs-", 1:ncol(se), sep = "")
    rownames(se) <- paste("var-", 1:nrow(se), sep = "")

    expect_error(writeH5AD(se, se1_file), NA)
  }
})
