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

# anndata writes pandas string columns and indices (the default with
# pandas 3) as nullable-string-array groups, and None with the null encoding
as_nullable_strings <- function(file, name, mask) {
  values <- as.vector(rhdf5::h5read(file, name))
  rhdf5::h5delete(file, name)
  rhdf5::h5createGroup(file, name)
  rhdf5::h5write(values, file, paste(name, "values", sep = "/"))
  rhdf5::h5write(as.integer(mask), file, paste(name, "mask", sep = "/"))
  h5 <- rhdf5::H5Fopen(file)
  group <- rhdf5::H5Gopen(h5, name)
  rhdf5::h5writeAttribute("nullable-string-array", group, "encoding-type")
  rhdf5::h5writeAttribute("0.1.0", group, "encoding-version")
  rhdf5::H5Gclose(group)
  rhdf5::H5Fclose(h5)
}

add_null <- function(file, name) {
  h5 <- rhdf5::H5Fopen(file)
  space <- rhdf5::H5Screate("H5S_NULL")
  dataset <- rhdf5::H5Dcreate(h5, name, "H5T_IEEE_F32LE", space)
  rhdf5::h5writeAttribute("null", dataset, "encoding-type")
  rhdf5::h5writeAttribute("0.1.0", dataset, "encoding-version")
  rhdf5::H5Dclose(dataset)
  rhdf5::H5Sclose(space)
  rhdf5::H5Fclose(h5)
}

test_that("nullable string arrays and null values are read", {
  se <- SummarizedExperiment(matrix(rnorm(20), ncol = 4))
  rownames(colData(se)) <- paste("obs-", 1:ncol(se), sep = "")
  rownames(se) <- paste("var-", 1:nrow(se), sep = "")
  rowData(se)$gene <- paste("gene-", 1:nrow(se), sep = "")
  metadata(se) <- list("origin" = "SummarizedExperiment")
  file <- tempfile(fileext = ".h5ad")
  writeH5AD(se, file)

  as_nullable_strings(file, "var/_index", rep(FALSE, nrow(se)))
  as_nullable_strings(file, "var/gene", c(FALSE, TRUE, FALSE, FALSE, FALSE))
  add_null(file, "uns/unset")

  expect_error(se_ <- readH5AD(file), NA)

  # Index and column are plain vectors; masked values are NA
  expect_equal(rownames(se_), rownames(se))
  expect_equal(rowData(se_)$gene, c("gene-1", NA, "gene-3", "gene-4", "gene-5"))
  expect_null(dim(rowData(se_)$gene))

  # None is read as NULL; other metadata is recovered
  expect_true("unset" %in% names(metadata(se_)))
  expect_null(metadata(se_)$unset)
  expect_equal(metadata(se_)$origin, "SummarizedExperiment")
})