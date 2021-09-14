library(MuData)
library(MultiAssayExperiment)

mae1_file <- tempfile()

test_that("a MAE object with a single matrix can be written to H5MU", {
  se <- SummarizedExperiment(matrix(rnorm(100), ncol = 5))
  rownames(colData(se)) <- paste("obs-", 1:ncol(se), sep = "")
  mae <- MultiAssayExperiment(list(x = se))

  expect_error(WriteH5MU(mae, mae1_file), NA)
})
