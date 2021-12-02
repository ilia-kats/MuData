library(fs)

fileh5mu <- paste0(file_temp(), ".h5mu")
test_that("COAD example", {
    coad <- curatedTCGAData::curatedTCGAData("COAD", "*", version = "2.0.1", dry.run = FALSE)
    w <- capture_warnings(WriteH5MU(coad, fileh5mu))
    expect_match(w, "Objects of class RaggedExperiment are currently unsupported, skipping...", fixed=TRUE, all=TRUE)
})
