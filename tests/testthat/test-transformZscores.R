# Generate a ClareExperiment object with random z-score data
data <- matrix(rnorm(25), nrow=5)
rownames(data) <- paste0("peptide", 1:5)
colnames(data) <- paste0("sample", 1:5)
sample_metadata <- data.frame(
  sample_id = colnames(data),
  subject_id = paste0("subject", 1:5),
  collection_date = as.Date(
    c("2023-01-01", "2023-02-01", "2023-03-01", "2023-04-01", "2023-05-01")
  )
)
peptide_metadata <- data.frame(
  peptide_id = rownames(data),
  species = c("speciesA", "speciesB", "speciesC", "speciesD", "speciesE")
)
ce <- ClareExperimentFromZscores(zscores=data, colData=sample_metadata, rowData=peptide_metadata)

# Test that transformZscores correctly transforms the data using default parameters
test_that("transformZscores transforms data with default parameters", {
  ce_transformed <- transformZscores(ce)
  transformed_data <- assay(ce_transformed, "zscores_log2")
  expected_data <- log2(data + 8) - 3
  expect_equal(transformed_data, expected_data)
})

# Test that transformZscores correctly transforms the data using custom parameters
test_that("transformZscores transforms data with custom parameters", {
  ce_transformed <- transformZscores(ce, base=10, offset=1)
  transformed_data <- assay(ce_transformed, "zscores_log2")
  expected_data <- log10(data + 10) - 1
  expect_equal(transformed_data, expected_data)
})

# Test that transformZscores returns a ClareExperiment object
test_that("transformZscores returns a ClareExperiment object", {
  ce_transformed <- transformZscores(ce)
  expect_true(class(ce_transformed) == "ClareExperiment")
})

# Test that transformZscores does not modify the input object
test_that("transformZscores does not modify the input object", {
  ce_copy <- ce
  transformZscores(ce_copy)
  expect_equal(ce, ce_copy)
})
