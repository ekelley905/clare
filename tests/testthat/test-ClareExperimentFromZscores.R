
# check if an S4 object is output
test_that("Returns se object", {
  expect_s4_class(new_expt, "SummarizedExperiment")
})

# test that input zscores matrix is not empty
test_that("input zscores matrix is not empty", {
  zscores_mat <- matrix(runif(16), nrow = 4)

  sample_metadata <- data.frame(
    sample_id = factor(c("sample1", "sample2", "sample3", "sample4")),
    subject_id = factor(c("subject1", "subject1", "subject2", "subject2")),
    collection_date = as.Date(
      "2023-01-01", "2023-02-01", "2023-03-01", "20203-04-01"
    )
  )

  peptide_metadata <- data.frame(
    peptide_id = c(
      "PV2T_000001", "PV2T_000002",
      "PV2T_000003", "PV2T_000004"
    ),
    species = c("speciesA", "speciesB", "speciesC", "speciesD")
  )

  ce_obj <- ClareExperimentFromZscores(zscores_mat, colData = sample_metadata, rowData = peptide_metadata)
  expect_false(nrow(assay(ce_obj, "zscores")) == 0)
})

# test that input colData dataframe is not empty
test_that("input colData dataframe is not empty", {
  zscores_mat <- matrix(runif(16), nrow = 4)

  sample_metadata <- data.frame(
    sample_id = factor(c("sample1", "sample2", "sample3", "sample4")),
    subject_id = factor(c("subject1", "subject1", "subject2", "subject2")), collection_date = as.Date(
      "2023-01-01", "2023-02-01", "2023-03-01", "20203-04-01"
    )
  )

  peptide_metadata <- data.frame(
    peptide_id = c(
      "PV2T_000001", "PV2T_000002",
      "PV2T_000003", "PV2T_000004"
    ),
    species = c("speciesA", "speciesB", "speciesC", "speciesD")
  )

  ce_obj <- ClareExperimentFromZscores(zscores_mat, colData = sample_metadata, rowData = peptide_metadata)
  expect_false(nrow(colData(ce_obj)) == 0)
})

# test that input rowData dataframe is not empty
test_that("input rowData dataframe is not empty", {
  zscores_mat <- matrix(runif(16), nrow = 4)

  sample_metadata <- data.frame(
    sample_id = factor(c("sample1", "sample2", "sample3", "sample4")),
    subject_id = factor(c("subject1", "subject1", "subject2", "subject2")), collection_date = as.Date(
      "2023-01-01", "2023-02-01", "2023-03-01", "20203-04-01"
    )
  )

  peptide_metadata <- data.frame(
    peptide_id = c(
      "PV2T_000001", "PV2T_000002",
      "PV2T_000003", "PV2T_000004"
    ),
    species = c("speciesA", "speciesB", "speciesC", "speciesD")
  )

  ce_obj <- ClareExperimentFromZscores(zscores_mat, colData = sample_metadata, rowData = peptide_metadata)
  expect_false(nrow(rowData(ce_obj)) == 0)
})

# test that input zscores matrix has same number of columns as input colData dataframe has rows
test_that("input zscores matrix has same number of columns as input colData dataframe has rows", {
  zscores_mat <- matrix(runif(20), nrow = 4)

  sample_metadata <- data.frame(
    sample_id = factor(c("sample1", "sample2", "sample3", "sample4")),
    subject_id = factor(c("subject1", "subject1", "subject2", "subject2")), collection_date = as.Date(
      "2023-01-01", "2023-02-01", "2023-03-01", "20203-04-01"
    )
  )

  peptide_metadata <- data.frame(
    peptide_id = c(
      "PV2T_000001", "PV2T_000002",
      "PV2T_000003", "PV2T_000004"
    ),
    species = c("speciesA", "speciesB", "speciesC", "speciesD")
  )


  expect_error(
    ce_obj <- ClareExperimentFromZscores(zscores_mat, colData = sample_metadata, rowData = peptide_metadata),
    "the number of columns in the input matrix and number of rows in the colData must be the same"
  )
})

# test that input zscores matrix and input rowData dataframe have matching rownames
test_that("input zscores matrix and input rowData dataframe have matching rownames for peptides", {
  zscores_mat <- matrix(runif(16), nrow = 4, dimnames = list(c("PV2T_000001", "PV2T_000002", "PV2T_000003", "PV2T_000004"), c("sample1", "sample2", "sample3", "sample4")))

  sample_metadata <- data.frame(
    sample_id = factor(c("sample1", "sample2", "sample3", "sample4")),
    subject_id = factor(c("subject1", "subject1", "subject2", "subject2")), collection_date = as.Date(
      "2023-01-01", "2023-02-01", "2023-03-01", "20203-04-01"
    )
  )

  peptide_metadata <- data.frame(
    peptide_id = c(
      "PV2T_000001", "PV2T_000002",
      "PV2T_000003", "PV2T_000004"
    ),
    species = c("speciesA", "speciesB", "speciesC", "speciesD")
  )

  ce_obj <- ClareExperimentFromZscores(zscores_mat, colData = sample_metadata, rowData = peptide_metadata)
  expect_equal(rownames(ce_obj@assays@data@listData$zscores), rownames(rowData(ce_obj)))
})

# test that input zscores matrix and input rowData dataframe have matching rownames, even if out of order
test_that("input zscores matrix and input rowData dataframe have matching rownames, even if out of order", {
  zscores_mat <- matrix(runif(16), nrow = 4, dimnames = list(c("PV2T_000002", "PV2T_000001", "PV2T_000003", "PV2T_000004"), c("sample1", "sample2", "sample3", "sample4")))

  sample_metadata <- data.frame(
    sample_id = factor(c("sample1", "sample2", "sample3", "sample4")),
    subject_id = factor(c("subject1", "subject1", "subject2", "subject2")), collection_date = as.Date(
      "2023-01-01", "2023-02-01", "2023-03-01", "20203-04-01"
    )
  )

  peptide_metadata <- data.frame(
    peptide_id = c(
      "PV2T_000001", "PV2T_000002",
      "PV2T_000003", "PV2T_000004"
    ),
    species = c("speciesA", "speciesB", "speciesC", "speciesD")
  )

  ce_obj <- ClareExperimentFromZscores(zscores_mat, colData = sample_metadata, rowData = peptide_metadata)
  expect_equal(rownames(ce_obj@assays@data@listData$zscores), rownames(rowData(ce_obj)))
})
