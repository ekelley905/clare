## code to prepare `data.R` dataset goes here
zscores_mat <- matrix(runif(16), nrow = 4)

sample_metadata <- data.frame(
  samples = factor(c("sample1", "sample2", "sample3", "sample4")),
  collection_date = as.Date(
    "2023-01-01", "2023-02-01", "2023-03-01", "20203-04-01"
  )
)

peptide_metadata <- data.frame(
  peptide_id = c("PV2T_000001","PV2T_000002",
                 "PV2T_000003", "PV2T_000004"),
  species = c("speciesA", "speciesB", "speciesC", "speciesD")
)

new_expt <- ClareExperimentFromZscores(
  zscores = zscores_mat,
  colData = sample_metadata,
  rowData = peptide_metadata
)

usethis::use_data(data.R, overwrite = TRUE)
