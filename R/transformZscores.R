#' Description of the `transformZscores` S4 method
#'
#' @param object The input ClareExperiment object
#' @param base The base to use for the logarithmic transformation (default: 2)
#' @param offset The offset to use for the logarithmic transformation (default: 3)
#' @param ... other arguments passed
#'
#' @return The transformed ClareExperiment object
#'
#' @export
#' @examples
#' zscores_mat <- matrix(runif(16), nrow = 4)
#' sample_metadata <- data.frame(
#' sample_id = factor(c("sample1", "sample2", "sample3", "sample4")),
#' subject_id = factor(c("subject1", "subject1", "subject2", "subject2")),
#' collection_date = as.Date(
#' "2023-01-01", "2023-02-01", "2023-03-01", "20203-04-01"
#' )
#' )
#' peptide_metadata <- data.frame(
#' peptide_id = c("PV2T_005322","PV2T_013621",
#' "PV2T_013622", "PV2T_013623"),
#' species = c("speciesA", "speciesB", "speciesC", "speciesD")
#' )
#' clare_expt <- ClareExperimentFromZscores(
#' zscores = zscores_mat,
#' colData = sample_metadata,
#' rowData = peptide_metadata
#' )
#' se_trans <- transformZscores(clare_expt)
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom methods setMethod
setMethod("transformZscores", "ClareExperiment",
          function(object, base = 2, offset = 3, ...) {
            assayData <- assay(object, "zscores")
            transformed_data <- transform_zscores(assayData, base = base, offset = offset, ...)
            assay(object, "zscores_log2") <- transformed_data
            return(object)
          }
)

#' Description of the `transform_zscores` helper function
#'
#' @param zscores The input data to be transformed
#' @param base The base to use for the logarithmic transformation (default: 2)
#' @param offset The offset to use for the logarithmic transformation (default: 3)
#'
#' @return The transformed data
#'
transform_zscores <- function(zscores, base=2, offset=3){
  ###log transform data
  # base=2;offset=3
  # add a constant to all values (8)
  data1 <- base^offset+zscores
  # Set any values less than 1 equal to 1.
  data1[data1<1]=1
  # take the log and set the min to -3.
  trans_zscores<- log(data1,base)-offset
  return(trans_zscores)
}
