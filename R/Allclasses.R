#' An S4 class to represent a longitudinal highly-multiplexed serology experiment
#' @rdname ClareExperiment
#' @export ClareExperiment
#' @import SummarizedExperiment
#' @importFrom methods is new
#' @importFrom S4Vectors SimpleList
#' @importFrom lubridate is.Date as_date
#'
#'
#'
setClass("ClareExperiment",
  contains = "SummarizedExperiment"
)

## in set validity enforce:
## date of collection use: lubridate
## look into lubridate::interval-class (an S4 timespan class)
##  don't use ts since it requires equispaced sampling OR
## sequence of collection in the colData data frame (increasing integer values?)
## zscores in assays (a slot for qiime2 provenance data?)
## should constructor transform the z scores?
## later can extend this class to include model formulas
## need to define zscores accessor function called below.

setValidity("ClareExperiment", function(object) {
  if (!("zscores" %in% assayNames(object))) {
    return("the assays slot must contain a matrix named 'zscores'")
  }
  # if (!is.numeric(zscores(object)))
  #   return("the z score data provided is not numeric")
  if ("counts" %in% assayNames(object)) {
    message("A counts assay is present. The counts must be processed to obtain Z
            scores prior to clare analysis")
  }

  if (!("sample_id" %in% colnames(colData(object)) & "subject_id" %in% colnames(colData(object)) &
        ("collection_date" %in% colnames(colData(object)) | "collection_sequence" %in% colnames(colData(object))))) {
    return("the colData must contain columns named 'sample_id', 'subject_id', and either 'collection_date' or 'collection_sequence'")
  }

  if (!is.Date(colData(object)$collection_date)) {
    # try converting it to date format
    object@colData@listData$collection_date <- lubridate::parse_date_time(object@colData@listData$collection_date, orders = c("ymd", "mdy", "dmy"))
    message("A collection_date column is present but not formatted as date. Converting to date.")
  }

  TRUE
})


#' ClareExperiment object and constructors
#'
#' @param se A SummarizedExperiment object.
#' @param zscores A matrix of Z scores
#' @param colData A data frame of sample names and data. Number of rows must match the number of columns present in `zscores`
#' @param rowData A data frame of peptide id's and data. Number of rows must match number of rows in `zscores`
#' @param ... arguments provided to \code{SummarizedExperiment} including rowRanges and metadata.
#'
#' @return A ClareExperiment object
#' @export
#'
#' @examples
#'
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
#' new_expt <- ClareExperimentFromZscores(
#' zscores = zscores_mat,
#' colData = sample_metadata,
#' rowData = peptide_metadata
#' )
ClareExperiment <- function(se) {
  if (!is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object")
  }
  # TODO add a lot more checks
  object <- new("ClareExperiment", se)
  # save the package version (ref DESeq2 code); this doesn't work as the func metadata is not available
  # metadata(object)[["version"]] <- packageVersion("clare")
}

#' @rdname ClareExperiment
#' @export
ClareExperimentFromZscores <- function(zscores, colData, rowData, ...) {

  if ("Sequence name" %in% colnames(zscores)) {
    pepnames <- zscores$`Sequence name`
    zscores <- as.matrix(zscores[, -1, drop = FALSE])
    row.names(zscores) <- pepnames
  }
  if (!(ncol(zscores) == nrow(colData))) {
    stop(paste0("the number of columns in the input matrix and number of rows in the colData must be the same"))
  }

  ## TODO
  ## add more tests:
  ## do names match but out of order?
  ## handling for df vs matrix

  se <- SummarizedExperiment(assays = SimpleList(zscores = zscores), colData = colData, rowData = rowData, ...)

  ce_obj <- ClareExperiment(se)

  return(ce_obj)
}
