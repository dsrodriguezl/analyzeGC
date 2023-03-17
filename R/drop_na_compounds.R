#' @title Remove the peaks that do not contain compounds of interest
#'
#' @description
#' The function takes a list of data frames and removes peaks (rows) that do not
#' contain compounds of interest. Peaks are considered to not contain compounds
#' of interest if corresponding entry in the "Compound" column is NA.
#'
#' @param filtered_data
#' An aligned data set list, as obtained with [trace_comps].
#'
#' @import dplyr
#'
#' @examples
#' # Removing the NA compounds from a single data set
#' IW_filtered <- filtered_samples_list$`Winter_In-hive workers_A. m. mellifera`
#' drop_na_compounds(IW_filtered)
#'
#' # Removing NA compounds from sevreal data sets within a list
#' filtered_samples_list2 <- filtered_samples_list |>
#'   lapply(drop_na_compounds)
#'
#' @export
drop_na_compounds <- function(filtered_data) {
  # Use lapply to iterate through each data frame in filtered_data
  filtered_data <- lapply(filtered_data, function(df) {
    # Remove rows where Compound is NA
    df |> filter(!is.na(get("Compound")))
  })
  # Return modified list of data frames
  filtered_data
}
