#' @title Extract aligned retention time data frames
#'
#' @description A function to obtain the retention time (RT) aligned data frame
#' from an aligned data set, as produced by [align_chromatograms2].
#' It automatically labels each peak as P1 to Pn, where n  =  number of peaks
#' in the data set.
#'
#' @param aligned_data
#' An aligned data set as produced by [align_chromatograms2]
#'
#' @examples
#'
#' # Extract area data frame from a single aligned data set
#' Winter_IW <-
#'   aligned_samples_data_list$`Winter_In-hive workers_A. m. mellifera`
#'
#' RT_winter_IW <- extract_RT(Winter_IW)
#'
#' # Extract area data frame from several aligned data sets within a list
#' samples_list_RT <- aligned_samples_data_list |>
#'   lapply(extract_RT)
#'
#' @export
extract_RT <- function(aligned_data){
  names_check <- c("RT", "Area")
  if (identical(names(aligned_data), names_check)) {
    aligned_data_RT <- aligned_data$RT
  } else {
    aligned_data_RT <- aligned_data$aligned$RT
  }

  # rownames(aligned_data_RT) <- paste0("P", 1:nrow(aligned_data_RT))
  aligned_data_RT
}
