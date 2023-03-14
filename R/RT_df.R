#' @title Extract aligned retention time data frames
#'
#' @description A function to obtain the retention time (RT) aligned data frame
#' from an aligned data set, as produced by align_chromatogram2.
#' It automatically labels each peak as P1 to Pn, where n  =  number of peaks
#' in the data set.
#'
#' @param aligned_data An aligned data set as produced by align_chromatogram2
#'
#' @examples
#'
#' # Extract area data frame from a single aligned data set
#' Winter_IW <-
#'   aligned_samples_data_list$`Winter_In-hive workers_A. m. mellifera`
#'
#' RT_winter_IW <- RT_df(Winter_IW)
#'
#' # Extract area data frame from several aligned data sets within a list
#' samples_list_RT <- aligned_samples_data_list |>
#'   lapply(RT_df)
#'
#' @export
RT_df <- function(aligned_data){
  if (length(aligned_data) == 2) {
    aligned_data_RT <- aligned_data$RT
  } else {
    aligned_data_RT <- aligned_data$aligned$RT
  }

  rownames(aligned_data_RT) <- paste0("P", 1:nrow(aligned_data_RT))
  aligned_data_RT
}
