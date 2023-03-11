#' @title Extract aligned area data frames
#'
#' @description A function to obtain the area aligned data frame from an
#' aligned data set, as produced by align_chromatogram2.
#' It automatically labels each peak as P1 to Pn, where n  =  number of peaks
#' in the data set.
#'
#' @param aligned_data An aligned data set as produced by align_chromatogram2.
#'
#' @author
#' Daniel S. Rodríguez-Leon <72925497+dsrodriguezl@users.noreply.github.com>
#'
#' @export
area_df <- function(aligned_data){
  if (length(aligned_data) == 2) {
    aligned_data_area <- aligned_data$Area
  } else {
    aligned_data_area <- aligned_data$aligned$Area
  }

  rownames(aligned_data_area) <- paste0("P", 1:nrow(aligned_data_area))
  aligned_data_area
}
