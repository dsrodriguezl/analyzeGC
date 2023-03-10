#' @title Normalize aligned area data frame
#'
#' @description
#' Function to normalize the aligned area data frame produced by
#' align_chromatograms2.
#'
#' It is a wrapper around GCalignR::norm_peaks.
#'
#' Its main purpose is to facilitate the usage of the function with lapply.
#' By default set the rt_col_name to "RT", and conc_col_name to "Area".
#'
#' Automatically labels each peak as P1 to Pn, where n  =  number of peaks
#' in the data set.
#'
#' @param aligned_data An aligned data set as produced by align_chromatogram2
#'
#' @export
area_norm <- function(aligned_data){
  if (length(aligned_data) > 2) {
    df_area_norm <- GCalignR::norm_peaks(data = aligned_data
                               , rt_col_name = "RT"
                               , conc_col_name = "Area") |>
      t() |> as.data.frame()
  } else {
    df <- aligned_data$Area

    df_area_norm <- df[2] /
      colSums(df[2]) *
      100
  }

  rownames(df_area_norm) <- paste0("P", 1:nrow(df_area_norm))
  df_area_norm <- df_area_norm |>
    t() |> as.data.frame()
  df_area_norm
}
