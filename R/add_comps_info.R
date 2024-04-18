#' @title Add the compounds' information to the aligned data sets within a list
#'
#' @param samples.list List of aligned data sets as obtained with
#' [adjust_abundance]
#'
#' @param comps.info.list List of data frames with the compounds' information,
#' as obtained with [get_hc_info].
#'
#' @import dplyr
#' @import tibble

#' @export
add_comps_info <- function(samples.list, comps.info.list) {

  # Check if the inputs are lists
  if (!is.list(samples.list) | !is.list(comps.info.list)) {
    stop("Both inputs must be lists.")
  }

  # Check if the lists are named
  if (is.null(names(samples.list)) | is.null(names(comps.info.list))) {
    stop("Both lists must be named.")
  }

  for (aligned_group in names(samples.list)) {

    # Check if the aligned_group exists in comps.info.list
    if (!aligned_group %in% names(comps.info.list)) {
      stop(paste("Aligned group"
                 , aligned_group
                 , "not found in comps.info.list."))
    }

    aligned_group_list <- samples.list[[aligned_group]]

    for (data_table in names(aligned_group_list)) {
      samples.list[[aligned_group]][[data_table]] <-
        merge(x = comps.info.list[[aligned_group]] |>
                select(contains("Peak"), contains("Compound"))
              , y = aligned_group_list[[data_table]]
              , all = T
              , sort = F) |>
        as_tibble() |>
        select(contains("Peak")
               , contains("mean_RT")
               , contains("Compound")
               , everything())
    }
    samples.list[[aligned_group]][["comps.info"]] <-
      comps.info.list[[aligned_group]] |>
      bind_cols(samples.list[[aligned_group]][["RT"]] |>
                  select(contains("mean_RT"))) |>
      select(contains("Peak")
             , contains("Compound")
             , contains("mean_RT")
             , everything())
  }
  samples.list
}
