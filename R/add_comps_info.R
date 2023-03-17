#' @title Add the compounds' information to the aligned data sets within a list
#'
#' @param samples.list List of aligned data sets as obtained with
#' [adjust_abundance]
#'
#' @param comps.info.list List of data frames with the compounds' information,
#' as obtained with [get_chc_info].
#'
#' @import dplyr
#' @import tibble

#' @export
add_comps_info <- function(samples.list, comps.info.list) {

  for (aligned_group in names(samples.list)) {

    aligned_group_list <- samples.list[[aligned_group]]

    for (data_table in names(aligned_group_list)) {
      samples.list[[aligned_group]][[data_table]] <-
        merge(x = comps.info.list[[aligned_group]] |>
                select(contains("Peak"):contains("mean_RT"))
              , y = aligned_group_list[[data_table]]
              , all = T
              , sort = F) |>
        as_tibble()

      samples.list[[aligned_group]][["comps.info"]] <-
        comps.info.list[[aligned_group]]
    }
  }
  samples.list
}
