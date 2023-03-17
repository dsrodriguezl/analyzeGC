#' @title Shape the final group table
#'
#' @param samples_plus_ri An aligned data set list containing data frames with
#' information about compounds and their peak areas and retention times, as
#' obtained with [kovats_retention_index].
#'
#' @import dplyr
#' @import tidyr

#' @export
shape_group_table <- function(samples_plus_ri) {
  # Extract data frame with abundance data
  area_table <- samples_plus_ri[["Area"]]

  # Extract data frame with the information of the compounds
  comps_info <- samples_plus_ri[["comps.info"]]

  # Transform 0s to NAs
  area_table[area_table == 0] <- NA

  new_table <- merge(x = comps_info
                     , y = area_table
                     , all = T
                     , sort = F) |>
    arrange("RI") |>
    select(-contains("Peak")) |>
    as_tibble()

  new_table
}
