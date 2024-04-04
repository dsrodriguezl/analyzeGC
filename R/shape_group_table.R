#' @title Shape the final group table
#'
#' @param samples_plus_ri An aligned data set list containing data frames with
#' information about compounds and their peak areas and retention times, as
#' obtained with [kovats_retention_index].
#'
#' @import dplyr
#' @import tidyr
#' @import purrr

#' @export
shape_group_table2 <- function(samples_plus_ri) {
  # Extract data frame with abundance data
  area_table <- samples_plus_ri[["Area"]]

  # Extract data frame with RT data
  RT_table <- samples_plus_ri[["RT"]] |>
    select(-contains("mean_RT"))

  # Extract data frame with the information of the compounds
  comps_info <- samples_plus_ri[["comps.info"]]

  # Transform 0s to NAs
  area_table[area_table == 0] <- NA
  RT_table[][is.na(area_table)] <- NA

  samples_plus_ri[["Area"]] <- area_table
  samples_plus_ri[["RT"]] <- RT_table

  samples_plus_ri |>
    discard_at("comps.info") |>
    lapply(function (table) {
      merge(x = comps_info
            , y = table
            , all = T
            , sort = F) |>
        arrange("RI") |>
        select(-contains("Peak")) |>
        as_tibble()
    })
}
