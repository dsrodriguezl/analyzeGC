#' @title Filter out compounds of low abundance within samples
#'
#' @description
#' The function filters out low-abundance compounds within samples
#' from an aligned data set list. The user specifies a threshold for
#' the minimum percentage of a sample that a peak must represent. The
#' function removes peaks below this threshold and updates the group
#' area and RT data accordingly. The mean retention time is recalculated
#' and added to the data.
#'
#' @param unfiltered_group An aligned data set list that includes comps.info.
#'
#' @param threshold A numeric threshold, indicating the minimum % of a sample a
#' peak must represent to remain in the data set.
#'
#' @import dplyr
#'
#' @export
trace_comps <- function(unfiltered_group, threshold) {

  group_area <- unfiltered_group[["Area"]]

  group_daten <- group_area |>
    select(!contains("Peak"):contains("Compound"))

  samples_names <- colnames(group_daten)

  cat('\n')
  ### Verify the total abundance per sample before deleting trace compounds
  print("Total abundance per sample before deleting trace compounds")
  print(group_daten |> colSums())

  ### Calculate the relative abundance (%) of each compound per sample
  group_daten_percent <- group_daten |> t() / rowSums(group_daten |> t()) * 100
  group_daten_percent <- group_daten_percent |>
    t() |>
    as.data.frame()

  # cat('\n')
  # ### Verify that the sum of all relative abundances per sample is exactly 100
  # print("Total relative abundance per sample before deleting trace compounds")
  # print(group_daten_percent |> colSums())

  ### Delete every peak of a sample that represent less
  ### than the specified % threshold of the sample
  group_daten[group_daten_percent < threshold] <- NA

  cat('\n')
  ### Verify the total abundance per sample after deleting trace compounds
  print("Total abundance per sample after deleting trace compounds")
  print(group_daten |> colSums(na.rm = T))
  cat('\n')

  # Modify group_area
  group_area <- group_area |>
    select(contains("Peak"):contains("Compound")) |>
    bind_cols(group_daten)

  # Modify group_RT
  group_RT <- unfiltered_group[["RT"]]

  group_RT[samples_names][is.na(group_area[samples_names])] <- NA

  # Filter out peaks that are no longer present in any sample
  group_area <- group_area |>
    filter(rowSums(group_daten, na.rm = T) > 0)

  group_RT <- group_RT |>
    filter(rowSums(group_daten, na.rm = T) > 0)

  group_comps_info <- unfiltered_group[["comps.info"]] |>
    filter(get("Peak") %in% group_area$Peak)

  # recalculate meanRT
  group_RT <- group_RT |>
    mutate("mean_RT" = rowMeans(group_RT |>
                                select(all_of(samples_names))
                              , na.rm = T) |>
             round(digits = 3))

  group_area$mean_RT <- group_RT$mean_RT

  group_comps_info$mean_RT <- group_RT$mean_RT

  filtered_group <- list("RT" = group_RT
                         , "Area" = group_area
                         , "comps.info" =  group_comps_info)

  filtered_group
}
