#' @title Shape the final aligned table
#'
#' @description This function shapes the final aligned (group) table, before
#' assembling the master table. It applies removes trace compounds, removes or
#' label unidentified compounds (See drop_unidentified_comps below), and
#' calculates the Kováts retention index of the peaks.
#'
#' Trace compounds filter: Removes compounds that have low abundance within
#' samples, regarding an user specified threshold (See trace_comps_threshold
#' below).
#'
#' Kováts retention index: The retention indices for the peaks in a data set
#' are calculated using the Kováts method. This is performed implementing the
#' [kovats_retention_index] function under the hood, check its documentation for
#' more details.
#'
#' @param aligned_data An aligned data set list ()that includes comps.info.
#' comps.info must be a named item in the list corresponding to a data frame
#' (or tibble) with the information of the compounds. It must at least have the
#' columns Peak, Compound, mean_RT, and Class.
#' Check ?add_comps_info for more information.
#' If you are analyzing CHC's, it is recommended to generate the comps.info
#' data frame with [get_hc_info].
#'
#' @param trace_comps_threshold A numeric threshold, indicating the minimum %
#' of a sample a peak must represent to remain in the data set (default = 0).
#'
#' @param drop_unidentified_comps Logical value (default to TRUE) indicating
#' whether to remove the peaks with unidentified compounds, those with no entry
#' (NA) in the column Compound.
#' If it is set as FALSE, NAs in column Compound will be replaced by the string
#' "unidentified".
#'
#' @param std.info A data frame containing information about the standards
#' that will be used in the calculation of retention indices, as obtained with
#' [shape_hcstd_info].
#'
#' @import dplyr
#' @import tidyr
#' @import purrr

#' @export
shape_aligned_table <- function(aligned_data
                                , trace_comps_threshold = 0
                                , drop_unidentified_comps = T
                                , std.info) {

  # Filter out compounds of low abundance within samples
  ## The function filters out low-abundance compounds within samples
  ## from an aligned data set list. The user specifies a threshold for
  ## the minimum percentage of a sample that a peak must represent. The
  ## function removes peaks below this threshold and updates the group
  ## area and RT data accordingly. The mean retention time is recalculated
  ## and added to the data.
  trace_comps <- function(unfiltered_group # An aligned data set list including
                                           # comps.info
                          , threshold # A numeric threshold, indicating the
                                      # minimum % of a sample a peak must
                                      # represent to remain in the data set.
                          ) {

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

  filtered_data <- trace_comps(aligned_data
                               , threshold = trace_comps_threshold)

  if (drop_unidentified_comps == T) {
    filtered_data <- lapply(filtered_data, function(df) {
      # Remove rows where Compound is NA
      df |>
        filter(!is.na(get("Compound")))
    })
  } else {
    filtered_data <- lapply(filtered_data, function(df) {
      # Remove rows where Compound is NA
      df |>
        mutate(Compound = ifelse(is.na(get("Compound"))
                                 , "unidentified"
                                 , get("Compound")))
    })
  }

  samples_plus_ri <- kovats_retention_index(filtered_data, std.info)

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
    (function(l) {
      l <- l |>
        # discard_at("comps.info") |>
        lapply(function (table) {
          merge(x = comps_info
                , y = table
                , all = T
                , sort = F) |>
            arrange("RI") |>
            select(-contains("Peak")) |>
            as_tibble()
        })
      l[["comps.info"]] <- l[["comps.info"]] |>
        colnames()
      class(l) <- "aligned_table"
      l
    })()
}
