#' @title Add empty peak rows to an aligned data frame
#'
#' @description
#' Function to create empty peaks within an aligned area/RT data frame, given a
#' set of instructions encoded within list of data frames.
#'
#' @param aligned_data
#' aligned data set, as obtained with [align_chromatograms2]
#'
#' @param empty.peaks
#' List with instructions for adding empty peaks to the data frames.
#'
#' The names of the list entries correspond to the name of a sample present
#' in the aligned data frame, to use as reference. The function only add a peak
#' to the data frame, if finds the reference sample within it.
#' This facilitates the usage of the function on a list of aligned data frames
#' with lapply, as the function will add the empty peaks to the correct aligned
#' data frame within the list.
#'
#' The entries of the list should correspond to data.frames/tibbles with two
#' columns ("position.reference" and "direction"). Position.reference
#' indicates where the empty peak should be added to the aligned data frame.
#' The position reference indicate a peak (e.g. P10) as a location reference
#' within the aligned data frame. Direction ("before" or "after") indicates
#' whether the empty peak should be created as the row before or after the
#' reference peak.
#'
#' @import tibble
#' @importFrom rlang :=
#'
#' @examples
#'
#'
#' # Create a list to guide the addition of empty peaks
#' ## Sample 335 is an In-hive worker
#' ## Sample 339 is an Out-hive worker
#'
#' library(tibble)
#' empty_peaks <- list("335" = tribble(~position.reference, ~direction,
#'                                      "P100", "after")
#'                     , "339" = tribble(~position.reference, ~direction,
#'                                      "P100", "before"))
#'
#' # Add empty peaks to a single aligned area/RT data frame
#' IW <- aligned_samples_data_list$`Winter_In-hive workers_A. m. mellifera`
#' OW <- aligned_samples_data_list$`Winter_Out-hive workers_A. m. mellifera`
#'
#' ## In-hive workers data frame
#' IW <- IW |>
#'   add_empty_peaks(empty.peaks = empty_peaks)
#'
#' ## Out-hive workers data frame
#' OW <- OW |>
#'   add_empty_peaks(empty.peaks = empty_peaks)
#'
#' # Add empty peaks to several aligned area/RT data frames within a list
#' aligned_samples_data_list <- aligned_samples_data_list |>
#'   lapply(add_empty_peaks
#'          , empty.peaks = empty_peaks)
#'
#'
#' @export
add_empty_peaks <- function(aligned_data, empty.peaks) {

  if (length(aligned_data) > 2) {
    aligned_data <- aligned_data[["aligned"]]
  } else {
    aligned_data <- aligned_data
  }

  # tmp_aligned_data <- list()
  for (df_name in names(aligned_data)) {
    cat('\n')
    print(df_name)
    aligned_df <- aligned_data[[df_name]]

    rownames(aligned_df) <- paste0("P", 1:nrow(aligned_df))

    if ("mean_RT" %in% colnames(aligned_df)) {
      aligned_df <- aligned_df |>
        select(-contains("mean_RT")) |>
        t() |>
        as.data.frame()
    }

    if (nrow(aligned_df) > 1) {
      # Loop iterating through samples for the alignment correction
      for (sample in names(empty.peaks)) {

        if (sample %in% row.names(aligned_df)) {
          cat('\n')
          paste0("Adding empty peaks to sample ", sample) |>
            print()

          df <- empty.peaks[[sample]]
          for (row in row.names(df) |> as.integer()) {
            paste0("Empty peak number ", row) |>
              print()

            new_peak <- paste(df$direction[row]
                              , df$position.reference[row]
                              , sep = "_")
            paste0("Name of new peak: ", new_peak) |>
              print()

            if (df$direction[row] == "before") {
              aligned_df <- aligned_df |>
                add_column("{new_peak}" := 0
                           , .before = df$position.reference[row])
            }
            if (df$direction[row] == "after") {
              aligned_df <- aligned_df |>
                add_column("{new_peak}" := 0
                           , .after = df$position.reference[row])
            }
          }
        }
      }
    }

    # tmp_aligned_data[[df_name]] <- aligned_df
    aligned_data[[df_name]] <- aligned_df
  }
  aligned_data
}
