#' @title Add empty peak rows to an aligned data frame
#'
#' @description
#' Function to create empty peaks within an aligned area/RT data frame, given a
#' set of instructions encoded within list of data frames.
#'
#' @param aligned_df
#' aligned area/RT data frame, as obtained with area_df or RT_df functions.
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
#'
#' @author
#' Daniel S. Rodríguez-Leon <72925497+dsrodriguezl@users.noreply.github.com>
#'
#' @export
add_empty_peaks <- function(aligned_df, empty.peaks) {

  if (nrow(aligned_df) > 1) {
    # Loop iterating through samples for the alignment correction
    for (sample in empty.peaks |> names()) {

      if (sample %in% row.names(aligned_df)) {
        cat('\n')
        print(sample)
        df <- empty.peaks[[sample]]
        for (row in row.names(df) |> as.integer()) {
          print(row)
          new_peak <- paste(df$direction[row]
                            , df$position.reference[row]
                            , sep = "_")
          print(new_peak)

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
  aligned_df
}
