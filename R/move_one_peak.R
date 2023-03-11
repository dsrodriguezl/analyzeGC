#' @title Move one peak up or down in an aligned data set
#'
#' @description Function for "manually" moving one peak of one sample one row
#' bellow/above on the aligned data frame. It is called by correct_alignment for
#' manual displacements of peak values within samples in an aligned data frame,
#' during the alignment correction process.
#'
#' @param df Aligned data frame as obtained with the area_df or RT_df functions.
#'
#' @param Sample Character string indicating the name of the sample within
#' which to displace the peak ara/RT value.
#'
#' @param Peak Character string indicating the peak label of the row containing
#' the value to be displaced.
#'
#' @param Dir character string indicating the direction in which the area/RT
#' value of should be displaced within the column of the indicated sample
#' ("up" or "down").
#'
#' @import dplyr
#' @import tidyr
#'
#' @author Daniel S. Rodríguez-Leon <72925497+dsrodriguezl@users.noreply.github.com>
#'
#' @export
move_one_peak <- function(df, Sample, Peak, Dir){
  # Ensure sample is type character
  Sample <- Sample |> as.character()

  print(paste("Target peak:", Peak, sep = " "))
  # Move the selected peak for the selected sample
  ## If the movemvent should be toward one row before
  if (Dir == "up") {
    # Define the column where the compound will be moved, based on Dir value
    new_col <- ncol(df |> select(1:all_of(Peak))) - 1 |> as.integer()

    # Display section of the table for the sample before modification
    print("sample table around target peak before corrective displacement")

    if ((length(df |> select(1:all_of(new_col))) + 2) > length(df)) {
      print("Target peak is last peak of the table")
      df |> select(all_of(new_col):(all_of(new_col) + 1)) |>
        filter(rownames(df) == Sample) |> print()
    } else {
      df |> select(all_of(new_col):(all_of(new_col) + 2)) |>
        filter(rownames(df) == Sample) |> print()
    }


    df[Sample, new_col] <- df[Sample, Peak]
    df[Sample, Peak] <- 0

    # Display section of the table for the sample after modification
    print("Sample table around target peak after corrective displacement")

    if ((length(df |> select(1:all_of(new_col))) + 2) > length(df)) {
      print("Target peak is last peak of the table")
      df |> select(all_of(new_col):(all_of(new_col) + 1)) |>
        filter(rownames(df) == Sample) |> print()
    } else {
      df |> select(all_of(new_col):(all_of(new_col) + 2)) |>
        filter(rownames(df) == Sample) |> print()
    }
  }
  ## If the movemvent should be toward one row after
  if (Dir == "down") {
    new_col <- ncol(df |> select(1:all_of(Peak))) + 1 |> as.integer()

    # Display section of the table for the sample before modification
    print("sample table around target peak before corrective displacement")
    df |> select((all_of(new_col) - 2):all_of(new_col)) |>
      filter(rownames(df) == Sample) |> print()

    df[Sample, new_col] <- df[Sample, Peak]
    df[Sample, Peak] <- 0

    # Display section of the table for the sample after modification
    print("Sample table around target peak after corrective displacement")
    df |> select((all_of(new_col) - 2):all_of(new_col)) |>
      filter(rownames(df) == Sample) |> print()
  }

  df
}
