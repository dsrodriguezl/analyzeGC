#' @title
#'
#' @description
#'
#' @param aligned_df
#' @param empty.peaks
#'
#' @import tibble
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
