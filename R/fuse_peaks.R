#' @title Fuse several peaks together
#'
#' @description Function to fuse the indicated peaks into one. It is called by
#' [fuse_all_peaks] to perform each of the peak fusion procedures on a master
#' table.
#'
#' @param master.table A master table data frame, as obtained from
#' [build_master_table].
#'
#' @param peaks_to_fuse A character vector specifying the peaks to be fused.
#'
#' @import dplyr
#'
#' @export
fuse_peaks <- function(master.table, peaks_to_fuse){

  # Check if peaks_to_fuse is empty.
  if (length(peaks_to_fuse) == 0) {
    stop("Error: 'peaks_to_fuse' must not be empty.")
  }

  # Check if any of the specified peaks are not present in master.table.
  missing_peaks <- setdiff(peaks_to_fuse, unique(master.table$Peak))
  if (length(missing_peaks) > 0) {
    stop(paste("Error: The following peaks are not present in 'master.table':", paste(missing_peaks, collapse = ", ")))
  }

  # Check if any of the required columns are missing from master.table.
  required_columns <- c("Peak", "RI", "Compound")
  missing_columns <- setdiff(required_columns, colnames(master.table))
  if (length(missing_columns) > 0) {
    stop(paste("Error: The following columns are missing from 'master.table':", paste(missing_columns, collapse = ", ")))
  }

  # Create a new data frame called peaks_sum by filtering rows in
  # master.table where the value of Peak is in peaks_to_fuse, selecting all
  # columns except for those between (and including) Peak and Mod. position,
  # calculating column sums while ignoring missing values, transposing the
  # resulting vector, and converting it to a data frame.
  peaks_sum <- master.table |>
    filter(get("Peak") %in% peaks_to_fuse) |>
    select(!contains("Peak"):contains("Mod.position")) |>
    colSums(na.rm = T) |>
    t() |>
    as.data.frame()

  # Create an empty data frame called empty_peaks with the same number of
  # columns as peaks_sum and one fewer row than the length of peaks_to_fuse.
  empty_peaks <- matrix(NA
                        , ncol = length(peaks_sum)
                        , nrow = length(peaks_to_fuse) - 1) |>
    as.data.frame()

  # Set column names of empty_peaks to match those of peaks_sum.
  colnames(empty_peaks) <- colnames(peaks_sum)

  # Add rows from empty_peaks to the bottom of peaks_sum.
  peaks_sum <- rbind(peaks_sum
                     , empty_peaks)

  # Calculate median value in the RI column for rows where Peak is in
  # peaks_to_fuse.
  new_RI <- master.table |>
    filter(get("Peak") %in% peaks_to_fuse) |>
    pull("RI") |>
    stats::median() |>
    # Round RI and ensure that is an integer
    round(digits = 0) |>
    as.integer()

  # Create a new data frame called new_RI with only the RI column containing
  # new_RI values, repeated as many times as the length of peaks_to_fuse.
  new_RI <- data.frame("RI" = rep(new_RI
                                , length(peaks_to_fuse)))

  # Add columns from master.table and new_RI to peaks_sum.
  peaks_sum <- cbind(master.table |>
                       filter(get("Peak") %in% peaks_to_fuse) |>
                       select(contains("Peak")
                              , contains("Compound"):contains("Mod.position"))
                     , new_RI
                     , peaks_sum) |>
    relocate(contains("RI"), .before = contains("Compound")) |>
    as_tibble()

  # Concatenate the information of the compounds contained by the fused peaks,
  # if it differs, it avoids loosing that information.
  # It is triggered by a mismatch in the name of the compounds contained in
  # the peaks that are being fused
  if (length(unique(peaks_sum$Compound)) != 1) {
    peaks_sum$Compound <- peaks_sum$Compound |>
      paste(collapse = "|") |>
      rep(length(peaks_sum$Compound))

    if (length(unique(peaks_sum$Class)) != 1) {
      peaks_sum$Class <- peaks_sum$Class |>
        paste(collapse = "|") |>
        rep(length(peaks_sum$Class))
    }

    if (length(unique(peaks_sum$Mod.position)) != 1) {
      peaks_sum$Mod.position <- peaks_sum$Mod.position |>
        paste(collapse = "|") |>
        rep(length(peaks_sum$Mod.position))
    }
  }

  # Update rows in master_table where Peak matches values in peaks_sum
  master.table <- rows_update(master.table
                              , peaks_sum
                              , by = 'Peak')
  return(master.table)
}
