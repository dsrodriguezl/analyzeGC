#' @title Correct the peaks alignment
#'
#' @description Function to correct alignment of peaks across a data frame, as
#' indicated by a set of instructions encoded within list of data frames.
#' It iterates through the data frame samples, following the given instructions,
#'  and calls move_one_peak to perform each peak displacement.
#'
#' @param aligned_data
#' Aligned data frame as obtained with the area_df or RT_df functions.
#'
#' @param movements_list
#' List with instructions to displace the peak area/RT values within the sample
#' columns of the aligned data frame.
#'
#' The names of the list entries correspond to the name of the samples for which
#' the peak area/RT values should be displaced. The function only performs the
#' indicated peak displacements, if it finds the sample within the given data
#' frame.
#' This facilitates the usage of the function on a list of aligned data frames
#' with lapply, as the function will perform the correct peaks displacements to
#' the correct data frame within the list.
#'
#' The entries of the list should correspond to data.frames/tibbles with two
#' columns ("peaks_list" and "movement_dirs").
#' peaks_list indicates the peaks (e.g. P10, P12) that holds the values to be
#' displaced within the aligned data frame. movement_dirs indicates the
#' direction ("down" or "up") in which the peak value should be displaced along
#' the column of the indicated sample.
#'
#' @import dplyr
#'
#' @export
correct_alignment <- function(aligned_data, movements_list) {
  if (nrow(aligned_data) > 1) {
    # Loop iterating through samples for the alignment correction
    for (sample in movements_list |> names()) {
      if (sample %in% rownames(aligned_data)) {
        cat('\n')
        # Report which is the sample assigned to the current iteration
        paste("Sample:", sample, sep = " ") |>
          print()

        # Extract the vector listing the peaks to be displaced,
        # within the corresponding sample, from the movements_list
        peaks_list <- movements_list[[sample]] |>
          pull("peaks_list")

        # Extract the vector listing the displacements to be performed on the
        # peaks of the corresponding sample from the movements_list
        movement_dirs <- movements_list[[sample]] |>
          pull(movement_dirs)

        # Assemble data frame to guide alignment corrections
        peaks_movement <- data.frame(Dir = movement_dirs
                                     , Peaks = peaks_list)

        # Set iterations counting on 1
        p_count <- 1

        # Loop iterating through each peak to be displaced within the sample
        # to perfomr the displacement of its value
        for (p in peaks_movement$Peaks) {
          cat('\n')

          # Report which is the peak assigned to the current iteration
          paste("Peak No.", p_count, sep = " ") |>
            print()

          aligned_data <- move_one_peak(aligned_data
                                       , Peak = p
                                       , Dir = peaks_movement |>
                                         filter(get("Peaks") == p) |>
                                         pull("Dir")
                                       , Sample = sample)
          p_count <- p_count + 1
          cat('\n')
        }
        print(paste("Finished! The alignment of"
                    , p_count - 1
                    , "peaks was corrected"
                    , sep = " "))
        cat('\n')
      }
    }
  }
  aligned_data
}
