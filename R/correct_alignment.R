#' @title Correct the peaks alignment
#'
#' @description Function to correct alignment of peaks across an aligned data
#' set as indicated by instructions encoded within data frames' list(s).
#' It iterates through the data frame samples, following the given instructions,
#'  to add new peaks and displace values from one peak to another
#'
#' @param aligned_data
#' Aligned data set as obtained with [align_chromatograms2] function.
#'
#' @param new_peaks
#' List with instructions for adding empty peaks to the aligned data frames,
#'  if needed.
#' The default value is NULL, assuming no new peaks should be added.
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
#' @param peak_movements
#' List with instructions to displace the peak area/RT values within the sample
#' columns of the aligned data frame.
#'
#' The names of the list entries correspond to the name of the samples for which
#' the peak area/RT values should be displaced. The function only performs the
#' indicated peak displacements, if it finds the sample within the given data
#' set.
#' This facilitates the usage of the function on a list of aligned data sets
#' with [lapply], as the function will perform the correct peaks displacements
#' to the correct data set within the list.
#'
#' The entries of the list should correspond to data.frames/tibbles with two
#' columns ("peaks_origin" and "peaks_target").
#' peaks_origin indicates the peaks (e.g. P10, P12) that holds the values to be
#' displaced within the aligned data frame. peaks_target indicates the peak to
#' which which the value should be displaced along the column of the indicated
#' sample.
#'
#' @import dplyr
#' @import tibble
#' @importFrom rlang :=
#'
#' @examples
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
#' # Create a list to guide the displacement of peak values
#' ## Sample 350 is an In-hive worker
#' ## Sample 328 is an Out-hive worker
#'
#' peaks_movements_list <- list("350" = data.frame(peaks_origin =
#'                                                       c(paste0("P"
#'                                                                , c(106, 107
#'                                                                , 124)))
#'                                              , peaks_target =
#'                                              c(paste0("P"
#'                                                       , c(105, 106
#'                                                           , 123))))
#'                         , "328" = data.frame(peaks_origin =
#'                                                         c(paste0("P"
#'                                                                  , c(26, 35
#'                                                                      , 85
#'                                                                      , 124
#'                                                                      , 128)))
#'                                              , peaks_target =
#'                                                c(paste0("P"
#'                                                         , c(25, 34
#'                                                             , 84
#'                                                             , 123
#'                                                             , 127)))))
#' # Correct the alignment of a single aligned area/RT data set
#' IW <- aligned_samples_data_list$`Winter_In-hive workers_A. m. mellifera`
#' IW <- correct_alignment(aligned_data = IW
#'                              , peak_movements = peaks_movements_list)
#'
#' # Correct the alignment of several aligned area/RT data frames within a list
#' corrected_samples_list_area <- lapply(aligned_samples_data_list
#'                                       , correct_alignment
#'                                       , peak_movements = peaks_movements_list)
#'
#'
#' @export
correct_alignment <- function(aligned_data, new_peaks = NULL, peak_movements) {
  if (class(aligned_data) == "GCalign") {
    aligned_data <- aligned_data[["aligned"]]
  }

  for (df_name in names(aligned_data)) {
    cat('\n')
    print(df_name)
    aligned_df <- aligned_data[[df_name]]

    if ("mean_RT" %in% colnames(aligned_df)) {
      # rownames(aligned_df) <- paste0("P", 1:nrow(aligned_df))
      aligned_df <- aligned_df |>
        select(-contains("mean_RT")) |>
        t() |>
        as.data.frame()
    } else {
      aligned_df <- aligned_df |>
        t() |>
        as.data.frame()
    }

    # Verify that the provided alignment has more than one sample
    if (nrow(aligned_df) > 1) {
      # Add new peaks in case new_peaks has been provided (it is not NULL)
      if (!is.null(new_peaks)) {
        cat('\n')
        print("Adding new empty peaks")
        # Loop through new_peaks to add the corresponding new peaks to the
        # aligned data frame
        for (sample in names(new_peaks)) {

          if (sample %in% row.names(aligned_df)) {
            cat('\n')
            paste0("Adding empty peaks to table with reference sample ", sample) |>
              print()

            df <- new_peaks[[sample]]
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

      cat('\n')
      print("Moving peak values")
      # Loop iterating through peak_movements to displace the peak values
      # as instructed in the data frames of the list
      for (sample in names(peak_movements)) {
        if (sample %in% rownames(aligned_df)) {
          cat('\n')
          # Report which is the sample assigned to the current iteration
          paste("Sample:", sample, sep = " ") |>
            print()

          # Extract the df with the peaks' displacement instructions for the
          # corresponding sample, from the peak_movements list
          peaks_list <- peak_movements[[sample]] #|>
            # pull("peaks_list")

          # # Extract the vector listing the displacements to be performed on the
          # # peaks of the corresponding sample from the movements
          # movement_dirs <- peak_movements[[sample]] |>
          #   pull(movement_dirs)
          #
          # # Assemble data frame to guide alignment corrections
          # peaks_movement <- data.frame(Dir = movement_dirs
          #                              , Peaks = peaks_list)

          # Set iterations counting on 1
          p_count <- 1

          # Loop iterating through each peak to be displaced within the sample
          # to perform the displacement of its value
          for (p_origin in peaks_list$peaks_origin) {
            cat('\n')

            p_target <- peaks_list |>
              filter(get("peaks_origin") == p_origin) |>
              pull("peaks_target")

            # Report which is the peak assigned to the current iteration
            paste("Movement No.", paste0(p_count, ":")
                  , "Value in", p_origin, "will be moved to", p_target
                  , sep = " ") |>
              print()

            aligned_df[sample, p_target] <- aligned_df[sample, p_origin]
            aligned_df[sample, p_origin] <- 0

            # aligned_df <- move_one_peak(aligned_df
            #                             , Peak = p
            #                             , Dir = peaks_movement |>
            #                               filter(get("Peaks") == p) |>
            #                               pull("Dir")
            #                             , Sample = sample)
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
    } else {
      warning(paste("aligned_data contains only one sample!"
                    , "No alignment correction will be performed"))
      return(aligned_data)
    }
    aligned_data[[df_name]] <- aligned_df
  }
  class(aligned_data) <- "corrected_alignment"

  aligned_data <- recalculate_meanRT(aligned_data)
  print("The mean RT values have been corrected")

  aligned_data
}
