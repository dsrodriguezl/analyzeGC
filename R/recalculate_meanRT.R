#' @title recalculate the mean RT of a list of data frames after alignment
#' correction
#'
#' @param RT.list
#' List containing the RT aligned data frames after alignment correction with
#' correct_alignment.
#'
#' @param area.list
#' List containing the area aligned data frames after alignment correction with
#' correct_alignment.
#'
#' @param output.list
#'  A character string indicating which of the lists ("Area" or "RT") by the
#'  function, after adding to its data frames the newly calculated mean RTs.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @examples
#'
#' corrected_samples_list_area2 <-
#'   recalculate_meanRT(RT.list = corrected_samples_list_RT
#'                      , area.list = corrected_samples_list_area
#'                      , output.list = "Area")
#'
#' corrected_samples_list_RT2 <-
#'   recalculate_meanRT(RT.list = corrected_samples_list_RT
#'                      , area.list = corrected_samples_list_area
#'                      , output.list = "RT")
#'
#'
#' @export
recalculate_meanRT <- function(RT.list
                               , area.list
                               , output.list) {
  # Function to remove rows with empty peaks
  drop_empty_peaks <- function(table) {
    table <- table |>
      select(-all_of(colnames(table[colSums(table) == 0])))
    table
  }

  # Remove empty peaks from both area and RT data frames within the lists
  RT.list <- RT.list |>
    lapply(drop_empty_peaks)

  area.list <- area.list |>
    lapply(drop_empty_peaks)

  # Calculate the new mean RT for each data frame within the list
  mean.RT.list <- list()
  # Iterate through the data frames in the list
  for (RT.table in names(RT.list)) {
    # Extract the RT data frame in turn for the current iteration
    temp.RT.table <- RT.list[[RT.table]] |>
      t() |>
      as.data.frame()

    # Replace 0s with NAs
    temp.RT.table[temp.RT.table == 0] <- NA

    # Calculate the new mean RT for the given data frame and store it in the
    # mean.RT.list, with the corresponding name
    mean.RT.list[[RT.table]] <- temp.RT.table |>
      transmute(mean_RT = rowMeans(temp.RT.table
                                   ,na.rm = TRUE))
  }
  # mean.RT.list

  # Define which list will be returned by the function
  if (output.list == "Area") {
    list_2_return <- area.list
  }
  if (output.list == "RT") {
    list_2_return <- RT.list
  }

  # add new mean RT to the data frames in the list that will be returned
  for (table in names(list_2_return)) {
    # Extract the corresponding table form the list
    table.mean.RT <- mean.RT.list[[table]]

    # Shape the table
    tmp.table <- cbind.data.frame("Peak" = row.names(table.mean.RT)
                                  , "mean_RT" = table.mean.RT$mean_RT
                                  # Place samples as columns
                                  , list_2_return[[table]] |>
                                    t() |>
                                    as.data.frame() |>
                                    # Order them in descending order
                                    # , regarding their total abundance
                                    select(all_of(area.list[[table]] |>
                                                    rowSums() |>
                                                    sort(decreasing = T) |>
                                                    names()
                                                  )
                                           )
                                  ) |>
      # Add column with peak label
      mutate("Peak" = paste0("P", 1:ncol(list_2_return[[table]]))) |>
      # Place the peak label as the first column
      select(contains("Peak"), everything()) |>
      # turn the data frame into a tibble
      as_tibble()

    # replace the table in the list to return with its reshaped version
    list_2_return[[table]] <- tmp.table
  }
  list_2_return
}
