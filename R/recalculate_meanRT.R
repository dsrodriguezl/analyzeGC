#' @title recalculate the mean RT of data set after alignment correction
#'
#' @param aligned_data Aligned data set as obtained with [correct_alignment]
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @examples
#'
#' # Recalculate the mean RT for a single data set
#' corrected_IW <- corrected_samples_list$`Winter_In-hive workers_A. m. mellifera`
#'
#' corrected_IW <-recalculate_meanRT(corrected_IW)
#'
#' # Recalculate the mean RT for several data sets within a list
#' corrected_samples_list2 <- lapply(corrected_samples_list, recalculate_meanRT)
#'
#' @export
recalculate_meanRT <- function(aligned_data) {

  # Function to remove rows with empty peaks
  drop_empty_peaks <- function(table) {
    if (length(colnames(table[colSums(table) == 0])) > 0) {
      table <- table |>
        select(-all_of(colnames(table[colSums(table) == 0])))
    }
    table
  }

  # Iterate through RT and Area data frames within aligned_data
  for (df_name in names(aligned_data)) {
    # cat('\n')
    # print(df_name)
    aligned_df <- aligned_data[[df_name]]

    # Remove empty peaks
    aligned_df <- drop_empty_peaks(aligned_df)

    # Overwrite old version of aligned_df inside aligned_data
    aligned_data[[df_name]] <- aligned_df
  }

  # Extract the RT data frame
  temp.RT.table <- aligned_data[["RT"]] |>
    t() |>
    as.data.frame()

  # Replace 0s with NAs
  temp.RT.table[temp.RT.table == 0] <- NA

  # Get a vector listing the samples in descending order, regarding their total
  # abundance
  samples_order <- aligned_data[["Area"]] |>
    rowSums() |>
    sort(decreasing = T) |>
    names()

  # Calculate the new mean RT for the given data frame and store it in the
  # mean.RT.list, with the corresponding name
  aligned_data[["RT"]] <- temp.RT.table |>
    mutate("Peak" = paste0("P", 1:nrow(temp.RT.table))
           , mean_RT = rowMeans(temp.RT.table
                                ,na.rm = TRUE) |>
             round(digits = 3)) |>
    select(contains("Peak"), contains("mean_RT")) |>
    bind_cols(temp.RT.table |>
                # Order the samples
                select(all_of(samples_order))) |>
    as_tibble()

  aligned_data[["Area"]] <- aligned_data[["RT"]] |>
    select(contains("Peak"), contains("mean_RT")) |>
    bind_cols(aligned_data[["Area"]] |>
                t() |>
                as.data.frame() |>
                # Order the samples
                select(all_of(samples_order)))

  aligned_data
}
