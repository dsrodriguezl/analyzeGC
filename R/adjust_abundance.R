#' @title Adjust abundance values within samples and display diagnostic plots
#'
#' @param aligned_data Aligned data set, as obtained with [recalculate_meanRT].
#'
#' @param std.info Data frame with the standards information, as obtained
#' from [shape_hcstd_info].
#' It must contains columns "mean_RT" and "area_correction".
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @examples
#'
#' # Adjusting the abundance for a single data set
#' IW_data <- corrected_samples_list2$`Winter_In-hive workers_A. m. mellifera`
#'   adjusted_IW <- adjust_area_and_plot(aligned_data = IW_data
#'                                         , std_info = std_info)
#'
#' # Adjusting the abundance for several data sets within a list
#' filtered_samples_list <- unfiltered_samples_list |>
#'   lapply(adjust_area_and_plot, std_info = std_info)
#'
#' @export
adjust_abundance <- function(aligned_data, std.info) {

  area_table <- aligned_data$Area

  # Define a function that takes an input value x and divides it by a correction
  # factor
  adjust.area <- function(x, correction =  area.adjustment) (x / correction)

  area_table2 <- tibble()

  # Loop over each row in std.info, except for the last row
  for (i in 1:(nrow(std.info) - 1)) {
    # Extract the mean_RT values defining the RT interval within which
    # to correct the abundance of peaks in the current iteration
    interval <- std.info$mean_RT[i:(i + 1)]

    # If this is the first iteration of the loop
    if (i == 1) {
      # Extract the area correction values for compounds that eluted before
      # the first standard
      area.adjustment <- std.info |>
        filter(get("mean_RT") == min(interval)) |>
        pull("area_correction")

      # Adjust the abundance of peaks, within every sample, with RT <= mean_Rt
      # of the first standard
      tmp <- area_table |>
        filter(get("mean_RT") <= min(interval)) |>
        mutate_at(area_table |>
                    select(-(contains("Peak"):contains("mean_RT"))) |>
                    colnames()
                  , adjust.area)

      # Add the adjusted rows to area_table2
      area_table2 <-  area_table2 |>
        bind_rows(tmp)
    }

    # Extract the area correction values for compounds that eluted before
    # the standard with the highest mean_RT in the current iteration
    # and after the standard with the lowest mean_RT in the current iteration
    area.adjustment <- std.info |>
      filter(get("mean_RT") == max(interval)) |>
      pull("area_correction")

    # Adjust the abundance of peaks, within every sample, with RT <= mean_Rt
    # of the standard with the highest mean_RT and RT > mean_RT with the lowest
    # mean_RT
    tmp <- area_table |>
      filter(min(interval) < get("mean_RT") &
               get("mean_RT") <= max(interval)) |>
      mutate_at(area_table |>
                  select(-(contains("Peak"):contains("mean_RT"))) |>
                  colnames()
                , adjust.area)

    # Add the adjusted rows to area_table2
    area_table2 <-  area_table2 |>
      bind_rows(tmp)

    # If this is the last iteration of the loop
    if (i == (nrow(std.info) - 1)) {
      # Adjust the abundance of peaks, within every sample, with RT > mean_Rt
      # of the last standard
      tmp <- area_table |>
        filter(max(interval) <  get("mean_RT")) |>
        mutate_at(area_table |>
                    select(-(contains("Peak"):get("mean_RT"))) |>
                    colnames()
                  , adjust.area)

      # Add the adjusted rows to area_table2
      area_table2 <-  area_table2 |>
        bind_rows(tmp)
    }
  }
  corrected_area_table <- area_table2

  area_table <- area_table |>
    pivot_longer(-(contains("Peak"):contains("mean_RT"))
                 , names_to = "sample"
                 , values_to = "area")

  area_table2 <- area_table2 |>
    pivot_longer(-(contains("Peak"):contains("mean_RT"))
                 , names_to = "sample"
                 , values_to = "corrected_area")

  area_table <-  area_table |>
    bind_cols(corrected_area = area_table2$corrected_area)

  # Loop over each sample
  for (muestra in unique(area_table$sample)) {
    long_sample_table <- area_table |>
      pivot_longer(-(contains("Peak"):contains("sample"))
                   , names_to = "ab_type"
                   , values_to = "abundance") |>
      filter(get("sample") == muestra)

    # Create a bar plot showing the abundance on y-axis and mean_RT on x-axis
    # Use facet the plot, to show the un-adjusted and adjusted versions of the
    # peaks' abundance within the sample defined by the current iteration
    p <- long_sample_table |>
      ggplot(aes(y = get("abundance")
                 , x = get("mean_RT"))) +
      geom_col(color = "black") +
      facet_wrap(vars(get("ab_type"))
                 , ncol = 1
                 , scales = "free_y") +
      theme_classic() +
      labs(title = paste0("Sample "
                          , muestra
                          , " before and after abundance adjustment"))
    print(p)
  }

  aligned_data[["Area"]] <- corrected_area_table

  return(aligned_data)
}
