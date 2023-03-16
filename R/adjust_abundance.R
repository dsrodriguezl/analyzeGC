#' @title Adjust the abundance of the peaks within samples
#'
#' @param aligned_data Aligned data set, as obtained with [recalculate_meanRT]
#' @param std.info Data frame with the standards information, as obtained
#' from [shape_hcstd_info]
#'
#' @export
adjust_abundance <- function(aligned_data, std.info) {

  area_table <- aligned_data$Area

  adjust.area <- function(x, correction =  area.adjustment) (x / correction)

  area_table2 <- tibble()
  for (i in 1:(nrow(std.info) - 1)) {
    interval <- std.info$mean_RT[i:(i + 1)]

    if (i == 1) {
      area.adjustment <- std.info |>
        filter(mean_RT == min(interval)) |>
        pull(area_correction)

      tmp <- area_table |>
        filter(mean_RT <= min(interval)) |>
        mutate_at(area_table |>
                    select(-(Peak:mean_RT)) |>
                    colnames()
                  , adjust.area)
      area_table2 <-  area_table2 |>
        bind_rows(tmp)
    }

    area.adjustment <- std.info |>
      filter(mean_RT == max(interval)) |>
      pull(area_correction)

    tmp <- area_table |>
      filter(min(interval) < mean_RT & mean_RT <= max(interval)) |>
      mutate_at(area_table |>
                  select(-(Peak:mean_RT)) |>
                  colnames()
                , adjust.area)
    area_table2 <-  area_table2 |>
      bind_rows(tmp)

    if (i == (nrow(std.info) - 1)) {
      tmp <- area_table |>
        filter(max(interval) <  mean_RT) |>
        mutate_at(area_table |>
                    select(-(Peak:mean_RT)) |>
                    colnames()
                  , adjust.area)
      area_table2 <-  area_table2 |>
        bind_rows(tmp)
    }
  }
  corrected_area_table <- area_table2

  area_table <- area_table |>
    pivot_longer(-(Peak:mean_RT)
                 , names_to = "sample"
                 , values_to = "area")

  area_table2 <- area_table2 |>
    pivot_longer(-(Peak:mean_RT)
                 , names_to = "sample"
                 , values_to = "corrected_area")

  area_table <-  area_table |>
    bind_cols(corrected_area = area_table2$corrected_area)

  for (muestra in unique(area_table$sample)) {
    long_sample_table <- area_table |>
      pivot_longer(-(Peak:sample)
                   , names_to = "ab_type"
                   , values_to = "abundance") |>
      filter(sample == muestra)

    p <- long_sample_table |>
      ggplot(aes(y = abundance
                 , x = mean_RT)) +
      geom_col(color = "black") +
      facet_wrap(vars(ab_type)
                 , ncol = 1
                 , scales = "free_y") +
      theme_classic() +
      labs(title = paste0("Sample "
                          , muestra
                          , " before and after abundance adjustment"))
    print(p)
  }

  return(corrected_area_table)
}
