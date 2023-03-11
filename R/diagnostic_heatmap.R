#' @title Plot to evaluate peak alignment
#'
#' @description
#' Heatmap + dendrogram of the chemical composition of samples within an aligned
#' data frame.
#'
#' It uses gplots::heatmap.2 to produce the plot.
#'
#' If the received data frames contains onyl on sample, it returns a heatmap,
#' using ggplot2
#'
#' @param df_area_norm
#' Normalized area data frame, as produced by area_norm
#'
#' @param title Character string with the text for the title of the plot
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#'
#' @author
#' Daniel S. Rodríguez-Leon <72925497+dsrodriguezl@users.noreply.github.com>
#'
#' @export
diagnostic_heatmap <- function(df_area_norm, title) {
  heatmap_colors <- viridis::turbo(200)

  if (length(df_area_norm) > 2) {
    p <- gplots::heatmap.2(as.matrix(df_area_norm |> log1p())
                       , main = title
                       , srtCol = 90
                       , strRow = 90
                       , offsetRow = 0.02
                       , offsetCol = 0.02
                       , dendrogram = "row"
                       # Ensures that the columns are not ordered
                       , Colv = "NA"
                       , trace = "none"
                       , key.xlab = "log(1 + x) of relative abundance (%)"
                       , col = heatmap_colors)
  } else {
    p <- df_area_norm |>
      mutate(sample = row.names(df_area_norm)) |>
      pivot_longer(cols = starts_with("P"),
                   names_to = "Peak",
                   values_to = "rel.abundance")  |>
      ggplot(aes(x = Peak, y = sample)) +
      geom_raster(aes(fill = log1p(rel.abundance))) +
      scale_fill_viridis_c(option = "turbo") +
      ggtitle(title) +
      theme_classic()
  }
  print(p)
}
