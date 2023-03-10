#' @title Diagnostic
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
#' @export
diagnostic_heatmap <- function(df_area_norm, title) {
  heatmap_colors <- viridis::turbo(200)

  if (length(df_area_norm) > 2) {
    p <- gplots::heatmap.2(as.matrix(df_area_norm |> log1p())
                       , main = title
                       , srtCol = 90
                       , dendrogram = "row"
                       # Ensures that the columns are not ordered
                       , Colv = "NA"
                       , trace = "none"
                       , key.xlab = "log1p of relative abundance (%)"
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
