#' @title group samples within a list
#'
#' @description A function to Group samples gcms integration data frames based
#' on unique group labels, within the list containing the data frames.
#'
#' @param sample.info A data frame with a column
#' @param group.label
#' @param samples.data.list
#'
#' @import dplyr
#' @import purrr
#'
#' @export
mg_list <- function(sample.info, group.label, samples.data.list){
  mg_list <- list()
  for (i in unique(sample.info[[group.label]])) {
    tmp <- sample.info |>
      filter(get(group.label) == i) |>
      pull(`Individual`)
    tmp <- keep(samples.data.list
                , samples.data.list |>
                  names() %in% tmp) |>
      as.list()
    mg_list[[i]] <- tmp
  }
  samples.data.list <- mg_list |>
    # Keep only non-empty entries of the list
    compact()
  samples.data.list
}
