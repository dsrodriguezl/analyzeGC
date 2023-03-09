#' @title group samples within a list
#'
#' @description A function to Group samples gcms integration data frames based
#' on unique group labels, within the list containing the data frames.
#'
#' @param sample.info
#' @param samples.data.list
#'
#' @import dplyr
#' @import purrr
#'
#' @export
mg_list <- function(sample.info, samples.data.list){
  mg_list <- list()
  for (i in unique(sample.info$group_label)) {
    tmp <- sample.info |>
      filter(group_label == i) |>
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
