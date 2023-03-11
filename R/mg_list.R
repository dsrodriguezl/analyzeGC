#' @title group samples within a list
#'
#' @description A function to Group samples gcms integration data frames based
#' on unique group labels, within the list containing the data frames.
#'
#' @param sample.info
#' A data frame indicating the groups to which the samples belong
#'
#' @param group.label
#' A character indicating the name of the column of sample.info that contains
#' the group names by which the samples should be grouped
#'
#' @param samples.data.list
#' List of data frames with integration results per sample, as produced by
#' import_gcms_data
#'
#' @import dplyr
#' @import purrr
#'
#' @author
#' Daniel S. Rodríguez-Leon <72925497+dsrodriguezl@users.noreply.github.com>
#'
#' @export
mg_list <- function(sample.info, group.label, samples.data.list){
  mg_list <- list()
  for (i in unique(sample.info[[group.label]])) {
    tmp <- sample.info |>
      filter(get(group.label) == i) |>
      pull("Individual")
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
