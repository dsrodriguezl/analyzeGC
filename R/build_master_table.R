#' @title Build a master table from several group tables
#'
#' @description Function to build a master table from a list of group tables.
#' It relies on the base::merge function to fuse the group tables together
#' into the master table.
#'
#' @param tables.list List of group tables.
#' The group tables must contain an integer RI column, indicating the retention
#' index of the peaks (rows).
#'
#' @import dplyr
#' @import purrr
#'
#' @export
build_master_table <- function(tables.list) {
  master.table <- tables.list |>
    reduce(merge
           , all = T
           , sort = F) |>
    arrange(get("RI"))
  master.table
}
