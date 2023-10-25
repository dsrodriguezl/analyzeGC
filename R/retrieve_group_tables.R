#' @title Retrieve the group tables from the master table
#'
#' @param group.label Character indicating the column of grouping.info
#' that contains the labels of the groups that define the group tables.
#' It is meant to be use din a pipe line with [assess_duplicated_compounds] or
#' [group_frequency_filter].
#'
#' @param master.table master table data frames, as obtained with
#' [build_master_table].
#'
#' @param grouping.info Data frame with information of the samples, within which
#' the group labels can be found.
#'
#' @import dplyr
#' @import tidyr
#'
#' @export
retrieve_group_tables <- function(group.label
                                 , master.table
                                 , grouping.info){
  # Vector listing all unique group labels
  group_labels <- grouping.info[[group.label]] |>
    unique()

  # Vector listing the name of columns with information of the compounds
  comps.vars <- master.table |>
    select(-all_of(grouping.info$Individual)) |>
    colnames()

  group.tables.list <- list()
  for (group in group_labels) {
   samples <- grouping.info |>
     filter(get(group.label) == group) |>
      mutate("Individual" = as.character(get("Individual"))) |>
      pull("Individual")

   group.table <- master.table |>
     select(all_of(c(comps.vars, samples)))

   group.table[is.na(group.table)] <- 0

   group.table <- group.table |>
     mutate("present" = ifelse(group.table |>
                               select(-all_of(comps.vars)) |>
                               rowSums() > 0
                             , T
                             , F)) |>
     relocate(contains("present"), .after = all_of(comps.vars))

   group.tables.list[[group]] <- group.table
  }
  return(group.tables.list)
}
