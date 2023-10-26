#' @title Build a master table from several group tables
#'
#' @description Function to build a master table from a list of group tables.
#' It relies on the base::merge function to fuse the group tables together
#' into the master table.
#'
#' @param tables.list List of group tables, as obtained from
#' [shape_group_table].
#' The group tables must contain an integer RI column, indicating the retention
#' index of the peaks (rows).
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import purrr
#'
#' @export
build_master_table2 <- function(tables.list) {

  # Function to remove the "Peak" column from a data frame, if it exists
  remove_peak <- function(df) {
    if ("Peak" %in% colnames(df)) {
      df <- df |>
        select(-contains("Peak"))
    }
    df
  }

  # Apply the function to each data frame in the list
  tables.list <- tables.list |>
    lapply(lapply, remove_peak)


  # Merge all data frames in the list into one master table
  master.table <- tables.list |>
    (function(l){
      RT_list <- l |>
        lapply(function(sl){
          sl |>
            keep(names(sl) |>
                          str_detect("RT"))
        })
      area_list <- l |>
        lapply(function(sl){
          sl |>
            keep(names(sl) |>
                   str_detect("Area"))
        })
      list("RT" = RT_list
           , "Area" = area_list)
    })() |>
    lapply(function(l){
      df <- l |>
        reduce(merge
               , all = T
               , sort = F) |>
        arrange(across(contains("RI")))
      df |>
        magrittr::set_colnames(df |>
                                 colnames() |>
                                 str_remove_all("RT.|Area."))
    })

  # Add a "Peak" column and move it before "RI" column
  master.table <- master.table |>
    lapply(function(df){
      df |>
        mutate("Peak" =  paste0("P", seq_len(nrow(df)))) |>
        relocate(contains("Peak"), .before = contains("RI")) |>
        as_tibble()
    })

  return(master.table)
}
