#' @title Import GC integration results from CSV  files
#'
#' @description A function to import the integration results data frames,
#' as obtained from MassHunter.
#'
#' @param data.path.list list containing the paths of the CSV files from which
#' to import the integration results of each sample.
#'
#' @param patterns_2_delete a character type vector, listing the text  patterns
#' to be removed from the CSV files' paths, to extract the names of the samples
#' from them.
#'
#' @returns A list of tibble data frames with the integration results
#' (RT: retention time of each peak, and Area: area under the curve of each
#' peak).
#' Each tibble on the list corresponds to the CSV file of an individual
#' sample.
#'
#' @import stringr
#' @import readr
#' @import dplyr
#'
#' @author Daniel S. Rodríguez-Leon <72925497+dsrodriguezl@users.noreply.github.com>
#'
#' @export

import_gcms_data <- function(data.path.list
                             , patterns_2_delete){

  names_list <- str_split(data.path.list
                          , "/"
                          , simplify = T) |>
    str_subset(".CSV") |>
    str_remove_all(patterns_2_delete)

  data.list <- data.path.list |>
    lapply(read_csv, comment = "#"
           , col_select = c("Center X"
                            , "Area")) |>
    rename(RT = "Center X")

  names(data.list) <- names_list
  data.list
}
