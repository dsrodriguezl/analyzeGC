#' @title Import GC integration results from CSV  files
#'
#' @description A function to import the integration results data frames,
#' as obtained from MassHunter.
#'
#' @param data.path.list list containing the paths of the CSV files from which
#' to import the integration results of each sample.
#'
#' @param patterns_2_delete a character type vector, listing the text  patterns
#' to be removed from the CSV files' names, to extract the names of the samples
#' from them. This step is performed via str_remove_all.
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
#' @examples
#'
#' samples_path_data <- list.files(path = system.file("extdata"
#'                                                    , package = "analyzeGC")
#'                                 #  Get all CSV files in the folder
#'                                 , pattern = ".CSV|.csv"
#'                                 , full.names = TRUE) |>
#'   # Do not include standards in the samples' list
#'   stringr::str_subset('STD', negate = TRUE)
#'
#' import_mh_data(samples_path_data, patterns_2_delete = "DR_")
#'
#' @export

import_mh_data <- function(data.path.list
                             , patterns_2_delete = " "){

  names_list <- str_split(data.path.list
                          , "/"
                          , simplify = T) |>
    str_subset(".CSV|.csv") |>
    str_remove_all(paste(patterns_2_delete
                    , ".CSV"
                    , ".csv"
                    , sep = "|"))

  data.list <- data.path.list |>
    lapply(read_csv, comment = "#"
           , col_select = c("Center X"
                            , "Area")
           , show_col_types = F) |>
    lapply(rename, "RT" = "Center X")

  names(data.list) <- names_list
  data.list
}
