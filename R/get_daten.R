#' @title Extract the abundance data from  a finished master table
#'
#' @description A function to extract the data frame with the peak
#' abundance per sample from the final master table
#'
#' @param master.table A finished master table data frame
#'
#' @export
get_daten <- function(master.table){

  master.daten <- master.table |>
    select(-(contains("Peak"):contains("Mod.position"))) |>
    t() |>
    as.data.frame() |>
    magrittr::set_colnames(master.table$Peak)

  master.daten
}
