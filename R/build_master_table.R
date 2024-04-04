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
    # Nest data frames into two sublists, regarding their type (i.e. RT or Area)
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
    # Assemble the master table
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
    }) |>
    # The merging process could artificially duplicate some peaks within some
    # samples. The following code finds and remove such duplicates.
    (function(mt){
      ri_dup <- mt |>
        pluck("RT") |>
        filter(duplicated(RI)) |>
        pull(RI) |>
        unique()

      RT_df <- mt |>
        pluck("RT")

      RT_corrections <- lapply(ri_dup, function(ri) {
        df_ri <- RT_df |>
          filter(RI == ri) |>
          mutate_at(vars(-(Peak:Mod.position)), function(column) {
            ifelse(column |> duplicated()
                   , NA
                   , column)
          })
      }) |>
        reduce(merge
               , all = T
               , sort = F) |>
        as_tibble()

      RT_df <- RT_df |>
        rows_update(RT_corrections)

      Area_df <- mt |>
        pluck("Area")

      samples <- RT_df |>
        select(-(Peak:Mod.position)) |>
        colnames()

      Area_df[samples][is.na(RT_df[samples])] <- NA

      mt[["RT"]] <- RT_df
      mt[["Area"]] <- Area_df

      mt
    })()

  return(master.table)
}
