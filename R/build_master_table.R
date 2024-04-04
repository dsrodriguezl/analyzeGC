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
build_master_table <- function(tables.list) {

  # Function to remove the "Peak" column from a data frame, if it exists
  remove_peak <- function(df) {
    if ("Peak" %in% colnames(df)) {
      df <- df |>
        select(-contains("Peak"))
    }
    df
  }

  if (tables.list |> pluck(1) |> (function(x) {
    is.list(x) &
      (!is_tibble(x) | !is.data.frame(x))
  })()) {
    list_type <- "nested"
  } else {
    list_type <- "simple"
  }

  # Apply the function to each data frame in the list
  if(list_type == "nested") {
    tables.list <- tables.list |>
      lapply(lapply, remove_peak)
  }

  if(list_type == "simple") {
    tables.list <- tables.list |>
      lapply(remove_peak)
  }

  # Merge all data frames in the list into one master table
  if(list_type == "nested") {
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
  }

  if(list_type == "simple") {
    master.table <- tables.list |>
      reduce(merge
             , all = T
             , sort = T) |>
      arrange(get("RI"))
  }

  # Add a "Peak" column and move it before "RI" column
  if(list_type == "nested") {
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
          filter(duplicated(get("RI"))) |>
          pull("RI") |>
          unique()

        if (ri_dup |> length() == 0) {
          master.table <- mt
          return(master.table)
        }

        RT_df <- mt |>
          pluck("RT")

        RT_corrections <- lapply(ri_dup, function(ri) {
          df_ri <- RT_df |>
            filter(get("RI") == ri) |>
            mutate_at(vars(-(contains("Peak"):contains("Mod.position")))
                      , function(column) {
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
          select(-(contains("Peak"):contains("Mod.position"))) |>
          colnames()

        Area_df[samples][is.na(RT_df[samples])] <- NA

        mt[["RT"]] <- RT_df
        mt[["Area"]] <- Area_df

        mt
      })()
  }

  if(list_type == "simple") {
    master.table <- master.table |>
      mutate("Peak" =  paste0("P", seq_len(nrow(master.table)))) |>
      relocate(contains("Peak"), .before = contains("RI")) |>
      as_tibble()
  }

  return(master.table)
}
