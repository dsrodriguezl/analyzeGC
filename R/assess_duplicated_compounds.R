#' @title Assess the presence of duplicated compounds among groups
#'
#' @description This function takes in a list of group tables extracted from a
#' master table and assesses the presence of duplicated compounds among groups.
#'
#' @param group.tables.list List of group tables extracted from the master table
#'  with [retrieve_group_tables].
#'
#' @param plot Logical value indicating whether to generate density distribution
#' plots for each compound (default: TRUE)
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr slice
#' @importFrom dplyr contains
#' @importFrom dplyr mutate
#' @importFrom dplyr everything
#' @importFrom dplyr arrange
#' @importFrom purrr pluck
#' @importFrom purrr reduce
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_label
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 labs
#' @importFrom ggdist stat_slab
#' @importFrom magrittr set_colnames
#'
#' @export
assess_duplicated_compounds <- function(group.tables.list, plot = T) {
  # Find duplicated compound names in the group tables

  ## If group.tables.list is a list of data frames
  if(group.tables.list |> pluck(1) |> is.data.frame()) {
    duplicated_compounds <- group.tables.list |>
      pluck(1) |>
      filter(duplicated(get("Compound"))) |>
      pull("Compound") |>
      unique()
  }

  ## If group.tables.list is a list of lists of data frames
  if(group.tables.list |> pluck(1) |> pluck(1) |> is.data.frame()) {
    duplicated_compounds <- group.tables.list |>
      pluck(1) |>
      lapply(filter, duplicated(get("Compound"))) |>
      lapply(select, "Compound") |>
      lapply(unique) |>
      reduce(merge, sort = F, all = T) |>
      pull("Compound")
  }

  # Empty list to register the presence data of each duplicated compound
  group_presence <- list()

  # Loop over each compound name in duplicated_compounds
  for (compound_name in duplicated_compounds) {


    compound_table <- group.tables.list |>
      lapply(function(group.table) {
        ## If group.table is a data frame
        if(group.table |> is.data.frame()) {
          group.table <- group.table |>
            # Filter rows where the "Compound" column is equal to the current
            # `compound_name`
            filter(get("Compound") == compound_name) |>
            # select columns from "Peak" to "present"
            select(contains("Peak"):contains("present")) |>
            # Add a new column called "count" that is equal to 1 if "present"
            # is TRUE, otherwise it is equal to 0
            mutate("count" = ifelse(get("present") == T, 1, 0)
                   , .keep = "unused")
        }

        ## If group.table is a list of data frames
        if(group.table |> pluck(1) |> is.data.frame()) {
          group.table <- group.table |>
            pluck("Area") |>
            # Filter rows where the "Compound" column is equal to the current
            # `compound_name`
            filter(get("Compound") == compound_name) |>
            # select columns from "Peak" to "present"
            select(contains("Peak"):contains("present")) |>
            # Add a new column called "count" that is equal to 1 if "present"
            # is TRUE, otherwise it is equal to 0
            mutate("count" = ifelse(get("present") == T, 1, 0)
                   , .keep = "unused")
        }

        group.table

      }) |>
      reduce(rbind) |>
      arrange("RI")

    # Loop over each unique value in the "Peak"column of compound_table
    for (peak in unique(compound_table$Peak)) {
      # Sum the values in the "count" column where "Peak" corresponds to that of
      # the current iteration (peak), and store the result as an integer in the
      # "count" column.
      compound_table["count"][compound_table["Peak"] == peak] <-
        compound_table |>
        filter(get("Peak") == peak) |>
        pull("count") |>
        sum() |>
        as.integer()
    }

    # Remove duplicated rows based on the "Peak" column
    compound_table <- compound_table |>
      filter(!duplicated(get("Peak")))

    # If plot is TRUE, generate a density distribution plot for the compound
    if (plot == T) {
      p <- compound_table |>
        ggplot(aes(y = get("RI"))) +
        stat_slab() +
        geom_label(aes(x = 0.5
                       , label = paste(paste(get("Peak")
                                             , get("RI")
                                             , sep = " | ")
                                       , get("count")
                                       , sep = ": "))
                   , fontface = "bold"
                   , alpha = 0.6) +
        scale_y_reverse(breaks = compound_table$RI) +
        theme_classic() +
        labs(x =  "density"
             , y = "RI"
             , title = paste0("Density distribution of the retention index"
                              , " of peaks containing "
                              , compound_name)
             , subtitle = paste0("For each peak "
                                 , "the following information is diplayed: "
                                 , "peak label | Retention index: "
                                 , "count of groups in which the peak "
                                 , "was detected"))
      print(p)
    }

    # Filter rows  where the "Compound"column is equal to compound_name and
    # select the "present" column
    compound_table <- group.tables.list |>
      (function(l) {
        ## If l is a list of data frames
        if(l |> pluck(1) |> is.data.frame()) {
          l <- l |>
            lapply(filter, get("Compound") == compound_name) |>
            lapply(select, contains("present")) |>
            # Combine all tables into one
            reduce(cbind) |>
            # Set column names to the names of group.tables.list
            set_colnames(names(l)) |>
            # Add new columns called "Peak" and "RI"
            mutate("Peak" = l |>
                     lapply(filter, get("Compound") == compound_name) |>
                     pluck(1) |>
                     pull("Peak")
                   , "RI" = l  |>
                     lapply(filter, get("Compound") == compound_name) |>
                     pluck(1) |>
                     pull("RI")) |>
            # Set columns order
            select(contains("Peak"), contains("RI"), everything())
        }

        ## If l is a list of lists of data frames
        if(l |> pluck(1) |> pluck(1) |> is.data.frame()) {
          l <- l |>
            lapply(pluck, "Area") |>
            lapply(filter, get("Compound") == compound_name) |>
            lapply(select, contains("present")) |>
            # Combine all tables into one
            reduce(cbind) |>
            # Set column names to the names of group.tables.list
            set_colnames(names(l)) |>
            # Add new columns called "Peak" and "RI"
            mutate("Peak" = l |>
                     lapply(pluck, "Area") |>
                     lapply(filter, get("Compound") == compound_name) |>
                     pluck(1) |>
                     pull("Peak")
                   , "RI" = l  |>
                     lapply(pluck, "Area") |>
                     lapply(filter, get("Compound") == compound_name) |>
                     pluck(1) |>
                     pull("RI")) |>
            # Set columns order
            select(contains("Peak"), contains("RI"), everything())
        }

        l
      })()

    # Set row names to the vlaues in the "Peak" column
    row.names(compound_table) <- compound_table$Peak

    # Remove "Peak" column and transpose the data frame
    compound_table <- compound_table |>
      select(-contains("Peak")) |>
      t() |> as.data.frame()

    # Store compound_table into group_presence, defining the key with
    # compound_name
    group_presence[[compound_name]] <- compound_table
  }
  return(group_presence)
}
