#' @title Remove rare compounds from a group table
#'
#' @description The function removes from a group table, the compounds that are
#'  rare to that group, regarding to a minimum frequency threshold.
#'  Any compound that has a frequency within the group table below the given
#'  threshold is erased from the data frame of the group table.
#'
#' @param group_table A group table data frame, as obtained with
#' [retrieve_group_tables].
#'
#' @param threshold A numeric value, between 0 and 1 (default = 0.5), indicating
#' the minimum frequency threshold, a compound should have within the group
#' table, to not be erased.
#'
#' @import dplyr
#' @import tidyr
#'
#' @export
group_frequency_filter <- function(group_table, threshold = 0.5){
  # Check that group_table is a data frame
  if (!is.data.frame(group_table)) {
    stop("Error: group_table must be a data frame")
  }

  # Check that threshold is a numeric value between 0 and 1
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("Error: threshold must be a numeric value between 0 and 1")
  }

  group_table <- group_table |>
    select(-contains("present"))

  group_table[group_table == 0] <- NA

  # Select the columns corresponding to the samples data
  data_cols <- group_table |>
    select(!contains("Peak"):contains("Mod.position"))

  # Calculate the proportion of NAs for each peak
  na_proportions <- rowSums(is.na(data_cols)) / ncol(data_cols)

  # Identify rare compounds
  rare_compounds <- group_table$Peak[1 - na_proportions < threshold]

  # Calculate frequency for rare compounds only
  rare_comp_freq <- 1 - na_proportions[1 - na_proportions < threshold]

  cat('\n')
  print(paste0("The following "
              , length(rare_compounds)
              , " compounds are rare (frequency < "
              , threshold
              , ") within the group"))
  print("Therefore, they will be deleted from the data set of that group")
  data.frame("compound" = rare_compounds, "frequency" = rare_comp_freq) |>
    print()
  cat('\n')

  # Erase, from the group data, the peaks containing compounds that are rare to the group
  group_table <- group_table |>
    filter(!get("Peak") %in% rare_compounds)

  group_table[is.na(group_table)] <- 0

  cat('\n')
  print(" All rare compounds were removed from the corresponding group data")
  cat('\n')

  # Return the modified group_table
  return(group_table)
}
