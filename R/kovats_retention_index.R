#' @title Calculate the Kováts retention index
#'
#' @description
#' This function calculates the retention indices for the peaks in a data set
#' using the Kováts method. The retention indices are calculated based on the
#' mean retention times of n-alkanes within the samples or in a standard mixture
#' in case of an n-alkane that is not present in the samples.
#'
#' @param filtered_data An aligned data set list containing data frames with
#' information about compounds and their peak areas and retention times, as
#' obtained with [drop_na_compounds] or [trace_comps].
#'
#' @param std.info A data frame containing information about the standards
#' that will be used in the calculation of retention indices, as obtained with
#' [shape_hcstd_info].
#'
#' @import dplyr
#'
#' @examples
#' # Calculate the retention index for a single data set
#' IW_filtered <- filtered_samples_list2$`Winter_In-hive workers_A. m. mellifera`
#' ri_IW <- kovats_retention_index(IW_filtered, std.info = std_info)
#'
#' # Calculate the retention index for several data sets within a list
#' ri_samples_list <- filtered_samples_list2 |>
#'   lapply(kovats_retention_index, std.info = std_info)
#'
#' @export
kovats_retention_index <- function(filtered_data, std.info) {
  comps_info <- filtered_data[["comps.info"]]

  tmp.std.info <- std.info |>
    filter(!get("Chain.length") %in%
             (comps_info |>
                filter(get("Class") == "Alkane") |>
                pull("Chain.length")))
  # print("Unnecessary standards were deleted")

  comps_info <- merge(tmp.std.info
                      , comps_info
                      , all = T
                      , sort = F) |>
    arrange("mean_RT") |>
    select(contains("Peak")
           , everything()) |>
    select((contains("Peak"):contains("Mod.position"))) |>
    as_tibble()
  # print("Standards were added to the group table")

  cat('\n')

  # New data frame to avoid modifications in the original
  comps_info2 <- comps_info

  # list of peaks containing compounds of interest
  target_peaks <- comps_info2 |>
    # In case the table contains compounds without ID (NAs),
    # it gives them the "NA" string as ID, to avoid them from being left out
    # in the coming filter step
    mutate(Class = ifelse(is.na(Class)
                          , "NA"
                          , Class)) |>
    filter(get("Class") != "STD") |>
    pull("Peak")

  print(paste("The RI of"
              , length(target_peaks)
              , "compounds will be calculated"
              , sep = " "))

  # Standards are considered alkanes
  comps_info2["Class"][comps_info2["Class"] == "STD"] <- "Alkane"

  # A data frame listing the n-alkanes
  n_alkanes <- comps_info2 |>
    filter(get("Class") == "Alkane") |>
    pull("Compound")

  # Empty vector to store the retention indexes
  ri_list <- c()

  # set the lap count on 1
  lap_count = 1

  # Loop iterating over target_peaks
  for (n in target_peaks) {
    # Data frame with the info of the compound
    # that is defined by current iteration
    current_compound <- comps_info2 |>
      filter(get("Peak") == n)

    # Report lap_count number and compound
    cat('\n')
    print(paste("Compound", lap_count, current_compound[, "Compound"], sep = " "))

    # Calculate the RI depending on whether the compound is or not an alkane
    if (!current_compound$Compound %in% n_alkanes) {

      # Report current compound name
      # and that it is not considered as an alkane
      print(paste("Compound"
                  , current_compound[, "Compound"]
                  , "is not an Alkane"
                  , sep = " "))

      # get rt of current compound
      rt <- current_compound |>
        pull("mean_RT")

      # Extract previous alkanes data
      prev_alkane <- comps_info2 |>
        filter(get("mean_RT") < rt) |>
        filter(get("Class") == "Alkane")

      # If there is an alkane before the target peak, RI is calculated.
      # if not, then RI is equal to lap_count
      if(prev_alkane |> nrow() > 0) {
        prev_alkane <- prev_alkane |>
          filter(get("Chain.length") == max(prev_alkane[["Chain.length"]]))

        # Extract next alkane data
        next.alka <- comps_info2 |>
          filter(get("mean_RT") > rt) |>
          filter(get("Class") == "Alkane")

        if(next.alka |> nrow() > 0) {
          next.alka <- next.alka |>
            filter(get("Chain.length") == min(next.alka[["Chain.length"]]))

          # Calculate the retention index for the target compound
          retention_index <- 100 *
            (((current_compound$mean_RT - prev_alkane$mean_RT) /
                (next.alka$mean_RT - prev_alkane$mean_RT)) +
               prev_alkane$Chain.length)

          # Round up the RI to make it an integer
          retention_index <- retention_index |> round(digits = 0)
        } else {
          retention_index <- ri_list[[length(ri_list)]] + rt
        }

        # Store the RI in ri_list
        ri_list <- c(ri_list, retention_index)

      } else {

        retention_index <- lap_count

        # Store the RI in ri_list
        ri_list <- c(ri_list, retention_index)
      }

      # Report RI
      print(paste("RI:", retention_index, sep = " "))
    } else {
      # Report current compound/peak name
      # and that it is considered as an alkane
      print(paste(current_compound[, "Compound"], "is an Alkane", sep = " "))

      # Calculate the retention index for the target compound
      retention_index <- current_compound$Chain.length * 100

      # Store the retention index in the list
      ri_list <- c(ri_list, retention_index)

      # Report RI
      print(paste("RI:", retention_index, sep = " "))
    }
    lap_count = lap_count + 1 # adjust the lap count
  }

  # Ensure that ri_list contains integers
  ri_list <- as.integer(ri_list)

  # print(ri_list)

  # Remove standards
  comps_info2 <- comps_info2 |>
    filter(!get("Compound") %in% (comps_info |>
                                    filter(get("Class") == "STD") |>
                                    pull("Compound")))

  # Add RI column with all the calculated RI to the comps_info data frame
  # replacing the mean_RT column
  comps_info <- comps_info2 |>
    select(-contains("mean_RT")) |>
    bind_cols("RI" = ri_list) |>
    mutate("Peak" = paste0("P", 1:nrow(comps_info2))) |>
    select(contains("Peak"), contains("RI"), everything())

  rm(comps_info2)

  cat('\n')
  print("The RI calculation of all compounds has finished!")
  cat('\n')

  filtered_data[["comps.info"]] <- comps_info

  filtered_data[["Area"]] <- filtered_data[["Area"]] |>
    bind_cols("RI" = ri_list) |>
    mutate("Peak" = comps_info$Peak) |>
    select(contains("Peak")
           , contains("RI")
           , everything()
           , -contains("mean_RT"))

  filtered_data[["RT"]] <- filtered_data[["RT"]] |>
    bind_cols("RI" = ri_list) |>
    mutate("Peak" = comps_info$Peak) |>
    select(contains("Peak")
           , contains("RI")
           , contains("mean_RT")
           , everything())

  print("The data set has been modified")

  return(filtered_data)
}


