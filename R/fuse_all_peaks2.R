
#' @title Fuse peaks
#'
#' @description The function fuses peaks within a master table as specified by
#' a list.
#'
#' @param df A master table data frame, as obtained from
#' [build_master_table], OR an aligned data frame as obtained from
#' [add_comps_info]
#'
#' @param fusion.list List of data frames specifying the fusions to be performed.
#'
#' @param df_type type of data frame given as df; either "master.table" OR
#' "group.table".
#'
#' @import dplyr
#' @import tidyr
#'
#' @export
fuse_all_peaks2 <- function(df, fusion.list, df_type){

  # Handle case where fusion.list is empty
  if (length(fusion.list) == 0) {
    warning("fusion.list is empty; no fusions will be performed")
    return(df)
  }

  # Handle case where fusion.list has only one item
  if (length(fusion.list) == 1 ) {
    warning("fusion list has only one item; only one peak will be fused.")
    return(df)
  }

  if (df_type == "master.table") {master.table <- df}
  if (df_type == "group.table") {group.table <- df}

  # Initialize fusion count
  f_count <- 1

  # Loop through each item in the fusion list
  for (peaks_2_fuse in fusion.list) {
    cat('\n')

    # Report the current fusion number and peaks to be fused
    print(paste("Fusion No.", f_count, sep = " "))
    print(paste(length(peaks_2_fuse)
                ,"peaks to be fused:"
                , paste(peaks_2_fuse, collapse = " ")
                , sep = " "))
    cat('\n')

    if (df_type == "master.table") {
      fuse_peaks_mt <- function(master.table, peaks_to_fuse){

        # Check if peaks_to_fuse is empty.
        if (length(peaks_to_fuse) == 0) {
          stop("Error: 'peaks_to_fuse' must not be empty.")
        }

        # Check if any of the specified peaks are not present in master.table.
        missing_peaks <- setdiff(peaks_to_fuse, unique(master.table$Peak))
        if (length(missing_peaks) > 0) {
          stop(paste("Error: The following peaks are not present in 'master.table':"
                     , paste(missing_peaks, collapse = ", ")))
        }

        # Check if any of the required columns are missing from master.table.
        required_columns <- c("Peak", "RI", "Compound")
        missing_columns <- setdiff(required_columns, colnames(master.table))
        if (length(missing_columns) > 0) {
          stop(paste("Error: The following columns are missing from 'master.table':"
                     , paste(missing_columns, collapse = ", ")))
        }

        # Create a new data frame called peaks_sum by filtering rows in
        # master.table where the value of Peak is in peaks_to_fuse, selecting all
        # columns except for those between (and including) Peak and Mod. position,
        # calculating column sums while ignoring missing values, transposing the
        # resulting vector, and converting it to a data frame.
        peaks_sum <- master.table |>
          filter(get("Peak") %in% peaks_to_fuse) |>
          select(!contains("Peak"):contains("Mod.position")) |>
          colSums(na.rm = T) |>
          t() |>
          as.data.frame()

        # Create an empty data frame called empty_peaks with the same number of
        # columns as peaks_sum and one fewer row than the length of peaks_to_fuse.
        empty_peaks <- matrix(NA
                              , ncol = length(peaks_sum)
                              , nrow = length(peaks_to_fuse) - 1) |>
          as.data.frame()

        # Set column names of empty_peaks to match those of peaks_sum.
        colnames(empty_peaks) <- colnames(peaks_sum)

        # Add rows from empty_peaks to the bottom of peaks_sum.
        peaks_sum <- rbind(peaks_sum
                           , empty_peaks)

        # Calculate median value in the RI column for rows where Peak is in
        # peaks_to_fuse.
        new_RI <- master.table |>
          filter(get("Peak") %in% peaks_to_fuse) |>
          pull("RI") |>
          stats::median() |>
          # Round RI and ensure that is an integer
          round(digits = 0) |>
          as.integer()

        # Create a new data frame called new_RI with only the RI column containing
        # new_RI values, repeated as many times as the length of peaks_to_fuse.
        new_RI <- data.frame("RI" = rep(new_RI
                                        , length(peaks_to_fuse)))

        # Add columns from master.table and new_RI to peaks_sum.
        peaks_sum <- cbind(master.table |>
                             filter(get("Peak") %in% peaks_to_fuse) |>
                             select(contains("Peak")
                                    , contains("Compound"):contains("Mod.position"))
                           , new_RI
                           , peaks_sum) |>
          relocate(contains("RI"), .before = contains("Compound")) |>
          as_tibble()

        # Concatenate the information of the compounds contained by the fused peaks,
        # if it differs, it avoids loosing that information.
        # It is triggered by a mismatch in the name of the compounds contained in
        # the peaks that are being fused
        if (length(unique(peaks_sum$Compound)) != 1) {
          peaks_sum$Compound <- peaks_sum$Compound |>
            stats::na.omit() |>
            paste(collapse = "|") |>
            rep(length(peaks_sum$Compound))

          if (length(unique(peaks_sum$Class)) != 1) {
            peaks_sum$Class <- peaks_sum$Class |>
              stats::na.omit() |>
              paste(collapse = "|") |>
              rep(length(peaks_sum$Class))
          }

          if (length(unique(peaks_sum$Mod.position)) != 1) {
            peaks_sum$Mod.position <- peaks_sum$Mod.position |>
              stats::na.omit() |>
              paste(collapse = "|") |>
              rep(length(peaks_sum$Mod.position))
          }
        }

        # Update rows in master_table where Peak matches values in peaks_sum
        master.table <- rows_update(master.table
                                    , peaks_sum
                                    , by = 'Peak')
        return(master.table)
      }

      # Print the peaks status before fusion
      print("Peaks before fusion")
      master.table |> filter(get("Peak") %in% peaks_2_fuse) |>
        as.data.frame() |>  print()

      # Fuse peaks using fuse_peaks
      master.table <- fuse_peaks_mt(master.table = master.table
                                 , peaks_to_fuse = peaks_2_fuse)

      cat('\n')
      # Print the peaks status after fusion
      print("Peaks after fusion")
      master.table |> filter(get("Peak") %in% peaks_2_fuse) |>
        as.data.frame() |>  print()
    }

    if (df_type == "group.table") {
      fuse_peaks_gt <- function(group.table, peaks_to_fuse){

        # Check if peaks_to_fuse is empty.
        if (length(peaks_to_fuse) == 0) {
          stop("Error: 'peaks_to_fuse' must not be empty.")
        }

        # Modify Area data frame ----
        group_area <- group.table[["Area"]]

        samples_names <- colnames(group_area)

        # Check if any of the specified peaks are not present in group.table.
        missing_peaks <- setdiff(peaks_to_fuse, unique(group_area$Peak))
        if (length(missing_peaks) > 0) {
          stop(
            paste(
              "Error: The following peaks are not present in Area data frame:"
              , paste(missing_peaks, collapse = ", ")))
        }

        # Check if any of the required columns are missing from group.table.
        required_columns <- c("Peak", "mean_RT", "Compound")
        missing_columns <- setdiff(required_columns, colnames(group_area))
        if (length(missing_columns) > 0) {
          stop(
            paste(
              "Error: The following columns are missing from teh Area data frame:"
              , paste(missing_columns, collapse = ", ")))
        }

        # Create a new data frame called peaks_sum by filtering rows in
        # group.table where the value of Peak is in peaks_to_fuse, selecting all
        # columns except for those between (and including) Peak and Compound,
        # calculating column sums while ignoring missing values, transposing the
        # resulting vector, and converting it to a data frame.
        peaks_sum <- group_area |>
          filter(get("Peak") %in% peaks_to_fuse) |>
          select(!contains("Peak"):contains("Compound")) |>
          colSums(na.rm = T) |>
          t() |>
          as.data.frame()

        # Create an empty data frame called empty_peaks with the same number of
        # columns as peaks_sum and one fewer row than the length of peaks_to_fuse.
        empty_peaks <- matrix(NA
                              , ncol = length(peaks_sum)
                              , nrow = length(peaks_to_fuse) - 1) |>
          as.data.frame()

        # Set column names of empty_peaks to match those of peaks_sum.
        colnames(empty_peaks) <- colnames(peaks_sum)

        # Add rows from empty_peaks to the bottom of peaks_sum.
        peaks_sum <- rbind(peaks_sum
                           , empty_peaks)

        # Add columns from group.table and new_RI to peaks_sum.
        peaks_sum <- cbind(group_area |>
                             filter(get("Peak") %in% peaks_to_fuse) |>
                             select(contains("Peak"):contains("Compound"))
                           , peaks_sum) |>
          as_tibble()

        # Concatenate the information of the compounds contained by the fused
        # peaks, if it differs, it avoids loosing that information.
        # It is triggered by a mismatch in the name of the compounds contained
        # in the peaks that are being fused
        if (length(unique(peaks_sum$Compound)) != 1) {
          peaks_sum$Compound <- peaks_sum$Compound |>
            stats::na.omit() |>
            paste(collapse = "|") |>
            rep(length(peaks_sum$Compound))
        }

        # Update rows in group_area where Peak matches values in peaks_sum
        group_area <- rows_update(group_area
                                    , peaks_sum
                                    , by = 'Peak')

        # Modify RT data frame ----
        group_RT <- group.table[["RT"]]

        samples_names <- colnames(group_RT)

        # Check if any of the specified peaks are not present in group_RT.
        missing_peaks <- setdiff(peaks_to_fuse, unique(group_RT$Peak))
        if (length(missing_peaks) > 0) {
          stop(
            paste(
              "Error: The following peaks are not present in RT data frame:"
              , paste(missing_peaks, collapse = ", ")))
        }

        # Check if any of the required columns are missing from group_RT
        required_columns <- c("Peak", "mean_RT", "Compound")
        missing_columns <- setdiff(required_columns, colnames(group_RT))
        if (length(missing_columns) > 0) {
          stop(
            paste(
              "Error: The following columns are missing from the RT data frame:"
                     , paste(missing_columns, collapse = ", ")))
        }

        # Create a new data frame called peaks_sum by filtering rows in
        # group_RT where the value of Peak is in peaks_to_fuse, selecting all
        # columns except for those between (and including) Peak and Compound,
        # calculating column medians while ignoring missing values.
        peaks_sum <- group_RT |>
          filter(get("Peak") %in% peaks_to_fuse) |>
          select(!contains("Peak"):contains("Compound")) |>
          summarize(across(everything(), ~median(., na.rm = T)))

        # Create an empty data frame called empty_peaks with the same number of
        # columns as peaks_sum and one fewer row than the length of peaks_to_fuse.
        empty_peaks <- matrix(NA
                              , ncol = length(peaks_sum)
                              , nrow = length(peaks_to_fuse) - 1) |>
          as.data.frame()

        # Set column names of empty_peaks to match those of peaks_sum.
        colnames(empty_peaks) <- colnames(peaks_sum)

        # Add rows from empty_peaks to the bottom of peaks_sum.
        peaks_sum <- rbind(peaks_sum
                           , empty_peaks)

        # Add columns from group.table and new_RI to peaks_sum.
        peaks_sum <- cbind(group_RT |>
                             filter(get("Peak") %in% peaks_to_fuse) |>
                             select(contains("Peak"):contains("Compound"))
                           , peaks_sum) |>
          as_tibble()

        # Concatenate the information of the compounds contained by the fused
        # peaks, if it differs, it avoids loosing that information.
        # It is triggered by a mismatch in the name of the compounds contained
        # in the peaks that are being fused
        if (length(unique(peaks_sum$Compound)) != 1) {
          peaks_sum$Compound <- peaks_sum$Compound |>
            stats::na.omit() |>
            paste(collapse = "|") |>
            rep(length(peaks_sum$Compound))
        }

        peaks_sum <- peaks_sum |>
          mutate(mean_RT =
                   rowMeans(peaks_sum |>
                              select(!contains("Peak"):contains("Compound"))
                            , na.rm = T))

        # Update rows in group_area where Peak matches values in peaks_sum
        group_RT <- rows_update(group_RT
                                  , peaks_sum
                                  , by = 'Peak')

        # Modify comps.info data frame ----
        group_comps <- group.table[["comps.info"]]

        # Check if any of the specified peaks are not present in group_RT.
        missing_peaks <- setdiff(peaks_to_fuse, unique(group_comps$Peak))
        if (length(missing_peaks) > 0) {
          stop(
            paste(
              "Error: The following peaks are not present in comps.info data frame:"
              , paste(missing_peaks, collapse = ", ")))
        }

        # Check if any of the required columns are missing from group_RT
        required_columns <- c("Peak", "Compound", "mean_RT")
        missing_columns <- setdiff(required_columns, colnames(group_comps))
        if (length(missing_columns) > 0) {
          stop(
            paste(
              "Error: The following columns are missing from the comps.info data frame:"
              , paste(missing_columns, collapse = ", ")))
        }

        peaks_sum <- group_comps |>
          filter(get("Peak") %in% peaks_to_fuse) |>
          select(!all_of(c("Peak", "mean_RT"))) |>
          mutate(across(everything()
                        , function(x){
                          x <- x |>
                            stats::na.omit() |>
                            paste(collapse = "|")
                          x <- ifelse(x == "", NA, x)
                        })) |>
          bind_cols(group_comps |>
                      filter(get("Peak") %in% peaks_to_fuse) |>
                      select(all_of(c("Peak", "mean_RT")))) |>
          select(contains("Peak")
                 , contains("Compound")
                 , contains("mean_RT")
                 , everything())

        if ("Chain.length" %in% colnames(peaks_sum)) {
          peaks_sum <- peaks_sum |>
            mutate(Chain.length = get("Chain.length") |>
                     as.integer())
        }

        # Update rows in group_area where Peak matches values in peaks_sum
        group_comps <- rows_update(group_comps
                                , peaks_sum
                                , by = 'Peak')

        # Replace the data frames of the group.table object
        group_area <- group_area |>
          mutate(mean_RT = group_RT$mean_RT
                 , Compound = group_comps$Compound)

        group.table[["Area"]] <- group_area

        group_RT <- group_RT |>
          mutate(Compound = group_comps$Compound)

        group.table[["RT"]] <- group_RT

        group_comps <- group_comps |>
          mutate(mean_RT = group_RT$mean_RT)

        group.table[["comps.info"]] <- group_comps

        return(group.table)
      }

      # Print the peaks status before fusion
      print("Peaks before fusion")
      group.table |>
        lapply(function(x){
          x |>
            filter(get("Peak") %in% peaks_2_fuse) |>
            as.data.frame()
      }) |>  print()

      # Fuse peaks using fuse_peaks
      group.table <- fuse_peaks_gt(group.table = group.table
                                 , peaks_to_fuse = peaks_2_fuse)

      cat('\n')
      # Print the peaks status after fusion
      print("Peaks after fusion")
      group.table |>
        lapply(function(x){
          x |>
            filter(get("Peak") %in% peaks_2_fuse) |>
            as.data.frame()
        }) |>  print()

    }

    # Increment the fusion count
    f_count <- f_count + 1
    cat('\n')
  }

  # Report the number of fusions performed
  print(paste("Finished!"
              , f_count - 1
              , "fusions were performed"
              , sep = " "))

  if (df_type == "master.table") {
    # Remove empty peaks that were created during the peaks fusing procedure
    master.table <- master.table |>
      filter(master.table |>
               select(!contains("Peak"):contains("Mod.position")) |>
               rowSums(na.rm = T) > 0)

    #  Correct peak numbering and transform into a tibble
    df <- master.table |>
      # Correct the peaks numbering
      mutate("Peak" = paste0("P", 1:nrow(master.table))) |>
      as_tibble()
  }

  if (df_type == "group.table") {
    group.table <- group.table |>
      lapply(filter, !is.na(get("mean_RT")))

    df <- group.table |>
      lapply(mutate, Peak = paste0("P", 1:nrow(group.table$Area))) |>
      as_tibble()

  }

  # Report that the empty peaks were removed
  print(paste0("All empty peaks that were created during the fusion "
               , "of the peaks were removed"))

  df
}
