
#' @title Fuse peaks to correct miss-alignemnets in a master table
#'
#' @description The function fuses peaks within a master table as specified by
#' a list. It uses a [fuse_peaks] to perform each fusion.
#'
#' @param master.table A master table data frame, as obtained from
#' [build_master_table].
#'
#' @param fusion.list List of data frames specifying the fusions to be performed.
#'
#' @export
fuse_all_peaks <- function(master.table, fusion.list){

  # Handle case where fusion_list is empty
  if (length(fusion_list) == 0) {
    warning("fusion_list is empty; no fusions will be performed")
    return(master_table)
  }

  # Handle case where fusion_list has only one item
  if (length(fusion_list) == 1 ) {
    warning("fusion list has only one item; only one peak will be fused.")
    return(master.table)
  }

  # Initialize fusion count
  f_count <- 1

  # Loop through each item in the fusion list
  for (i in fusion.list) {
    cat('\n')
    # Report the current fusion number and peaks to be fused
    print(paste("Fusion No.", f_count, sep = " "))
    print(paste(length(i)
                ,"peaks to be fused:"
                , paste(i, collapse = " ")
                , sep = " "))
    cat('\n')
    # Print the peaks status before fusion
    print("Peaks before fusion")
    master.table |> filter(Peak %in% i) |> print()

    # Fuse peaks using fuse_peaks
    master.table <- fuse_peaks(master.table = master.table, peaks_to_fuse = i)

    cat('\n')
    # Print the peaks status after fusion
    print("Peaks after fusion")
    master.table |> filter(Peak %in% i) |> print()

    # Increment the fusion count
    f_count <- f_count + 1
    cat('\n')
  }

  # Report the numbe rof fusions performed
  print(paste("Finished!"
              , f_count - 1
              , "fusions were performed"
              , sep = " "))

  # Remove empty peaks that were created during the peaks fusing procedure
  master.table <- master.table |>
    filter(master.table |>
             select(!Peak:RI) |>
             rowSums(na.rm = T) > 0)

  # Report that the empty peaks were removed
  print(paste0("All empty peaks that were created during the fusion "
               , "of the peaks were removed"))

  #  Correct peak numbering and transform into a tibble
  master.table <- master.table |>
    # Correct the peaks numbering
    mutate(Peak = paste0("P", 1:nrow(master.table))) |>
    as_tibble()

  return(master.table)
}
