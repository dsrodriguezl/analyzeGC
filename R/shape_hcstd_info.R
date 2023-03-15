#' @title Shape the standards information data frame
#'
#' @description Function to shape the data frame holding the compound
#' information for n-alkane standards
#'
#' @param comps_id.STD
#' data frame with the compound information for the n-alkane standards.
#'
#' The data frame should have a "Compound" and a "mean_RT" columns, plus a
#' column for every standard run.
#'
#' The Compound column should have the name of the n-alkane standards with the
#' following format: C"numeric chain length"_ane_NA.
#' The rest of the data frame, is just as it can be obtained via the RT_df
#' function, or after alignment correction.
#'
#' @param aligned_std
#' An object holding the aligned standards, as obtained from
#' [align_chromatograms2], or [recalculate_meanRT] in case their alignment
#' required to be corrected.
#'
#' @param short_std_pattern
#' Character string pattern to identify the short chain-length standards (until
#' C20). The function relies on the dplyr::starts_with selection helper to do
#' so.
#'
#' @param long_std_pattern
#' Character string pattern to identify the long chain-length standards (from
#' C21). The function relies on the dplyr::starts_with selection helper to do
#' so.
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#'
#' @examples
#'
#' std_info <- shape_hcstd_info(comps_id.STD = comps_id_std
#'                              , aligned_std = aligned_standards
#'                              , short_std_pattern = "L"
#'                              , long_std_pattern = "H")
#'
#' @export
shape_hcstd_info <- function(comps_id.STD
                              , aligned_std
                              , short_std_pattern
                              , long_std_pattern) {

  if (length(aligned_std) > 2) {
    aligned_std <- aligned_std[["aligned"]]
  }
  std_area <- aligned_std[["Area"]]

  # Temporal data frame
  std_df <- comps_id.STD

  # Extract low standards
  short_std <- std_df |>
    select(contains("mean_RT"), starts_with(short_std_pattern))

  # Extract high standards
  long_std <- std_df |>
    select(contains("mean_RT"), starts_with(long_std_pattern))

  # End of low standards
  short_end <- std_df |>
    filter(get("Compound") == "C20_ane_NA") |>
    pull("mean_RT")

  # Beginning of high standards
  long_start <- std_df |>
    filter(get("Compound") == "C21_ane_NA") |>
    pull("mean_RT")

  # Erase RT entries of peaks that are out of the range of each standards' set
  ## Low
  for (col in colnames(short_std |> select(-contains("mean_RT")))) {
    short_std[col][short_std["mean_RT"] > short_end] <- NA
  }

  ## High
  for (col in colnames(long_std |> select(-contains("mean_RT")))) {
    long_std[col][long_std["mean_RT"] < long_start] <- NA
  }

  # Replace the columns of the standard runs with the corrected versions
  std_df <- std_df |>
    select(contains("mean_RT"):contains("Compound")) |>
    bind_cols(long_std |> select(-contains("mean_RT"))
              , short_std |> select(-contains("mean_RT")))

  # Calculate the correct mean RT
  std_df <- std_df |>
    mutate("mean_RT" =
             rowMeans(std_df |>
                        select(-(contains("mean_RT"):contains("Compound")))
                      , na.rm = T)) |>
    #  Place Compound as the first column
    select(contains("Compound"), everything())
  std_df

  # Create std.info data frame
  std.info <- cbind(std_df
                    # Extract from the names the chain length, class,
                    # and the position of any given unsaturation.
                    # The later will be full of NAs as the standards
                    # are all n-alkanes
                    , data.frame(row.names =
                                   rownames(comps_id.STD)
                                 , t(as.data.frame(
                                   strsplit(comps_id.STD$Compound
                                            , "_"))))) |>
    as_tibble()

  # Define columns names
  colnames(std.info) <-
    c(colnames(std.info)[1:(length(colnames(std.info)) - 3)]
      , "Chain.length", "Class", "Mod.position")

  # Correct entries format
  ## The peaks/compounds within the standard runs, must be differentiated from
  ## the alkanes in the samples
  ## Class
  std.info['Class'][std.info['Class'] == "ane"] <- "STD"
  std.info

  ## Compound names
  # Store the compound names in a vector, where they will be altered into their
  # final format
  std_names <- std.info$Compound
  # Set the name of the compound inside the standard samples
  std_names[!is.na(std_names)] <- paste(std.info |>
                                          filter(!is.na(get(
                                            "Compound"))) |>
                                          pull("Class")
                                        , std.info |>
                                          filter(!is.na(get(
                                            "Compound"))) |>
                                          pull("Chain.length")
                                        , sep = "-")
  # Change the compound names in the std.info data frame to the correct format
  # names
  std.info$Compound <- std_names
  std.info

  ## Chain length
  # It needs to be defined as an integer
  std.info$Chain.length <- std.info$Chain.length |>
    str_remove("C") |>
    as.integer()

  std_df <- std.info |>
    select(starts_with(long_std_pattern)
           , starts_with(short_std_pattern))

  for (col in colnames(std_df)) {
    std_df[col][std_df[col] == 0] <- NA
    std_area[col][is.na(std_df[col])] <- NA
    std.info[col] <- std_area[col]
  }

  columns_2_omit <- std.info |>
    select(-all_of(colnames(std_df))) |>
    colnames()
  std.info <- std.info |>
    mutate("area" = rowMeans(std.info |>
                             select(-all_of(columns_2_omit))
                           , na.rm = T))
  std.info["mean_RT"][std.info["mean_RT"] == 0] <- NA
  std.info["area"][std.info["area"] == 0] <- NA

  std.info <- std.info |>
    relocate("area", .after = "mean_RT") |>
    select(contains("Compound"):contains("area")
           , contains("Chain.length"):contains("Mod.position")) |>
    filter(!is.na(get("Compound")))

  std.info <- std.info |>
    mutate("median_area" = stats::median(get("area"))
           , "area_correction" = get("area") / get("median_area")
           , "corrected_area" = get("area") / get("area_correction"))

  p <- std.info |>
    ggplot(aes(y = get("area")
               , x = get("mean_RT"))) +
    geom_vline(aes(xintercept = get("mean_RT"))
               , linetype = "dotted") +
    geom_step(direction = "vh"
              , linewidth = 1
              , color = "orange") +
    geom_point(color = "black"
              , fill = "red"
              , shape = 21
              , size = 4) +
    geom_step(aes(y = get("corrected_area"))
              , direction = "vh"
              , size = 1
              , color = "green") +
    geom_point(aes(y = get("corrected_area"))
               , color = "black"
                 , fill = "black"
                 , shape = 21
               , size = 4) +
    geom_richtext(aes(x = get("mean_RT")
                      , y = max(get("area")) + 1500000
                      , label = paste0("C", get("Chain.length")))
                  , size = 3.5
                  , angle = 270
                  , fontface = "bold") +
    ggside::geom_ysideboxplot(orientation = "x") +
    theme_classic() +
    labs(title = paste0("Abundance of standards (observed an corrected) vs"
                        , " mean retention time")
         , x = "mean RT (minutes)"
         , y = "Abundance (area)"
         , subtitle = "**Observed abundance:**
         <span style = 'color:red'>**red dots**</span> connected by an
         <span style = 'color:orange'>**orange line**</span>. <br>
         **Corrected abundance:**
         <span style = 'color:black'>**black dots**</span> connected by a
         <span style = 'color:green'>**green line**</span>") +
    theme(plot.subtitle = element_markdown())
  print(p)

  return(std.info)
}
