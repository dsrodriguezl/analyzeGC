load_all()


# linear regression for projecting mean RTs of standards ----
comps_id.STD = comps_id_std
aligned_std = aligned_standards
short_long_splitted = TRUE
short_std_pattern = "L"
long_std_pattern = "H"
project_std = c(4, 5, 6, 7, 8, 9)

# If the aligned_std object is a GCalignR object, extract the alignment
if (length(aligned_std) > 2) {
  aligned_std <- aligned_std[["aligned"]]
}
std_area <- aligned_std[["Area"]]
std_RT <- aligned_std[["RT"]]

# Temporal data frame
std_df <- comps_id.STD |>
  select(contains("Peak")
         , contains("Compound")) |>
  bind_cols(std_RT)

std_df[std_df == 0] <- NA

if (short_long_splitted == TRUE) {
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

}

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
  select(all_of(std_RT |>
                  select(-contains("mean_RT")) |>
                  colnames()))

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


if(!is.null(project_std)) {

  project_std_df <- data.frame("Compound" = paste0("STD-C"
                                                   , project_std)
                               , "mean_RT" = as.numeric(NA)
                               , "Chain.length" = project_std |>
                                 as.integer()
                               , "Class" = "STD"
                               , "Mod.position" = "NA")

  m_dif  <- std.info |>
    pull("mean_RT") |>
    diff() |>
    median()

  project_std_df <- project_std |>
    lapply(function(x) {
      bigger <- std.info |>
        filter(Chain.length > x)

      smaller <- std.info |>
        filter(Chain.length < x)

      if (nrow(smaller) > 0) {
        bigger_exist <- nrow(bigger) > 0

        df <- smaller |>
          slice(nrow(smaller)) |>
          (function(y) {
            project_std_df |>
              filter(get("Chain.length") == x) |>
              mutate(mean_RT =
                       (y |>
                          pull("mean_RT") +
                          (abs(y |>
                              pull("Chain.length") - x) * m_dif)))
          })()
      }

      if (nrow(bigger) > 0) {
        df <- bigger |>
          slice(1) |>
          (function(y) {
            project_std_df |>
              filter(get("Chain.length") == x) |>
              mutate(mean_RT =
                       (y |>
                          pull("mean_RT") -
                          ((y |>
                              pull("Chain.length") - x) * m_dif)))

          })()
      }
      df
    }) |>
    reduce(rbind)


  std.info <- std.info |>
    rows_insert(project_std_df) |>
    arrange(get("mean_RT"))
}


# RI-wise alignment trial ----
## Import data ----
## Samples
# Import the CSV files with the samples' integration results
samples_path_data <- list.files(path = system.file("extdata/gcms_integration"
                                                   , package = "analyzeGC")
                                #  Get all CSV files in the folder
                                , pattern = ".CSV|.csv"
                                , full.names = T) |>
  # Do not include standards
  str_subset('STD', negate = T)


samples_data_list <- import_mh_data(samples_path_data
                                    , patterns_2_delete = "DR_")

## Standards
standards_path_data <- list.files(path = system.file("extdata/gcms_integration"
                                                     , package = "analyzeGC")
                                  #  Get all CSV files in the folder
                                  , pattern = ".CSV|.csv"
                                  , full.names = T) |>
  # Only include standards
  str_subset('STD')


standards_data_list <- import_mh_data(standards_path_data
                                      , patterns_2_delete = "STD")

## nest samples with standards ----
nested_list <- samples_data_list |>
  names() |>
  lapply(function (muestra) {
    list(samples_data_list |>
      keep_at(muestra)
      , standards_data_list) |>
      list_flatten()
  }) |>
  set_names(names(samples_data_list))

aligned_nested <- nested_list |>
  lapply(align_chromatograms2
         , blanks = NULL
         , linear_shift_criteria = 0.02
         , partial_alignment_threshold = 0.05
         , row_merging_threshold = 0.15)

aligned_standards <- align_chromatograms2(standards_data_list
                                          , blanks = NULL
                                          , linear_shift_criteria = 0.02
                                          , partial_alignment_threshold = 0.05
                                          , row_merging_threshold = 0.15)


## comps_id_std ----
comps_id_std <- here::here("data-raw", "std_compounds-id.csv") |>
  readr::read_csv() |>
  rename("Peak" = contains(".1"))

std_info <- shape_hcstd_info(comps_id.STD = comps_id_std
                             , aligned_std = aligned_standards
                             , short_std_pattern = "L"
                             , long_std_pattern = "H")

## Add std names to nested alignment ----
aligned_nested2 <- aligned_nested  |>
  lapply(function(x) {
    x |>
      correct_alignment(movements_list = list()) |>
      recalculate_meanRT()
  })

aligned_nested3 <-
  add_comps_info(samples.list = aligned_nested2
                 , comps.info.list = aligned_nested2 |>
                   lapply(function(x) {
                     x |>
                       (function(l) {
                         df_RT <- l |>
                           extract_RT() |>
                           mutate(Compound = c(NA) |> as.character()) |>
                           rows_update(comps_id_std |>
                                         select(contains("L"), Compound) |>
                                         filter(!is.na(Compound)) |>
                                         mutate_if(is.double, function(x){
                                           ifelse(x == 0, NA, x)
                                           }) |>
                                         drop_na()) |>
                           rows_update(comps_id_std |>
                                         select(contains("H"), Compound) |>
                                         filter(!is.na(Compound)) |>
                                         mutate_if(is.double, function(x){
                                           ifelse(x == 0, NA, x)
                                           }) |>
                                         drop_na()) |>
                           select(Peak, Compound, everything())
                         df_Area <- l |>
                           pluck("Area") |>
                           mutate(Compound = c(NA) |>
                                    as.character()) |>
                           rows_update(df_RT |>
                                         select(Peak:Compound)) |>
                           select(Peak, Compound, everything())

                         l[["RT"]] <- df_RT
                         l[["Area"]] <- df_Area

                         l
                        })() |>
                       extract_RT() |>
                       select(Peak:Compound)
                     }) |>
                   lapply(get_chc_info))

aligned_nested_ri <- aligned_nested3 |>
  lapply(kovats_retention_index, std.info = std_info)

align_chromatograms_ri <- function(data2align
                                 , blanks = NULL
                                 , linear_shift_criteria
                                 , partial_alignment_threshold
                                 , row_merging_threshold){
  # if (length(data2align) == 1) {
  #   warning(paste("data2align contains only one sample!"
  #                 , "The data will be formated in a list with the RT"
  #                 , "and Area values separated, but no alignment procedure"
  #                 , "will be performed"))
  #   nombre <- names(data2align)
  #   row_names <- paste0("P", 1:nrow(data2align[[1]]))
  #
  #   RT <- data2align[[1]] |>
  #     mutate("mean_RT" = get("RT")) |>
  #     select(contains("mean_RT"), contains("RT")) |>
  #     as.data.frame()
  #   colnames(RT) <- c("mean_RT", nombre)
  #   row.names(RT) <- row_names
  #
  #   Area <- data2align[[1]] |>
  #     mutate("mean_RT" = get("RT")) |>
  #     select(contains("mean_RT"), contains("Area")) |>
  #     as.data.frame()
  #   colnames(Area) <- c("mean_RT", nombre)
  #   row.names(Area) <- row_names
  #
  #   df <- list("RT" = RT, "Area" = Area)
  # }

  if (length(data2align) > 1) {
    withr::local_seed(12345)
    df <- align_chromatograms(data = data2align
                              , rt_col_name = "RI"
                              , max_linear_shift = linear_shift_criteria
                              , max_diff_peak2mean = partial_alignment_threshold
                              , min_diff_peak2peak = row_merging_threshold
                              , blanks = blanks
    )

    row.names(df[["aligned"]][["RI"]]) <-
      paste0("P", 1:nrow(df[["aligned"]][["RI"]]))
    row.names(df[["aligned"]][["Area"]]) <-
      paste0("P", 1:nrow(df[["aligned"]][["Area"]]))
  }

  df
}


ri_aligned <- aligned_nested_ri |>
  imap(function(s_data, s_id) {
    s_data |>
      pluck("Area") |>
      select(RI, all_of(s_id)) |>
      mutate(Area = get(s_id)
             , .keep = "unused") |>
      as.data.frame()
  }) |>
  align_chromatograms_ri(linear_shift_criteria = 1
                         , partial_alignment_threshold = 3
                         , row_merging_threshold = 2)

ri_aligned |>
  diagnostic_heatmap(title = "linear_shift_criteria = 1
                         , partial_alignment_threshold = 3
                         , row_merging_threshold = 2"
                     , alignment.type = "automatic")



