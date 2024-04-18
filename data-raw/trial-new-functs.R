load_all()

# Import data ----
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

# nest samples with standards ----
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


# comps_id_std ----
comps_id_std <- here::here("data-raw", "std_compounds-id.csv") |>
  readr::read_csv() |>
  rename("Peak" = contains(".1"))

std_info <- shape_hcstd_info(comps_id.STD = comps_id_std
                             , aligned_std = aligned_standards
                             , short_std_pattern = "L"
                             , long_std_pattern = "H")

# Add std names to nested alignment ----
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



