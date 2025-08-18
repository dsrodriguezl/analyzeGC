## code to prepare data sets for the package
load_all()

# Samples_data_list ----
data_list <-
  import_mh_integration(path =
                   system.file("extdata/gcms_integration/MS_integration.csv"
                               , package = "analyzeGC")
                 , patterns_2_delete = "DR_")

samples_data_list <- data_list |>
  # Do not include standards
  discard_at(str_detect(names(data_list), 'STD'))

# Create the samples_data _list data file for the package
use_data(samples_data_list, overwrite = TRUE)

# grouping_info ----
# Load the samples list
grouping_info  <- here::here("data-raw", "samples-list.csv") |>
  readr::read_csv() |>
  mutate("Individual" = as.character(get("Individual")))

grouping_info <- grouping_info |>
  mutate_at(vars(!matches("Individual"))
            , as.factor)

# Create the samples_info data file for the package
use_data(grouping_info, overwrite = TRUE)

# standards_data_list ----
# Import the CSV files with the samples' integration results
standards_data_list <-data_list |>
  # Do not include standards
  keep_at(str_detect(names(data_list), 'STD'))

# Create the standards_data_list data file for the package
use_data(standards_data_list, overwrite = TRUE)

# aligned_samples_data_list ----
grouping_info <- grouping_info |>
  unite(group_label
        , where(is.factor)
        , sep = "_"
        , remove = F)

samples_data_list <- mg_list(sample.info = grouping_info
                             , group.label = "group_label"
                             , samples.data.list = samples_data_list)

aligned_samples_data_list <-
  samples_data_list$`Winter_In-hive workers_A. m. mellifera` |>
  align_chromatograms2(blanks = NULL
                       , linear_shift_criteria = 0.02
                       , partial_alignment_threshold = 0.05
                       , row_merging_threshold = 0.15)

aligned_samples_data_list <- samples_data_list |>
  lapply(align_chromatograms2
         , blanks = NULL
         , linear_shift_criteria = 0.02
         , partial_alignment_threshold = 0.05
         , row_merging_threshold = 0.15)

# Create the aligned_samples_data_list data file for the package
use_data(aligned_samples_data_list, overwrite = TRUE)

# Export alignments
# for (df in names(aligned_samples_data_list)) {
#   write.csv(aligned_samples_data_list[[df]][["aligned"]][["RT"]]
#             , here::here("data-raw"
#                    , paste0("aligned-RT_", df, ".csv")))
#   rm(df)
# }

# aligned_standards ----
aligned_standards <- align_chromatograms2(standards_data_list
                                          , blanks = NULL
                                          , linear_shift_criteria = 0.02
                                          , partial_alignment_threshold = 0.05
                                          , row_merging_threshold = 0.15)

aligned_standards |>
  diagnostic_heatmap(title = "Alignment of standards")

# Create the aligned_standards data file for the package
use_data(aligned_standards, overwrite = TRUE)

# Export CSV to make the comps_id data sets
## It is commented to avoid overwriting the file after adding to it the
## compounds' ids
# write.csv(aligned_standards$aligned$RT
#           , here::here("data-raw", "std_compounds-id.csv"))

# Trying functions to diagnose alignment ----
# samples_area_norm_list <- aligned_samples_data_list |>
#   lapply(area_norm)

aligned_samples_data_list$`Winter_In-hive workers_A. m. mellifera` |>
  diagnostic_heatmap(title = "Alignment of IW CHCs")

pdf(here::here("data-raw"
         , "uncorrected-alignment-plots.pdf")
    , width = 30
    , height = 15)
for (df in names(aligned_samples_data_list)) {
  diagnostic_heatmap(aligned_samples_data_list[[df]]
                     , title = paste0("Alignment of "
                                      , df))
}

diagnostic_heatmap(aligned_standards
                   , title = paste0("Alignment of standards"))
dev.off()

# corrected_samples_list ----
 ## Trying code for adding empty peaks
empty_peaks <- list("335" = tribble(~position.reference, ~direction,
                                        "P100", "after")
                    , "339" = tribble(~position.reference, ~direction,
                                          "P100", "before"))

# Generate the data set
peaks_movements <- list("350" = data.frame(peaks_origin = c(paste0("P"
                                                                 , c(108
                                                                     , 126)))
                                           , peaks_target = c(paste0("P"
                                                                     , c(107
                                                                         , 125))))
                        , "351" = data.frame(peaks_origin = c(paste0("P"
                                                                   , c(107
                                                                       , 108
                                                                       , 119
                                                                       , 120
                                                                       , 126)))
                                             , peaks_target = c(paste0("P"
                                                                       , c(106
                                                                           , 107
                                                                           , 118
                                                                           , 119
                                                                           , 125))))
                        , "352" = data.frame(peaks_origin = c(paste0("P"
                                                                   , c(107
                                                                       , 108
                                                                       , 126)))
                                             , peaks_target = c(paste0("P"
                                                                       , c(106
                                                                           , 107
                                                                           , 125))))
                        , "333" = data.frame(peaks_origin = c(paste0("P"
                                                                   , c(26)))
                                             , peaks_target = c(paste0("P"
                                                                       , c(27))))
                        , "345" = data.frame(peaks_origin = c(paste0("P"
                                                                   , c(106
                                                                       , 107)))
                                             , peaks_target = c(paste0("P"
                                                                       , c(105
                                                                           , 106))))
                        )

corrected_samples_list <- lapply(aligned_samples_data_list
                                 , correct_alignment
                                 , new_peaks = empty_peaks
                                 , peak_movements = peaks_movements)

use_data(corrected_samples_list, overwrite = TRUE)

# Trying code to diagnose alignment correction
pdf(here::here("data-raw"
               , "corrected-alignment-plots.pdf")
    , width = 30
    , height = 15)
for (df in names(corrected_samples_list)) {
  diagnostic_heatmap(corrected_samples_list[[df]]
                     , title = paste0("corrected alignment of "
                                      , df)
                     # , alignment.type = "corrected"
                     )
}
dev.off()

# comps_id_std ----
comps_id_std <- here::here("data-raw", "std_compounds-id.csv") |>
  readr::read_csv() |>
  rename("Peak" = contains(".1"))

use_data(comps_id_std, overwrite = T)

# std_info ----
std_info <- shape_hcstd_info(comps_id.STD = comps_id_std
                             , aligned_std = aligned_standards
                             , short_std_pattern = "L"
                             , long_std_pattern = "H")

use_data(std_info, overwrite = T)

# comps_id_samples ----
comps_id_paths <- list.files(here::here("data-raw")
           , pattern = "compounds-id"
           , full.names = T) |>
  stringr::str_subset("std", negate = T)

comps_id_list <- comps_id_paths |>
  lapply(readr::read_csv)

names(comps_id_list) <- comps_id_paths |>
  stringr::str_split("/"
                     , simplify = T) |>
  stringr::str_subset(".CSV|.csv") |>
  stringr::str_split("-i"
                     , simplify = T) |>
  stringr::str_subset("Winter") |>
  stringr::str_split("_c"
                     , simplify = T) |>
  stringr::str_subset("Winter")

use_data(comps_id_list, overwrite = T)

# comps_info_list ----
comps_info_list <- comps_id_list |>
  lapply(get_hc_info)

use_data(comps_info_list, overwrite = T)

# adjusted_samples_list ----
pdf(here::here("data-raw"
         , "samples_correction-plots.pdf")
    , width = 12
    , height = 6)
adjusted_samples_list <- corrected_samples_list |>
  lapply(adjust_abundance, std.info = std_info)

dev.off()

use_data(adjusted_samples_list, overwrite = T)

# unfiltered_samples_list ----
unfiltered_samples_list <- add_comps_info(samples.list = adjusted_samples_list
                                          , comps.info.list = comps_info_list)

use_data(unfiltered_samples_list, overwrite = T)

# group_tables_list ----
group_tables_list <- unfiltered_samples_list |>
  lapply(shape_aligned_table
         , trace_comps_threshold = 0.01
         , drop_unidentified_comps = T
         , std.info = std_info)

use_data(group_tables_list, overwrite = T)

# master_table ----
master_table <- build_master_table(group_tables_list)

write.csv(master_table[["Area"]]
          , here::here("data-raw", "master_table.csv")
          , row.names = F)

use_data(master_table, overwrite = T)

# trying retrieve_group_tables ----
grouping_info <- grouping_info |>
  unite(group_label
        , where(is.factor)
        , sep = "_"
        , remove = FALSE)

retrieve_group_tables(group.label = "group_label"
                      , master.table = master_table
                      , grouping.info = grouping_info)

# duplicated_compounds_presence ----
pdf(here::here("data-raw"
         , 'density-distribution_duplicated-compounds.pdf')
    , width = 10, height = 8)
duplicated_compounds_presence <-
  retrieve_group_tables(group.label = "group_label"
                        , master.table = master_table
                        , grouping.info = grouping_info) |>
  assess_duplicated_compounds()
dev.off()

use_data(duplicated_compounds_presence, overwrite = T)

# master_table2 ----
fusion_list <- list(c(paste0("P"
                             , c(1, 2)))
                    , c(paste0("P"
                               , c(7, 8)))
                    , c(paste0("P"
                               , c(9, 10)))
                    , c(paste0("P"
                               , c(13, 14)))
                    , c(paste0("P"
                               , c(15, 16)))
                    , c(paste0("P"
                               , c(19, 20)))
                    , c(paste0("P"
                               , c(24, 25)))
                    , c(paste0("P"
                               , c(27, 28)))
                    , c(paste0("P"
                               , c(29, 30)))
                    , c(paste0("P"
                               , c(32, 33)))
                    , c(paste0("P"
                               , c(38, 39)))
                    , c(paste0("P"
                               , c(40, 41)))
                    , c(paste0("P"
                               , c(42, 43)))
                    , c(paste0("P"
                               , c(44, 45)))
                    , c(paste0("P"
                               , c(46, 47)))
                    , c(paste0("P"
                               , c(49, 50)))
                    , c(paste0("P"
                               , c(52, 53)))
                    , c(paste0("P"
                               , c(58, 59)))
                    , c(paste0("P"
                               , c(60, 61)))
                    , c(paste0("P"
                               , c(62, 63)))
                    , c(paste0("P"
                               , c(66, 67)))
                    , c(paste0("P"
                               , c(69, 70)))
                    , c(paste0("P"
                               , c(72, 73)))
                    , c(paste0("P"
                               , c(74, 75)))
                    , c(paste0("P"
                               , c(76, 77, 78)))
                    , c(paste0("P"
                               , c(80, 81)))
                    , c(paste0("P"
                               , c(82, 83)))
                    , c(paste0("P"
                               , c(84, 85)))
                    , c(paste0("P"
                               , c(86, 87)))
                    , c(paste0("P"
                               , c(88, 89))))
fusion_list

master_table2 <- fuse_all_peaks(master.table = master_table
                                , fusion.list = fusion_list)

use_data(master_table2, overwrite = T)

# group_tables_list2 ----
grouping_info <- grouping_info |>
  unite(group_label
        , where(is.factor)
        , sep = "_"
        , remove = FALSE)

group_tables_list2 <- retrieve_group_tables(group.label = "group_label"
                                            , master.table = master_table2
                                            , grouping.info = grouping_info) |>
  lapply(group_frequency_filter)

use_data(group_tables_list2, overwrite = T)

# master_table_reassembled ----
master_table_reassembled <- build_master_table(group_tables_list2)

use_data(master_table_reassembled, overwrite = T)

# master_table_transformed ----

#  Trying the scale transformation
# pdf(here::here("data-raw", "calibration_curves.pdf"))
# abundance_transformation(master.table = master_table_reassembled
#                          , transformation = "scale"
#                          , internal_standard_peak = "P4"
#                          , internal_standard_amount = 250
#                          , calibration_plot = T)
# dev.off()

master_table_final <- abundance_transformation(master_table_reassembled)

use_data(master_table_final, overwrite = T)

#  ----
