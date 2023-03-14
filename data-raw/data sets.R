## code to prepare `data sets` dataset goes here

library(here)
# load the functions of the package
load_all()

# Samples_data_list ----
# Import the CSV files with the samples' integration results
samples_path_data <- list.files(path = system.file("extdata/gcms_integration"
                                                   , package = "analyzeGC")
                                #  Get all CSV files in the folder
                                , pattern = ".CSV|.csv"
                                , full.names = T) |>
  # Do not include standards
  str_subset('STD', negate = T)

samples_data_list <- import_gcms_data(samples_path_data
                                      , patterns_2_delete = "DR_")

# Create the samples_data _list data file for the package
use_data(samples_data_list, overwrite = TRUE)

# Samples_info ----
# Load the samples list
grouping_info  <- here("data-raw", "samples-list.csv") |>
  readr::read_csv()

grouping_info <- grouping_info |>
  mutate_at(vars(!matches("Individual"))
            , as.factor)

# Create the samples_info data file for the package
use_data(grouping_info, overwrite = TRUE)

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

# samples_list_RT and samples_list_area ----

samples_list_RT <- aligned_samples_data_list |>
  lapply(RT_df)
use_data(samples_list_RT, overwrite = TRUE)

samples_list_area <- aligned_samples_data_list |>
  lapply(area_df)
use_data(samples_list_area, overwrite = TRUE)

  ## Trying functions to diagnose alignment
samples_area_norm_list <- aligned_samples_data_list |>
  lapply(area_norm)

samples_area_norm_list$`Winter_In-hive workers_A. m. mellifera` |>
  diagnostic_heatmap(title = "Alignment of IW CHCs"
                     , alignment.type = "automatic")

for (df in names(samples_area_norm_list)) {
  dev.new()
  diagnostic_heatmap(samples_area_norm_list[[df]]
                     , title = paste0("Alignment of "
                                      , df)
                     , alignment.type = "automatic")
}

# corrected_samples_list_area and  corrected_samples_list_RT ----
 ## Trying code for adding empty peaks
empty_peaks <- list("335" = tribble(~position.reference, ~direction,
                                        "P100", "after")
                    , "339" = tribble(~position.reference, ~direction,
                                          "P100", "before"))

area_IW <- samples_list_area$`Winter_In-hive workers_A. m. mellifera`
area_OW <- samples_list_area$`Winter_Out-hive workers_A. m. mellifera`

area_IW <- area_IW |>
  add_empty_peaks(empty.peaks = empty_peaks)

area_OW <- area_OW |>
  add_empty_peaks(empty.peaks = empty_peaks)

samples_list_area |>
  lapply(add_empty_peaks
         , empty.peaks = empty_peaks)

# Generate the data set
peaks_movements <- list("350" = data.frame(peaks_list = c(paste0("P"
                                                                 , c(106, 107
                                                                     , 124)))
                                           , movement_dirs = c('up', 'up'
                                                               , 'up'))
                        , "351" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(106, 107
                                                                       , 124)))
                                             , movement_dirs = c('up','up'
                                                                 ,'up'))
                        , "352" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(106, 107
                                                                       , 124
                                                                       , 148)))
                                             , movement_dirs = c('up','up'
                                                                 ,'up'
                                                                 ,'up'))
                        , "354" = data.frame(peaks_list =  'P144'
                                             , movement_dirs = 'up')
                        , "328" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(26, 35, 85
                                                                       , 124, 128)))
                                             , movement_dirs = c('up', 'up', 'up'
                                                                 , 'up', 'up'))
                        , "331" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(85, 128)))
                                             , movement_dirs = c('up', 'up'))
                        , "333" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(52, 128)))
                                             , movement_dirs = c('up', 'up'))
                        , "337" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(26, 52
                                                                       , 128)))
                                             , movement_dirs = c('up', 'up'
                                                                 , 'up'))
                        , "338" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(52)))
                                             , movement_dirs = c('up'))
                        , "341" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(52, 124
                                                                       , 128)))
                                             , movement_dirs = c('up','up'
                                                                 ,'up'))
                        , "345" = data.frame(peaks_list = c(paste0("P"
                                                                   , c(106, 107)))
                                             , movement_dirs = c('up', 'up'))
                        )

corrected_samples_list_area <- lapply(samples_list_area
                                      , correct_alignment
                                      , movements_list = peaks_movements)
use_data(corrected_samples_list_area, overwrite = TRUE)

corrected_samples_list_RT <- lapply(samples_list_RT
                                    , correct_alignment
                                    , movements_list = peaks_movements)
use_data(corrected_samples_list_RT, overwrite = TRUE)

# Trying code to diagnose alignment correction
for (df in names(corrected_samples_list_area)) {
  dev.new()
  diagnostic_heatmap(corrected_samples_list_area[[df]]
                     , title = paste0("corrected alignment of "
                                      , df)
                     , alignment.type = "corrected")
}

# corrected_samples_list_area2 and corrected_samples_list_RT2 ----
corrected_samples_list_area2 <-
  recalculate_meanRT(RT.list = corrected_samples_list_RT
                     , area.list = corrected_samples_list_area
                     , output.list = "Area")
use_data(corrected_samples_list_area2, overwrite = TRUE)

corrected_samples_list_RT2 <-
  recalculate_meanRT(RT.list = corrected_samples_list_RT
                     , area.list = corrected_samples_list_area
                     , output.list = "RT")
use_data(corrected_samples_list_RT2, overwrite = TRUE)









