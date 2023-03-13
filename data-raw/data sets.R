## code to prepare `data sets` dataset goes here

library(here)
# load the functions of the package
load_all()

# Samples_data_list ----
# Import the CSV files with the samples'integration results
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

#  ----

samples_list_RT <- aligned_samples_data_list |>
  lapply(RT_df)

samples_list_area <- aligned_samples_data_list |>
  lapply(area_df)

#  ----
samples_area_norm_list <- aligned_samples_data_list |>
  lapply(area_norm)

#  ----

samples_area_norm_list$`Winter_In-hive workers_A. m. mellifera` |>
  diagnostic_heatmap(title = "trial"
                     , alignment.type = "automatic")






