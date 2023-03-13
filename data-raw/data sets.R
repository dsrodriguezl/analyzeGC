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
samples_info  <- here("data-raw", "samples-list.csv") |>
  readr::read_csv()

# Create the samples_info data file for the package
use_data(samples_info, overwrite = TRUE)

#  ----

