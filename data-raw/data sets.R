## code to prepare `data sets` dataset goes here

# load the fucntions of the package
load_all()

# Samples_data_list ----
# Import the CSV files with the samples'integration results
samples_path_data <- list.files(path = system.file("extdata"
                                                   , package = "analyzeGC")
                                #  Get all CSV files in the folder
                                , pattern = ".CSV|.csv"
                                , full.names = T) |>
  # Do not include standards
  str_subset('STD', negate = T)

samples_data_list <- import_gcms_data(samples_path_data)

# Create the samples_data _list data file for the package
use_data(samples_data_list, overwrite = TRUE)

# Samples_info ----
# Load the samples list
samples_info  <- here::here("data-raw", "samples-list.csv") |>
  readr::read_csv()

# Create the samples_info data file for the package
use_data(samples_info, overwrite = TRUE)

#  ----

