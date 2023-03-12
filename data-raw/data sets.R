## code to prepare `data sets` dataset goes here

load_all()
samples_path_data <- list.files(path = system.file("extdata"
                                                   , package = "analyzeGC")
                                #  Get all CSV files in the folder
                                , pattern = ".CSV|.csv"
                                , full.names = T) |>
  # Do not include standards
  str_subset('STD', negate = T)

samples_data_list <- analyzeGC::import_gcms_data(samples_path_data)

use_data(samples_data_list, overwrite = TRUE)

samples_info  <- here::here("data-raw", "samples-list.csv") |>
  readr::read_csv()

use_data(samples_info, overwrite = TRUE)

