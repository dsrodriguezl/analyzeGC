#' @title Import GC integration results from CSV  files
#'
#' @description A function to import the integration results data frames,
#' as obtained from MassHunter.
#'
#' @param path Path of the CSV file with the exported integration results of the
#'  samples.
#'
#' @param patterns_2_delete a character type vector, listing the text  patterns
#' to be removed from the CSV files' names, to extract the names of the samples
#' from them. This step is performed via str_remove_all.
#'
#' @param zip_export Logical value indicating whether to export a ZIP file
#' with all the separate CSV files of each sample (default: FALSE).
#'
#' @returns A list of tibble data frames with the integration results
#' (RT: retention time of each peak, and Area: area under the curve of each
#' peak).
#' Each tibble on the list corresponds to the CSV file of an individual
#' sample.
#'
#' @import stringr
#' @import readr
#' @import dplyr
#' @import tibble
#' @import tools
#'
#'
#' @export

import_mh_integration <- function(path
                             , patterns_2_delete = " "
                           , zip_export = F){

  # Temporary folder for CSVs
  tmpdir     <- file.path(tempdir(), "split_csv_tmp")
  zipfile     <- file.path(dirname(path),
                           paste0(tools::file_path_sans_ext(basename(path)),
                                  "_tables.zip"))

  # Aliases
  as_tibble <- tibble::as_tibble
  read_csv  <- function(txt) readr::read_csv(txt
                                             , show_col_types = FALSE
                                             , progress = FALSE)
  write_csv <- function(x, csv_path) readr::write_csv(x, csv_path)

  # Detect separator/blank rows (only commas or empty)
  is_sep <- function(x) {
    grepl("^\\s*(,\\s*)*$", x) | trimws(x) == ""
  }

  # Clean sample name: remove " + TIC Scan " and ".D..." and normalize
  clean_id <- function(firstline) {
    # remove "+ TIC Scan"
    id <- trimws(gsub("\\s*\\+\\s*TIC\\s*Scan\\s*", "", firstline))
    # remove ".D" and everything after
    id <- sub("\\.D.*", "", id)
    # drop any path
    id <- basename(id)
    id <- trimws(id)
    # Sanitize for filenames and list names
    id <- gsub("[^A-Za-z0-9_\\-]", "_", id)
    id
  }

  # Read stacked file
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  # Locate the start of each table using the marker "+ TIC Scan"
  starts <- grep("\\+\\s*TIC\\s*Scan", lines)

  if (length(starts) == 0) {
    stop(paste("No '+ TIC Scan' markers found."
               ,"Check the file content and marker pattern."))
  }

  # Determine the end line for each table (line before the next start; last
  # line for final table)
  ends <- c(starts[-1] - 1, length(lines))

  # Split, parse to tibbles, store in named list
  tables <- list()
  raw_names <- character(length(starts))

  for (i in seq_along(starts)) {
    # Extract chunk for this table
    chunk <- lines[starts[i]:ends[i]]

    # Trim leading/trailing separator/blank lines
    while (length(chunk) > 0 && is_sep(chunk[1])) chunk <- chunk[-1]
    while (length(chunk) > 0 && is_sep(chunk[length(chunk)])) chunk <-
        chunk[-length(chunk)]

    if (length(chunk) == 0) next  # skip empty chunks (defensive)

    # First line holds the sample marker/path
    firstline <- chunk[1]

    # Derive a friendly sample name
    id <- clean_id(firstline)
    if (identical(id, firstline) || !nzchar(id)) id <- sprintf("part_%02d", i)
    raw_names[i] <- id

    # Data lines = everything after the first line
    data_lines <- if (length(chunk) > 1) chunk[-1] else character(0)

    # Trim any leading/trailing separator/blank lines within data_lines
    while (length(data_lines) > 0 && is_sep(data_lines[1])) data_lines <-
      data_lines[-1]
    while (length(data_lines) > 0 &&
           is_sep(data_lines[length(data_lines)])) data_lines <-
      data_lines[-length(data_lines)]

    # Parse into a tibble
    if (length(data_lines) == 0) {
      tbl <- as_tibble(list())  # empty tibble
    } else {
      csv_text <- paste(data_lines, collapse = "\n")
      # Prefer readr::read_csv (tibble output); fallback to base if necessary
      tbl <- tryCatch(
        read_csv(csv_text),
        error = function(e) {
          con <- textConnection(csv_text)
          on.exit(close(con), add = TRUE)
          df <- utils::read.csv(con, stringsAsFactors = FALSE
                                , check.names = FALSE)
          as_tibble(df)
        }
      )
    }

    # Temporarily store; final naming (unique) applied after the loop
    tables[[i]] <- tbl
  }

  # Apply unique names to handle duplicates gracefully
  final_names <- raw_names |>
    str_remove_all(paste(patterns_2_delete
                         , ".CSV"
                         , ".csv"
                         , sep = "|"))
  names(tables) <- final_names



  # Create ZIP file and remove temp CSVs
  if (zip_export) {
    # Create temp folder
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
    outpaths <- character(0)
    for (nm in names(tables)) {
      path <- file.path(tmpdir, paste0(nm, ".csv"))
      write_csv(tables[[nm]], path)
      outpaths <- c(outpaths, path)
    }

    # Create ZIP and clean up
    if (file.exists(zipfile)) unlink(zipfile)
    utils::zip(zipfile, files = outpaths, flags = "-j")
    unlink(tmpdir, recursive = TRUE)

    message("ZIP file exported to: ", zipfile)
  }

  tables |>
    lapply(function(df)(df |>
                          select(all_of(c("RT","Area")))
    ))
}
