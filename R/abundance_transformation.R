#' Transform a numeric master_daten frame
#'
#' This function takes a master table and transforms the peak abundance values
#' within each sample, based on the specified transformation option.
#'
#' @param master.table A master table, for which the abundance values should be
#' transformed.
#'
#' @param transformation The type of transformation to perform.
#' Can be one of "percentage", "proportion", or "scale".
#' Defaults to "percentage".
#'
#' @param internal_standard A character specifying the Peak (e.g. "P2") that
#' contians the internal standard added to the samples, to be used
#' `transformation` is set to "scale".
#'
#' @import dplyr
#'
abundance_transformation <- function(master.table
                                     , transformation = c("percentage"
                                                          , "proportion"
                                                          , "scale")
                                     , internal_standard = NULL) {
  # Match the transformation argument to one of the allowed options
  transformation <- match.arg(transformation)

  master_daten <- master.table |>
    select_if(is.double)

  master_daten[is.na(master_daten)] <- 0

  compound_variables <- master.table |>
    select(-all_of(colnames(master_daten)))

  # If the transformation is "percentage", divide each value in a column by the
  # sum of that column and multiply by 100
  if (transformation == "percentage") {
    col_sums <- colSums(master_daten)
    master_daten <- sweep(master_daten, 2, col_sums, "/") * 100
  }

  # If the transformation is "proportion", divide each value in a column by the
  # sum of that column
  if (transformation == "proportion") {
    col_sums <- colSums(master_daten)
    master_daten <- sweep(master_daten, 2, col_sums, "/")
  }

  # If the transformation is "scale", use a calibration curve to transform the
  # values column-wise
  if (transformation == "scale") {
    # Check that a calibration row was provided
    if (is.null(internal_standard)) {
      stop("internal_standard must be provided for scale transformation")
    }

    # Obtain the internal_standard_index from the Peak column
    internal_standard_index <-
      which(compound_variables$Peak == internal_standard)

    # Create a calibration curve for each column using lm()
    calibration_curves <- lapply(master_daten, function(column) {
      lm(column ~ c(0, column[internal_standard_index]))
    })

    # Use predict() to transform each column based on its corresponding
    # calibration curve
    transformed_columns <- Map(function(column, calibration_curve) {
      predict(calibration_curve, as.data.frame(column))
    }, master_daten, calibration_curves)

    # Combine transformed columns into a single master_daten frame
    master_daten <- as.data.frame(transformed_columns)

    # Set column names to match original master_daten frame
    names(master_daten) <- names(transformed_columns)
  }

  # reassemble the master table
  master.table <- compound_variables |>
    bind_cols(master_daten)

  return(master.table)
}
