#' @title Transform peak abundance values in a master table
#'
#' @description This function takes a master table and transforms the peak
#' abundance values within each sample, based on the specified transformation
#' option.
#'
#' @param master.table A master table, for which the abundance values should be
#' transformed.
#'
#' @param transformation TA character string specifying the type of
#' transformation to perform. Possible values: "percentage" (default),
#' "proportion", or "scale".
#'
#' @param internal_standard_peak A character string specifying the Peak
#' (e.g., "P2") containing the internal standard added to the samples,
#' to be used when the `transformation` argument is set to "scale".
#'
#' @param internal_standard_amount A numeric value indicating the known amount
#' of internal standard added to the samples.
#'
#'@param calibration_plot A logical value indicating whether to generate
#' calibration plots for each sample column when `transformation` is set to
#' "scale".
#'
#' @import dplyr
#' @import purrr
#' @import ggplot2
#'
#' @examples
#'
#' # Percentage transformation
#' master_table_transformed <-
#' abundance_transformation(master_table_reassembled)
#'
#' # Proportion transformation
#' master_table_transformed <-
#' abundance_transformation(master_table_reassembled
#'                          , transformation = "proportion")
#'
#' # Scale transformation
#' master_table_transformed <-
#' abundance_transformation(master_table_reassembled
#'                          , transformation = "scale"
#'                          , internal_standard_peak = "P4"
#'                          , internal_standard_amount = 250
#'                          , calibration_plot = F)
#'
abundance_transformation <- function(master.table
                                     , transformation = c("percentage"
                                                          , "proportion"
                                                          , "scale")
                                     , internal_standard_peak = NULL
                                     , internal_standard_amount = NULL
                                     , calibration_plot = T) {
  # Match the transformation argument to one of the allowed options
  transformation <- match.arg(transformation)

  # Select only the columns that contain the abundance values of the peaks for
  # each sample
  master_daten <- master.table |>
    select_if(is.double)

  # Replace any NA values with 0
  master_daten[is.na(master_daten)] <- 0

  # Select the columns containing the information of the compunds in the peaks
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
    if (is.null(internal_standard_peak)) {
      stop("internal_standard_peak must be provided for scale transformation")
    }

    # Check that a calibration row was provided
    if (is.null(internal_standard_amount)) {
      stop("internal_standard_amount must be provided for scale transformation")
    }

    if (!is.numeric(internal_standard_amount)) {
      stop("internal_standard_amount must be a numeric value")
    }

    # Obtain the internal_standard_index from the Peak column
    internal_standard_index <-
      which(compound_variables$Peak == internal_standard_peak)

    #  Create a calibration curve for each column using lm()
    calibration_curves <-
      lapply(master_daten
             , function(sample) {
               # Create a data frame with two rows, one for the origin (0, 0)
               # and one for the internal standard amount at the corresponding
               # index of the sample
               data.frame(sample = c(0, sample[internal_standard_index])
                          , y = c(0, internal_standard_amount)) |>
                 # Use the lm() function to fit a linear model of y as a
                 # function of sample
                 lm(y ~ 0 + sample
                    , data = _)
    })

    if (calibration_plot == T) {
      create_plot <- function(model, index) {
        # Get the slope and intercept from the model
        slope <- coef(model)

        # Get the two points used to generate the model
        data <- data.frame(
          sample = c(0, model$model$sample[2]),
          y = c(0, model$model$y[2])
        )

        # Create the ggplot with an abline and the two points, and include the
        # index in the title
        ggplot(data, aes(x = get("sample"), y = get("y"))) +
          geom_point(size = 3) +
          geom_abline(slope = slope
                      , intercept = 0
                      , color = "blue") +
          labs(x = "Sample Value"
               , y = "Internal Standard Amount"
               , title = paste("Calibration curve for sample", index)) +
          theme_classic()
      }

      # Use map2() to create a list of ggplots, one for each model, with the
      # corresponding index included in the title
      calibration_plots <- map2(calibration_curves
                                , names(calibration_curves)
                                , create_plot)
      print(calibration_plots)
    }

    # Use predict() to transform each column based on its corresponding
    # calibration curve
    transformed_samples <- map2(.x = master_daten
                                , .y = calibration_curves
                                , ~ predict(.y
                                            , data.frame(sample = .x))) |>
      # Replace any negative value with 0
      map(~ ifelse(. < 0, 0, .))

    # Combine transformed columns into a single master_daten frame
    master_daten <- as.data.frame(transformed_samples)

    # Set column names to match original master_daten frame
    names(master_daten) <- names(transformed_samples)
  }

  # reassemble the master table
  master.table <- compound_variables |>
    bind_cols(master_daten)

  return(master.table)
}
