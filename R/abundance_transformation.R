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
#' @param internal_standard_peak A character specifying the Peak (e.g. "P2") that
#' contians the internal standard added to the samples, to be used
#' `transformation` is set to "scale".
#'
#' @param internal_standard_amount A numeric value, indicating the known amount
#' of internal standard that was added to the samples.
#'
#' @import dplyr
#' @import purrr
#'
abundance_transformation <- function(master.table
                                     , transformation = c("percentage"
                                                          , "proportion"
                                                          , "scale")
                                     , internal_standard_peak = NULL
                                     , internal_standard_amount = NULL) {
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
                 lm(y ~ sample
                    , data = _)
    })

    # Use predict() to transform each column based on its corresponding
    # calibration curve
    transformed_samples <- map2(.x = master_daten
                                , .y = calibration_curves
                                , ~ predict(.y
                                            , data.frame(sample = .x))) |>
      # Replace any negative value with 0
      map(~ ifelse(. < 0, 0, .))
      # map2(master_daten
      #      , calibration_curves
      #      , ~ {
      #        # Use the predict() function to transform the values of the
      #        # current sample based on the corresponding calibration curve
      #        predict(.y
      #                , as.data.frame(.x))
      #      })

    # Plot calibration curves and transformed samples
    calibration_plots <-
      map2(calibration_curves
           , names(calibration_curves)
           , ~ {
             # Create a scatter plot of the calibration data points (0, 0) and
             # (internal_standard_amount, corresponding sample value) using red
             # points. Add a blue line to the plot that represents the linear
             # regression line of the calibration curve. Set the title and
             # subtitle of the plot to indicate that it is a calibration curve
             # and to indicate the sample name
             ggplot() +
               geom_point(data = data.frame(sample = c(0
                                                       , .x$data$sample[2])
                                            ,
                                            y = c(0
                                                  , internal_standard_amount)
                                            )
                          , aes(sample, y)
                          , color = "red") +
               geom_smooth(aes(sample
                               , y)
                           , method = "lm"
                           , data = .x
                           , color = "blue") +
               labs(title = "Calibration Curve"
                    , subtitle = .y
                    , x = "Sample"
                    , y = "Amount")
           })

    transformed_plots <-
      map2(transformed_samples
           , names(transformed_samples)
           , ~ {
             # Create a scatter plot of the transformed data points
             # using black points. Set the title and subtitle of the plot to
             # indicate that it is a transformed column and to indicate the
             # sample name
             ggplot() +
               geom_point(data = data.frame(sample = master_daten[[.y]]
                                            , y = internal_standard_amount)
                          , aes(sample
                                , y)
                          , color = "black") +
               labs(title = "Transformed Sample"
                    , subtitle = .y
                    , x = "Sample"
                    , y = "Amount")
           })

    # FOR TRIALS
    # pdf(here::here("data-raw"
    #                , 'calibration_curves.pdf'))
    # calibration_plots
    # dev.off()
    #
    # pdf(here::here("data-raw"
    #                , 'transformed_samples.pdf'))
    # transformed_plots
    # dev.off()

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
