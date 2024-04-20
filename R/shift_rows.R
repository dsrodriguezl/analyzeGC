#' @title Shift rows in a data frame
#'
#' @description Shifts the values in specified columns of a data frame up or down for specified rows.
#' Rows can be specified using either row names or row numbers.
#'
#' @param df A data frame to modify.
#' @param cols_to_shift A tidy selection of columns to shift.
#' @param rows_to_shift A character vector of row names or a numeric vector of row numbers for which the values will be moved.
#' @param direction The direction to move the values; either "up" or "down".
#'
#' @return The modified data frame.
#'
#' @import dplyr
#' @importFrom rlang !!
#' @importFrom rlang enquo
#'
#' @examples
#' # Create a sample data frame with row names
#' df <- data.frame(col1 = c(1, 2, 3), col2 = c("a", "b", "c"), col3 = c(4, 5, 6))
#' rownames(df) <- c("row1", "row2", "row3")
#'
#' # Use the function to shift the values in col1 and col3 one row down for row2 (using row number)
#' df <- shift_rows(df, c(col1, col3), 2, direction = "down")
#'
#' @export
shift_rows <- function(df, cols_to_shift, rows_to_shift, direction) {
  # Check if df is a data frame
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }

  # Check if cols_to_shift is valid
  tryCatch({
    # Convert tidy selection to column names
    cols_to_shift <- df |>
      select({{cols_to_shift}}) |>
      colnames()
  }, error = function(e) {
    stop("cols_to_shift is not a valid tidy selection")
  })

  # Check if rows_to_shift is valid
  if (is.character(rows_to_shift)) {
    # Convert row names to row numbers
    rows_to_shift <- match(rows_to_shift, rownames(df))
    if (any(is.na(rows_to_shift))) {
      stop("One or more specified row names do not exist in df")
    }
  } else if (is.numeric(rows_to_shift)) {
    if (any(rows_to_shift < 1 | rows_to_shift > nrow(df))) {
      stop("One or more specified row numbers are outside the range of valid row numbers for df")
    }
  } else {
    stop("rows_to_shift must be a character vector (row names) or a numeric vector (row numbers)")
  }

  # Check if direction is valid
  if (!direction %in% c("up", "down")) {
    stop("direction must be either 'up' or 'down'")
  }

  # Check for edge cases
  if (direction == "up" && any(rows_to_shift == 1)) {
    stop("Cannot shift the first row up")
  }
  if (direction == "down" && any(rows_to_shift == nrow(df))) {
    stop("Cannot shift the last row down")
  }

  # Shift the values in the specified columns and rows
  if (direction == "up") {
    df[rows_to_shift - 1, cols_to_shift] <- df[rows_to_shift, cols_to_shift]
    df[rows_to_shift, cols_to_shift] <- 0

  } else if (direction == "down") {
    df[rows_to_shift + 1, cols_to_shift] <- df[rows_to_shift, cols_to_shift]
    df[rows_to_shift, cols_to_shift] <- 0
  }

  # Return the modified data frame
  df
}
