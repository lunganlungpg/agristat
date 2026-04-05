# Internal utility functions for agristat
#
# Shared helpers used across all six modules.

# ============================================================================
# Formatting helpers
# ============================================================================

#' Format a P-Value for Display
#'
#' Formats a numeric p-value as a character string suitable for publication
#' tables and console output.  Very small p-values are shown as `"< 0.001"`.
#'
#' @param p Numeric p-value (length 1).  `NA` returns `"NA"`.
#' @param digits Integer.  Number of decimal places when `p >= 0.001`
#'   (default `4`).
#'
#' @return A character string.
#'
#' @examples
#' format_p_value(0.0432)
#' format_p_value(0.00005)
#' format_p_value(NA)
#'
#' @export
format_p_value <- function(p, digits = 4L) {
  if (length(p) == 0L || is.na(p)) return("NA")
  if (!is.numeric(p)) return(as.character(p))
  if (p < 0.001) return("< 0.001")
  formatC(p, format = "f", digits = digits)
}

#' Format a Test Statistic for Display
#'
#' Formats a numeric statistic as a fixed-decimal character string.
#'
#' @param stat Numeric statistic (length 1).  `NA` returns `"NA"`.
#' @param digits Integer.  Number of decimal places (default `4`).
#'
#' @return A character string.
#'
#' @examples
#' format_statistic(3.14159)
#' format_statistic(NA)
#'
#' @export
format_statistic <- function(stat, digits = 4L) {
  if (length(stat) == 0L || is.na(stat)) return("NA")
  if (!is.numeric(stat)) return(as.character(stat))
  formatC(stat, format = "f", digits = digits)
}

#' Interpret a P-Value as Pass or Fail
#'
#' Compares a p-value against a significance level and returns a character
#' label `"PASS"` or `"FAIL"`.
#'
#' @param p Numeric p-value.
#' @param alpha Numeric significance level.
#' @param fail_below Logical.  If `TRUE` (default), the test FAILS when
#'   `p < alpha`; set to `FALSE` to invert the logic (used when a high
#'   p-value indicates failure, e.g., a normality test).
#'
#' @return Character string `"PASS"`, `"FAIL"`, or `"UNKNOWN"` if `p`
#'   is `NA`.
#'
#' @keywords internal
interpretation_helper <- function(p, alpha, fail_below = TRUE) {
  if (is.na(p)) return("UNKNOWN")
  if (fail_below) {
    if (p < alpha) "FAIL" else "PASS"
  } else {
    if (p >= alpha) "PASS" else "FAIL"
  }
}

# ============================================================================
# Validation helpers
# ============================================================================

#' Validate that an Object is a Non-Empty Data Frame
#'
#' @param df Object to check.
#' @param arg_name Character.  Name of the argument for error messages.
#'
#' @return `df` invisibly (throws an error if validation fails).
#'
#' @keywords internal
.validate_data_frame <- function(df, arg_name = "data") {
  if (!is.data.frame(df)) {
    rlang::abort(sprintf("`%s` must be a data frame.", arg_name))
  }
  if (nrow(df) == 0L) {
    rlang::abort(sprintf("`%s` must have at least one row.", arg_name))
  }
  invisible(df)
}

#' Check that Required Columns Exist in a Data Frame
#'
#' @param df Data frame.
#' @param cols Character vector of required column names.
#' @param arg_name Character.  Name of the data-frame argument.
#'
#' @return `df` invisibly (throws an error if a column is missing).
#'
#' @keywords internal
.check_column_exists <- function(df, cols, arg_name = "data") {
  missing_cols <- cols[!cols %in% names(df)]
  if (length(missing_cols) > 0L) {
    rlang::abort(sprintf(
      "`%s` is missing required column(s): %s.",
      arg_name,
      paste(sprintf("'%s'", missing_cols), collapse = ", ")
    ))
  }
  invisible(df)
}

#' Convert a P-Value to Significance Stars
#'
#' @param p Numeric p-value.
#'
#' @return Character string: `"***"`, `"**"`, `"*"`, `"."`, or `" "`.
#'
#' @keywords internal
.stars_pvalue <- function(p) {
  if (is.na(p)) return(" ")
  if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else if (p < 0.1)   "."
  else                " "
}

#' Check Whether a Numeric Value is a Whole Number
#'
#' @param x Numeric value.
#' @param tol Tolerance (default `.Machine$double.eps^0.5`).
#'
#' @return Logical.
#'
#' @keywords internal
.is_whole_number <- function(x, tol = .Machine$double.eps^0.5) {
  is.numeric(x) && all(abs(x - round(x)) < tol)
}
