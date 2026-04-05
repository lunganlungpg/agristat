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

# ============================================================================
# Step 21: Error-term helper functions (Phase 3)
# ============================================================================

#' Get the Error Structure of a Field Design
#'
#' Returns the `error_structure` list embedded in a `field_design` object,
#' providing degrees of freedom and descriptions for each error stratum.
#'
#' @param design A `field_design` object.
#'
#' @return A named list.  Each element describes one error stratum and
#'   contains at minimum `$df` (degrees of freedom) and optionally
#'   `$description`.
#'
#' @examples
#' d <- design_split_plot(
#'   whole_plot_factor = c("I+", "I-"),
#'   sub_plot_factor   = c("V1", "V2", "V3"),
#'   blocks = 4L, seed = 1L
#' )
#' get_error_structure(d)
#'
#' @seealso [compute_error_terms()], [extract_variance_components()]
#' @export
get_error_structure <- function(design) {
  if (!inherits(design, "field_design")) {
    rlang::abort("`design` must be a `field_design` object.")
  }
  design$error_structure
}

#' Extract Variance Components from a Mixed-Effects Model
#'
#' Extracts the random-effects variance components from a model fitted by
#' [lme4::lmer()].  Returns a named numeric vector of estimated variances
#' for each random-effects grouping factor plus the residual.
#'
#' @param model A `lmerMod` object (from [lme4::lmer()]) or an
#'   `agristat_aov` wrapping one.
#'
#' @return A named numeric vector.  Names correspond to grouping factors and
#'   `"Residual"`.  Returns `NA_real_` if the model has no random effects.
#'
#' @examples
#' \donttest{
#' d  <- design_split_plot(c("I+", "I-"), c("V1", "V2", "V3"),
#'                         blocks = 4L, seed = 1L)
#' df <- d$layout
#' df$yield <- rnorm(nrow(df), mean = 10, sd = 2)
#' m  <- analyze_design(d, response = "yield", data = df)
#' extract_variance_components(m)
#' }
#'
#' @seealso [compute_error_terms()], [get_error_structure()]
#' @export
extract_variance_components <- function(model) {
  # Strip agristat_aov wrapper so we can dispatch on the underlying class
  if (inherits(model, "agristat_aov")) {
    class(model) <- class(model)[class(model) != "agristat_aov"]
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    rlang::abort("Package 'lme4' is required for extract_variance_components().")
  }
  if (!inherits(model, "lmerMod")) {
    rlang::warn(
      paste0("extract_variance_components() is designed for lmerMod objects.",
             "  Returning NA.")
    )
    return(NA_real_)
  }
  vc     <- lme4::VarCorr(model)
  # Named vector: one entry per grouping factor + residual
  sigmas <- vapply(vc, function(x) attr(x, "stddev")^2,
                   FUN.VALUE = numeric(1L))
  names(sigmas) <- names(vc)
  sigmas[["Residual"]] <- attr(vc, "sc")^2
  sigmas
}

# ============================================================================
# Phase 4 (Step 31-33): Genetics helper utilities
# ============================================================================

#' Format an Effect Size for Display
#'
#' Formats a regression coefficient (effect size) with its standard error
#' as a combined string, suitable for GWAS result tables and reports.
#'
#' @param beta Numeric.  Estimated effect size (regression coefficient).
#' @param se Numeric.  Standard error of the estimate.
#' @param digits Integer.  Number of decimal places (default `4L`).
#'
#' @return A character string of the form `"beta (se)"`.
#'
#' @examples
#' format_effect_size(0.1234, 0.0456)
#'
#' @export
format_effect_size <- function(beta, se, digits = 4L) {
  if (anyNA(c(beta, se))) return("NA (NA)")
  sprintf(
    "%s (%s)",
    formatC(beta, format = "f", digits = digits),
    formatC(se,   format = "f", digits = digits)
  )
}

#' Bootstrap Confidence Interval for Heritability
#'
#' A thin wrapper that generates bootstrap confidence intervals for any
#' heritability ratio from raw variance component vectors.  Intended as a
#' reusable building block.
#'
#' @param v_g Numeric.  Genetic (numerator) variance.
#' @param v_e Numeric.  Residual (denominator) variance.
#' @param n_boot Integer.  Number of bootstrap replicates (default `1000L`).
#' @param ci_level Numeric.  Confidence level (default `0.95`).
#' @param seed Integer or `NULL`.  Random seed.
#'
#' @return A two-element numeric vector `c(lower, upper)`.
#'
#' @keywords internal
compute_heritability_ci <- function(v_g, v_e, n_boot = 1000L,
                                    ci_level = 0.95, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha <- (1 - ci_level) / 2
  boots <- replicate(n_boot, {
    bvg <- stats::rgamma(1, shape = max(v_g, 1e-9), rate = 1)
    bve <- stats::rgamma(1, shape = max(v_e, 1e-9), rate = 1)
    bvg / (bvg + bve)
  })
  stats::quantile(boots, c(alpha, 1 - alpha), names = FALSE)
}

#' Validate Marker Genotype Data
#'
#' Performs quality-control checks on a marker genotype matrix and returns
#' a summary QC report.  Checks include missing-data rate, minor allele
#' frequency, and monomorphic markers.
#'
#' @param marker_mat Numeric matrix (individuals x markers, 0/1/2 coding).
#' @param min_maf Numeric.  Minimum acceptable minor allele frequency
#'   (default `0.05`).
#' @param max_missing Numeric.  Maximum acceptable per-marker missing rate
#'   (default `0.20`).
#'
#' @return A list with components:
#'   \describe{
#'     \item{qc_table}{Data frame with per-marker statistics.}
#'     \item{n_fail_maf}{Number of markers below `min_maf`.}
#'     \item{n_fail_missing}{Number of markers above `max_missing`.}
#'     \item{n_monomorphic}{Number of monomorphic markers.}
#'     \item{pass_markers}{Character vector of markers passing all filters.}
#'   }
#'
#' @keywords internal
validate_marker_data <- function(marker_mat, min_maf = 0.05,
                                 max_missing = 0.20) {
  if (!is.matrix(marker_mat))
    rlang::abort("`marker_mat` must be a numeric matrix.")

  n_mrk      <- ncol(marker_mat)
  mrk_names  <- if (!is.null(colnames(marker_mat)))
                  colnames(marker_mat)
                else paste0("M", seq_len(n_mrk))

  miss_rate  <- apply(marker_mat, 2, function(x) mean(is.na(x)))
  allele_frq <- apply(marker_mat, 2, function(x) mean(x, na.rm = TRUE) / 2)
  maf        <- pmin(allele_frq, 1 - allele_frq)
  mono       <- apply(marker_mat, 2,
                      function(x) length(unique(x[!is.na(x)])) <= 1L)

  qc_table <- data.frame(
    marker       = mrk_names,
    missing_rate = miss_rate,
    maf          = maf,
    monomorphic  = mono,
    stringsAsFactors = FALSE
  )

  fail_maf     <- maf < min_maf
  fail_missing <- miss_rate > max_missing
  pass_markers <- mrk_names[!fail_maf & !fail_missing & !mono]

  list(
    qc_table        = qc_table,
    n_fail_maf      = sum(fail_maf, na.rm = TRUE),
    n_fail_missing  = sum(fail_missing, na.rm = TRUE),
    n_monomorphic   = sum(mono, na.rm = TRUE),
    pass_markers    = pass_markers
  )
}
