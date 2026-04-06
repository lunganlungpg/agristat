# Module 1 — Data Ingestion & Pre-Processing
#
# Implements hypothesis tracking, assumption testing, data transformation,
# Mendelian segregation testing, and field data import for the agristat
# package.

# ============================================================================
# Step 8: Hypothesis Tracking Functions
# ============================================================================

#' Declare a Statistical Hypothesis
#'
#' Creates a structured S3 object capturing a null hypothesis (\eqn{H_0}) and
#' an alternative hypothesis (\eqn{H_1}) together with the significance level
#' (\eqn{\alpha}).  The returned object records the declaration date/time and
#' can be printed or summarised for inclusion in reports.
#'
#' @param h0 A non-empty character string stating the null hypothesis
#'   (\eqn{H_0}).
#' @param h1 A non-empty character string stating the alternative hypothesis
#'   (\eqn{H_1}).
#' @param alpha Numeric significance level strictly between 0 and 1
#'   (default `0.05`).
#'
#' @return An S3 object of class `"hypothesis"` with components:
#'   \describe{
#'     \item{`h0_text`}{Character. The null hypothesis statement.}
#'     \item{`h1_text`}{Character. The alternative hypothesis statement.}
#'     \item{`alpha`}{Numeric. Significance level.}
#'     \item{`declared_date`}{`POSIXct`. Date and time the object was created.}
#'   }
#'
#' @examples
#' h <- declare_hypothesis(
#'   h0    = "No significant difference in yield between varieties",
#'   h1    = "At least one variety differs in yield",
#'   alpha = 0.05
#' )
#' print(h)
#' summary(h)
#'
#' @seealso [check_assumptions()], [auto_transform()]
#' @export
declare_hypothesis <- function(h0, h1, alpha = 0.05) {
  if (!is.character(h0) || length(h0) != 1L || nchar(trimws(h0)) == 0L) {
    rlang::abort("`h0` must be a non-empty character string.")
  }
  if (!is.character(h1) || length(h1) != 1L || nchar(trimws(h1)) == 0L) {
    rlang::abort("`h1` must be a non-empty character string.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    rlang::abort("`alpha` must be a single numeric value strictly between 0 and 1.")
  }

  structure(
    list(
      h0_text       = h0,
      h1_text       = h1,
      alpha         = alpha,
      declared_date = Sys.time()
    ),
    class = "hypothesis"
  )
}

#' Print a `hypothesis` Object
#'
#' Displays the null and alternative hypothesis statements together with the
#' significance level.
#'
#' @param x An object of class `"hypothesis"`.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Returns `x` invisibly.
#'
#' @examples
#' h <- declare_hypothesis("No effect", "There is an effect", alpha = 0.01)
#' print(h)
#'
#' @export
print.hypothesis <- function(x, ...) {
  cat("Hypothesis\n")
  cat(strrep("-", 40L), "\n", sep = "")
  cat(paste0("H\u2080: ", x$h0_text, "\n"))
  cat(paste0("H\u2081: ", x$h1_text, "\n"))
  cat(paste0("Significance level (\u03b1 = ", x$alpha, ")\n"))
  invisible(x)
}

#' Summarise a `hypothesis` Object
#'
#' Prints a detailed summary of the hypothesis including the declaration
#' date and time.
#'
#' @param object An object of class `"hypothesis"`.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Returns `object` invisibly.
#'
#' @examples
#' h <- declare_hypothesis("No effect", "There is an effect")
#' summary(h)
#'
#' @export
summary.hypothesis <- function(object, ...) {
  print(object, ...)
  cat(paste0(
    "Declared: ",
    format(object$declared_date, "%Y-%m-%d %H:%M:%S"),
    "\n"
  ))
  invisible(object)
}

# ============================================================================
# Step 9: Assumption Testing Engine
# ============================================================================

#' Check Statistical Assumptions
#'
#' Runs a battery of assumption tests on a numeric vector.  By default a
#' Shapiro-Wilk normality test is performed.  When a grouping variable is
#' supplied, Levene's test for homogeneity of variance is added.  The result
#' is an S3 object of class `"assumption_report"` containing test statistics,
#' pass/fail flags, and a plain-language recommendation.
#'
#' @param x A numeric vector (NAs are silently removed before testing).
#' @param groups An optional factor or character vector of the same length as
#'   `x` identifying group membership.  Required for Levene's test.
#' @param alpha Numeric significance level (default `0.05`).
#'
#' @return An S3 object of class `"assumption_report"` with components:
#'   \describe{
#'     \item{`test_results`}{Data frame with columns `test_name`,
#'       `statistic`, `p_value`, `passed_at_alpha`.}
#'     \item{`interpretation`}{Named character vector — `"PASS"` or `"FAIL"`
#'       for each test.}
#'     \item{`recommendation`}{Character string with a plain-language
#'       suggestion (e.g., apply a transformation or use a non-parametric
#'       test).}
#'     \item{`alpha`}{The significance level used.}
#'     \item{`n_groups`}{Integer.  Number of distinct groups (or `NA` if no
#'       grouping variable was supplied).}
#'   }
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(30, mean = 5, sd = 1)
#' report <- check_assumptions(x)
#' print(report)
#'
#' # With grouping variable (adds Levene's test)
#' grp <- rep(c("A", "B", "C"), each = 10)
#' report2 <- check_assumptions(x, groups = grp)
#' summary(report2)
#'
#' @seealso [auto_transform()], [suggest_nonparametric()]
#' @importFrom car leveneTest
#' @importFrom stats shapiro.test
#' @export
check_assumptions <- function(x, groups = NULL, alpha = 0.05) {
  # Validate inputs
  if (!is.numeric(x)) {
    rlang::abort("`x` must be a numeric vector.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    rlang::abort("`alpha` must be a single numeric value strictly between 0 and 1.")
  }

  # Remove NAs from x (and align groups)
  if (!is.null(groups)) {
    if (length(groups) != length(x)) {
      rlang::abort("`groups` must have the same length as `x`.")
    }
    keep <- !is.na(x)
    x_clean  <- x[keep]
    grp_clean <- groups[keep]
  } else {
    x_clean   <- x[!is.na(x)]
    grp_clean <- NULL
  }

  n <- length(x_clean)

  # Edge-case: too few observations
  if (n < 3L) {
    rlang::warn(paste0(
      "`x` has only ", n, " non-missing value(s). ",
      "Statistical tests require at least 3 observations. ",
      "Returning NA results."
    ))
    result <- .make_assumption_report(
      test_results = data.frame(
        test_name       = "Shapiro-Wilk",
        statistic       = NA_real_,
        p_value         = NA_real_,
        passed_at_alpha = NA,
        stringsAsFactors = FALSE
      ),
      interpretation = c("Shapiro-Wilk" = "UNKNOWN"),
      recommendation = "Insufficient data for assumption testing (n < 3).",
      alpha          = alpha,
      n_groups       = if (is.null(grp_clean)) NA_integer_ else
        length(unique(grp_clean))
    )
    return(result)
  }

  # Shapiro-Wilk normality test (max 5000 obs)
  sw_input <- if (n > 5000L) x_clean[seq_len(5000L)] else x_clean
  sw <- tryCatch(
    stats::shapiro.test(sw_input),
    error = function(e) {
      rlang::warn(paste("Shapiro-Wilk failed:", conditionMessage(e)))
      list(statistic = NA_real_, p.value = NA_real_)
    }
  )

  results_list <- list(
    data.frame(
      test_name        = "Shapiro-Wilk",
      statistic        = unname(sw$statistic),
      p_value          = sw$p.value,
      passed_at_alpha  = if (is.na(sw$p.value)) NA else (sw$p.value >= alpha),
      stringsAsFactors = FALSE
    )
  )

  interpretation <- c(
    "Shapiro-Wilk" = .pass_fail(sw$p.value, alpha, fail_below = FALSE)
  )

  # Levene's test (only when groups are provided)
  n_groups <- NA_integer_
  if (!is.null(grp_clean)) {
    grp_factor <- as.factor(grp_clean)
    n_groups   <- nlevels(grp_factor)

    if (n_groups < 2L) {
      rlang::warn("Levene's test requires at least 2 groups; skipping.")
    } else {
      lev <- tryCatch(
        car::leveneTest(x_clean, group = grp_factor, center = stats::median),
        error = function(e) {
          rlang::warn(paste("Levene's test failed:", conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(lev)) {
        lev_stat <- lev[["F value"]][1L]
        lev_p    <- lev[["Pr(>F)"]][1L]
        results_list[[2L]] <- data.frame(
          test_name        = "Levene",
          statistic        = lev_stat,
          p_value          = lev_p,
          passed_at_alpha  = if (is.na(lev_p)) NA else (lev_p >= alpha),
          stringsAsFactors = FALSE
        )
        interpretation["Levene"] <- .pass_fail(lev_p, alpha, fail_below = FALSE)
      }
    }
  }

  test_results <- do.call(rbind, results_list)
  rownames(test_results) <- NULL

  # Build recommendation
  normality_ok  <- isTRUE(test_results$passed_at_alpha[
    test_results$test_name == "Shapiro-Wilk"
  ])
  homogeneity_ok <- if (!is.null(grp_clean) &&
    "Levene" %in% test_results$test_name) {
    isTRUE(test_results$passed_at_alpha[test_results$test_name == "Levene"])
  } else {
    TRUE
  }

  recommendation <- .build_recommendation(normality_ok, homogeneity_ok, n_groups)

  .make_assumption_report(
    test_results   = test_results,
    interpretation = interpretation,
    recommendation = recommendation,
    alpha          = alpha,
    n_groups       = n_groups
  )
}

# Helper: create assumption_report
.make_assumption_report <- function(test_results, interpretation,
                                    recommendation, alpha, n_groups) {
  structure(
    list(
      test_results   = test_results,
      interpretation = interpretation,
      recommendation = recommendation,
      alpha          = alpha,
      n_groups       = n_groups
    ),
    class = "assumption_report"
  )
}

# Helper: pass/fail string
.pass_fail <- function(p, alpha, fail_below = TRUE) {
  if (is.na(p)) return("UNKNOWN")
  if (fail_below) {
    if (p < alpha) "FAIL" else "PASS"
  } else {
    if (p >= alpha) "PASS" else "FAIL"
  }
}

# Helper: build plain-language recommendation
.build_recommendation <- function(normality_ok, homogeneity_ok, n_groups) {
  msgs <- character(0L)
  if (!normality_ok) {
    msgs <- c(msgs,
      "Data appear non-normal. Consider applying a transformation (see auto_transform()) or using a non-parametric test (see suggest_nonparametric())."
    )
  }
  if (!homogeneity_ok) {
    msgs <- c(msgs,
      "Variances are not homogeneous across groups. Consider Welch's t-test (2 groups) or Welch ANOVA (3+ groups)."
    )
  }
  if (length(msgs) == 0L) {
    "Assumptions appear to be met. Proceed with standard parametric analysis."
  } else {
    paste(msgs, collapse = " ")
  }
}

#' Print an `assumption_report` Object
#'
#' @param x An object of class `"assumption_report"`.
#' @param ... Further arguments (ignored).
#' @return Returns `x` invisibly.
#' @export
print.assumption_report <- function(x, ...) {
  cat("Assumption Report\n")
  cat(strrep("-", 40L), "\n", sep = "")
  cat(sprintf("Significance level: \u03b1 = %g\n\n", x$alpha))

  cat("Test Results:\n")
  for (i in seq_len(nrow(x$test_results))) {
    row <- x$test_results[i, ]
    cat(sprintf(
      "  %-20s stat = %s  p = %s  [%s]\n",
      row$test_name,
      format_statistic(row$statistic),
      format_p_value(row$p_value),
      if (is.na(row$passed_at_alpha)) "UNKNOWN"
      else if (isTRUE(row$passed_at_alpha)) "PASS" else "FAIL"
    ))
  }

  cat("\nRecommendation:\n")
  cat(strwrap(x$recommendation, width = 70L, indent = 2L,
              exdent = 2L), sep = "\n")
  invisible(x)
}

#' Summarise an `assumption_report` Object
#'
#' @param object An object of class `"assumption_report"`.
#' @param ... Further arguments (ignored).
#' @return Returns `object` invisibly.
#' @export
summary.assumption_report <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

# ============================================================================
# Step 10: Data Transformation Utilities
# ============================================================================

#' Auto-Transform Data for Normality
#'
#' Tests a set of candidate transformations (log, sqrt, inverse, reciprocal,
#' Box-Cox) on a numeric vector and selects the one that yields the highest
#' Shapiro-Wilk p-value.  Returns an S3 object with the best transformation
#' name, before/after p-values, a formula suggestion, and a ggplot2 comparison.
#'
#' @param x A numeric vector.  NAs are removed before analysis.
#' @param response_name Character string used as the variable name in the
#'   `formula_suggestion` output (default `"y"`).
#'
#' @return An S3 object of class `"auto_transform_result"` with components:
#'   \describe{
#'     \item{`best_transform`}{Character. Name of the best transformation.}
#'     \item{`original_p_value`}{Numeric. Shapiro-Wilk p-value on raw data.}
#'     \item{`transformed_p_value`}{Numeric. Shapiro-Wilk p-value after
#'       transformation.}
#'     \item{`lambda`}{Numeric or `NA`. Box-Cox lambda when applicable.}
#'     \item{`formula_suggestion`}{Character. E.g., `"y ~ log(x)"`.}
#'     \item{`all_results`}{Data frame of all transformations tried.}
#'     \item{`plot`}{A `ggplot` object comparing the original and best
#'       transformed distributions.}
#'   }
#'
#' @examples
#' set.seed(1)
#' x <- exp(rnorm(50))   # log-normal data
#' result <- auto_transform(x)
#' print(result)
#'
#' @seealso [check_assumptions()], [suggest_nonparametric()]
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density after_stat
#'   facet_wrap vars labs theme_bw
#' @importFrom rlang .data
#' @importFrom stats shapiro.test optimise
#' @export
auto_transform <- function(x, response_name = "y") {
  if (!is.numeric(x)) rlang::abort("`x` must be a numeric vector.")
  if (!is.character(response_name) || length(response_name) != 1L) {
    rlang::abort("`response_name` must be a single character string.")
  }

  x_clean <- x[is.finite(x)]
  n <- length(x_clean)

  if (n < 3L) {
    rlang::abort("At least 3 finite values are required for auto_transform().")
  }

  # Shapiro-Wilk p-value helper (returns NA on error)
  .sw_p <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 3L) return(NA_real_)
    z_use <- if (length(z) > 5000L) z[seq_len(5000L)] else z
    tryCatch(
      stats::shapiro.test(z_use)$p.value,
      error = function(e) NA_real_
    )
  }

  original_p <- .sw_p(x_clean)

  # Define transformations
  transforms <- list(
    none      = list(fn = identity,        valid = function(v) TRUE),
    log       = list(fn = log,             valid = function(v) all(v > 0)),
    sqrt      = list(fn = sqrt,            valid = function(v) all(v >= 0)),
    inverse   = list(fn = function(v) 1 / v, valid = function(v) all(v != 0)),
    reciprocal = list(fn = function(v) 1 / v, valid = function(v) all(v != 0)),
    boxcox    = list(
      fn = function(v) {
        lam <- .boxcox_lambda(v)
        if (is.na(lam)) return(rep(NA_real_, length(v)))
        if (abs(lam) < 1e-10) log(v) else (v^lam - 1) / lam
      },
      valid = function(v) all(v > 0)
    )
  )
  # Remove duplicate reciprocal/inverse
  transforms[["reciprocal"]] <- NULL

  all_results <- vector("list", length(transforms))
  names(all_results) <- names(transforms)

  for (nm in names(transforms)) {
    tr <- transforms[[nm]]
    if (!tr$valid(x_clean)) {
      all_results[[nm]] <- data.frame(
        transform  = nm,
        p_value    = NA_real_,
        lambda     = NA_real_,
        applicable = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    z <- tryCatch(tr$fn(x_clean), error = function(e) rep(NA_real_, n))
    p <- .sw_p(z)

    lam <- NA_real_
    if (nm == "boxcox") lam <- .boxcox_lambda(x_clean)

    all_results[[nm]] <- data.frame(
      transform  = nm,
      p_value    = p,
      lambda     = lam,
      applicable = TRUE,
      stringsAsFactors = FALSE
    )
  }

  all_df <- do.call(rbind, all_results)
  rownames(all_df) <- NULL

  # Best: highest p-value (excluding "none")
  candidates <- all_df[all_df$applicable & all_df$transform != "none" &
                          !is.na(all_df$p_value), ]
  if (nrow(candidates) == 0L) {
    best_name <- "none"
    best_p    <- original_p
    best_lam  <- NA_real_
  } else {
    idx       <- which.max(candidates$p_value)
    best_name <- candidates$transform[idx]
    best_p    <- candidates$p_value[idx]
    best_lam  <- candidates$lambda[idx]
  }

  # Formula suggestion
  formula_str <- .formula_suggestion(best_name, response_name, best_lam)

  # Build plot
  tr_fn   <- transforms[[best_name]]$fn
  x_trans <- if (best_name == "none") x_clean else
    tryCatch(tr_fn(x_clean), error = function(e) x_clean)

  plot_df <- rbind(
    data.frame(panel = "Original",    value = x_clean,  stringsAsFactors = FALSE),
    data.frame(panel = paste0("Transformed (", best_name, ")"),
               value = x_trans[is.finite(x_trans)],     stringsAsFactors = FALSE)
  )
  plot_df$panel <- factor(plot_df$panel, levels = unique(plot_df$panel))

  p_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["value"]])) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 20L, fill = "steelblue", colour = "white", alpha = 0.6
    ) +
    ggplot2::geom_density(colour = "darkred", linewidth = 0.8) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]]), scales = "free") +
    ggplot2::labs(
      title = paste0("Auto-transform: best = '", best_name, "'"),
      subtitle = sprintf(
        "Original SW p = %s  |  Transformed SW p = %s",
        format_p_value(original_p), format_p_value(best_p)
      ),
      x = "Value", y = "Density"
    ) +
    ggplot2::theme_bw()

  structure(
    list(
      best_transform       = best_name,
      original_p_value     = original_p,
      transformed_p_value  = best_p,
      lambda               = if (is.na(best_lam)) NA_real_ else best_lam,
      formula_suggestion   = formula_str,
      all_results          = all_df,
      plot                 = p_plot
    ),
    class = "auto_transform_result"
  )
}

# Internal: compute Box-Cox lambda via profile log-likelihood
.boxcox_lambda <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 3L) return(NA_real_)
  n        <- length(x)
  log_x_sum <- sum(log(x))

  llik <- function(lam) {
    y <- if (abs(lam) < 1e-10) log(x) else (x^lam - 1) / lam
    mu <- mean(y)
    -n / 2 * log(sum((y - mu)^2) / n) + (lam - 1) * log_x_sum
  }

  opt <- tryCatch(
    stats::optimise(llik, interval = c(-3, 3), maximum = TRUE),
    error = function(e) list(maximum = NA_real_)
  )
  opt$maximum
}

# Internal: build formula suggestion string
.formula_suggestion <- function(transform_name, response_name, lambda) {
  predictor <- switch(
    transform_name,
    none      = "x",
    log       = "log(x)",
    sqrt      = "sqrt(x)",
    inverse   = "1/x",
    boxcox    = {
      if (!is.na(lambda) && abs(lambda) < 1e-10) {
        "log(x)"
      } else if (!is.na(lambda)) {
        sprintf("(x^%.3f - 1) / %.3f", lambda, lambda)
      } else {
        "boxcox(x)"
      }
    },
    "x"
  )
  paste0(response_name, " ~ ", predictor)
}

#' Print an `auto_transform_result` Object
#'
#' @param x An object of class `"auto_transform_result"`.
#' @param ... Further arguments (ignored).
#' @return Returns `x` invisibly.
#' @export
print.auto_transform_result <- function(x, ...) {
  cat("Auto-Transform Result\n")
  cat(strrep("-", 40L), "\n", sep = "")
  cat(sprintf("Best transformation: %s\n", x$best_transform))
  cat(sprintf("Original SW p-value: %s\n", format_p_value(x$original_p_value)))
  cat(sprintf("Transformed SW p-value: %s\n",
              format_p_value(x$transformed_p_value)))
  if (!is.na(x$lambda)) {
    cat(sprintf("Box-Cox lambda: %.4f\n", x$lambda))
  }
  cat(sprintf("Formula suggestion: %s\n", x$formula_suggestion))
  invisible(x)
}

#' Summarise an `auto_transform_result` Object
#'
#' @param object An object of class `"auto_transform_result"`.
#' @param ... Further arguments (ignored).
#' @return Returns `object` invisibly.
#' @export
summary.auto_transform_result <- function(object, ...) {
  print(object, ...)
  cat("\nAll transformations tried:\n")
  print(object$all_results, row.names = FALSE)
  invisible(object)
}

#' Suggest Non-Parametric Alternatives
#'
#' Examines an `assumption_report` and recommends non-parametric (or
#' variance-robust) alternatives when parametric assumptions are violated.
#'
#' @param assumption_report An object of class `"assumption_report"` produced
#'   by [check_assumptions()].
#' @param n_groups Integer.  Number of treatment groups.  If `NULL`
#'   (default), the value stored in `assumption_report$n_groups` is used.
#'
#' @return A named character vector of recommendations.  Names are
#'   `"normality"` and/or `"homogeneity"` for the relevant violations.
#'   When all assumptions pass, a single `"parametric"` element is returned.
#'
#' @examples
#' set.seed(99)
#' x <- c(rexp(20), rnorm(20, 10))   # non-normal
#' grp <- rep(c("A", "B"), each = 20)
#' report <- check_assumptions(x, groups = grp)
#' suggest_nonparametric(report)
#'
#' @seealso [check_assumptions()], [auto_transform()]
#' @export
suggest_nonparametric <- function(assumption_report, n_groups = NULL) {
  if (!inherits(assumption_report, "assumption_report")) {
    rlang::abort("`assumption_report` must be an object of class 'assumption_report'.")
  }

  n_grp <- if (!is.null(n_groups)) {
    as.integer(n_groups)
  } else {
    assumption_report$n_groups
  }

  tr <- assumption_report$test_results
  alpha <- assumption_report$alpha

  normality_failed <- tryCatch({
    sw_row <- tr[tr$test_name == "Shapiro-Wilk", ]
    nrow(sw_row) > 0L && !is.na(sw_row$passed_at_alpha) &&
      !isTRUE(sw_row$passed_at_alpha)
  }, error = function(e) FALSE)

  homogeneity_failed <- tryCatch({
    lev_row <- tr[tr$test_name == "Levene", ]
    nrow(lev_row) > 0L && !is.na(lev_row$passed_at_alpha) &&
      !isTRUE(lev_row$passed_at_alpha)
  }, error = function(e) FALSE)

  recs <- character(0L)

  if (normality_failed) {
    recs["normality"] <- if (!is.na(n_grp) && n_grp == 2L) {
      paste0(
        "Non-normality detected. Recommended: Mann-Whitney U test ",
        "(wilcox.test()) for two-group comparison."
      )
    } else if (!is.na(n_grp) && n_grp >= 3L) {
      paste0(
        "Non-normality detected. Recommended: Kruskal-Wallis test ",
        "(kruskal.test()) for three or more groups."
      )
    } else {
      paste0(
        "Non-normality detected. Recommended: Wilcoxon/Mann-Whitney ",
        "(2 groups) or Kruskal-Wallis (3+ groups)."
      )
    }
  }

  if (homogeneity_failed) {
    recs["homogeneity"] <- if (!is.na(n_grp) && n_grp == 2L) {
      paste0(
        "Unequal variances detected. Recommended: Welch's two-sample ",
        "t-test (t.test(var.equal = FALSE))."
      )
    } else {
      paste0(
        "Unequal variances detected. Recommended: Welch's one-way ANOVA ",
        "(oneway.test(var.equal = FALSE))."
      )
    }
  }

  if (length(recs) == 0L) {
    recs["parametric"] <- paste0(
      "All assumptions met at \u03b1 = ", alpha,
      ". Proceed with standard parametric analysis."
    )
  }

  recs
}

# ============================================================================
# Step 11: Mendelian Segregation Tester
# ============================================================================

#' Chi-Square Test for Mendelian Segregation
#'
#' Tests whether observed offspring counts are consistent with an expected
#' Mendelian segregation ratio using Pearson's chi-square goodness-of-fit test.
#'
#' The chi-square statistic is:
#' \deqn{\chi^2 = \sum_{i} \frac{(O_i - E_i)^2}{E_i}}
#' with degrees of freedom \eqn{df = k - 1} where \eqn{k} is the number of
#' classes.
#'
#' @param observed Numeric vector of observed counts (integers \eqn{\geq 0}).
#' @param expected_ratios Numeric vector of the same length as `observed`
#'   giving the expected ratio (e.g., `c(3, 1)` for a 3:1 ratio).
#' @param alpha Numeric significance level (default `0.05`).
#'
#' @return An S3 object of class `"chi_square_result"` with components:
#'   \describe{
#'     \item{`chi_sq_stat`}{Numeric. The \eqn{\chi^2} statistic.}
#'     \item{`p_value`}{Numeric. P-value from the chi-square distribution.}
#'     \item{`df`}{Integer. Degrees of freedom.}
#'     \item{`observed`}{Numeric vector. Supplied observed counts.}
#'     \item{`expected`}{Numeric vector. Expected counts derived from the
#'       total and the ratios.}
#'     \item{`ratio_string`}{Character. Human-readable ratio string
#'       (e.g., `"3:1"`).}
#'     \item{`interpretation`}{Character. Plain-language conclusion.}
#'     \item{`goodness_of_fit`}{Logical. `TRUE` if the data are consistent
#'       with the expected ratio at the given `alpha`.}
#'     \item{`alpha`}{Numeric. Significance level used.}
#'   }
#'
#' @examples
#' # 3:1 ratio (75 dominant : 25 recessive)
#' result <- chi_square_segregation(
#'   observed        = c(75, 25),
#'   expected_ratios = c(3, 1)
#' )
#' print(result)
#'
#' # 9:3:3:1 dihybrid ratio
#' result2 <- chi_square_segregation(
#'   observed        = c(90, 30, 30, 10),
#'   expected_ratios = c(9, 3, 3, 1)
#' )
#' print(result2)
#'
#' @seealso [check_assumptions()]
#' @importFrom stats pchisq
#' @export
chi_square_segregation <- function(observed, expected_ratios, alpha = 0.05) {
  # Input validation
  if (!is.numeric(observed) || any(is.na(observed)) || any(observed < 0)) {
    rlang::abort("`observed` must be a numeric vector of non-negative counts.")
  }
  if (!is.numeric(expected_ratios) || any(is.na(expected_ratios)) ||
      any(expected_ratios <= 0)) {
    rlang::abort("`expected_ratios` must be a numeric vector of positive values.")
  }
  if (length(observed) < 2L) {
    rlang::abort("`observed` must have at least 2 elements.")
  }
  if (length(observed) != length(expected_ratios)) {
    rlang::abort("`observed` and `expected_ratios` must have the same length.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    rlang::abort("`alpha` must be a single numeric value strictly between 0 and 1.")
  }

  n_total  <- sum(observed)
  ratio_sum <- sum(expected_ratios)
  expected  <- (expected_ratios / ratio_sum) * n_total

  if (any(expected < 1)) {
    rlang::warn(
      "Some expected counts are less than 1. Chi-square approximation may be unreliable."
    )
  } else if (any(expected < 5)) {
    rlang::warn(
      "Some expected counts are less than 5. Consider combining classes or using Fisher's exact test."
    )
  }

  df       <- length(observed) - 1L
  chi_sq   <- sum((observed - expected)^2 / expected)
  p_value  <- stats::pchisq(chi_sq, df = df, lower.tail = FALSE)

  ratio_string <- paste(as.integer(round(expected_ratios)), collapse = ":")
  goodness_of_fit <- p_value >= alpha

  interpretation <- if (goodness_of_fit) {
    sprintf(
      "Data are consistent with a %s ratio at \u03b1 = %g (p = %s).",
      ratio_string, alpha, format_p_value(p_value)
    )
  } else {
    sprintf(
      "Data deviate significantly from a %s ratio at \u03b1 = %g (p = %s).",
      ratio_string, alpha, format_p_value(p_value)
    )
  }

  structure(
    list(
      chi_sq_stat     = chi_sq,
      p_value         = p_value,
      df              = df,
      observed        = observed,
      expected        = expected,
      ratio_string    = ratio_string,
      interpretation  = interpretation,
      goodness_of_fit = goodness_of_fit,
      alpha           = alpha
    ),
    class = "chi_square_result"
  )
}

#' Print a `chi_square_result` Object
#'
#' @param x An object of class `"chi_square_result"`.
#' @param ... Further arguments (ignored).
#' @return Returns `x` invisibly.
#' @export
print.chi_square_result <- function(x, ...) {
  cat("Chi-Square Segregation Test\n")
  cat(strrep("-", 40L), "\n", sep = "")
  cat(sprintf("Expected ratio: %s\n", x$ratio_string))
  cat(sprintf("Observed:  %s\n", paste(x$observed, collapse = ", ")))
  cat(sprintf("Expected:  %s\n",
              paste(round(x$expected, 2L), collapse = ", ")))
  cat(sprintf(
    "\u03c7\u00b2(%d) = %s,  p = %s\n",
    x$df,
    format_statistic(x$chi_sq_stat),
    format_p_value(x$p_value)
  ))
  cat("\n", x$interpretation, "\n", sep = "")
  invisible(x)
}

# ============================================================================
# Step 12: Data Import Helpers
# ============================================================================

#' Import Field Trial Data
#'
#' Reads field trial data from a variety of file formats (CSV, TSV, Excel,
#' GeoJSON, Shapefile) with automatic type detection, missing-value
#' standardisation, and optional spatial object creation.
#'
#' @param file Character. Path to the data file.  Supported extensions:
#'   `.csv`, `.tsv`, `.txt` (delimited text), `.xlsx`, `.xls` (Excel),
#'   `.geojson`, `.json` (GeoJSON), `.shp` (Shapefile).
#' @param sheet Integer or character.  Sheet index or name for Excel files
#'   (default `1`).
#' @param na_strings Character vector of strings to treat as `NA`
#'   (default `c("NA", "", ".", "N/A")`).
#' @param ... Additional arguments forwarded to the underlying reader
#'   (`readr::read_csv()`, `readxl::read_excel()`, or `sf::st_read()`).
#'
#' @return
#' For delimited and Excel files: a [`tibble`][tibble::tibble] with
#' an `"agristat_metadata"` attribute (a data frame of column types).
#' Columns are coerced to numeric where possible; columns matching common
#' date patterns are parsed as `Date`; remaining character columns are
#' left as character.
#'
#' For spatial files (GeoJSON, Shapefile): an `sf` object.
#'
#' @examples
#' # Write a temporary CSV and import it
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(
#'   data.frame(plot = 1:4, yield = c(3.2, 4.1, 3.8, 4.5), trt = c("A","B","A","B")),
#'   tmp, row.names = FALSE
#' )
#' df <- import_field_data(tmp)
#' print(df)
#' attr(df, "agristat_metadata")
#'
#' @seealso [check_assumptions()], [declare_hypothesis()]
#' @importFrom readr read_csv
#' @importFrom utils read.csv
#' @importFrom tools file_ext
#' @export
import_field_data <- function(file, sheet = 1, na_strings = c("NA", "", ".", "N/A"), ...) {
  if (!is.character(file) || length(file) != 1L) {
    rlang::abort("`file` must be a single character string.")
  }
  if (!file.exists(file)) {
    rlang::abort(sprintf("File not found: '%s'", file))
  }

  ext <- tolower(tools::file_ext(file))

  # Dispatch by file type
  if (ext %in% c("csv", "tsv", "txt")) {
    result <- .import_delimited(file, na_strings = na_strings, ...)
  } else if (ext %in% c("xlsx", "xls")) {
    result <- .import_excel(file, sheet = sheet, na_strings = na_strings, ...)
  } else if (ext %in% c("geojson", "json", "shp", "gpkg")) {
    result <- .import_spatial(file, ...)
  } else {
    rlang::abort(sprintf(
      "Unsupported file extension '.%s'. Supported: csv, tsv, txt, xlsx, xls, geojson, json, shp, gpkg.",
      ext
    ))
  }

  result
}

# Internal: import delimited text files
.import_delimited <- function(file, na_strings, ...) {
  # Check for empty file
  file_size <- file.info(file)$size
  if (!is.na(file_size) && file_size == 0L) {
    rlang::abort(sprintf("File is empty: '%s'", file))
  }

  df <- tryCatch(
    readr::read_csv(
      file,
      na          = na_strings,
      show_col_types = FALSE,
      ...
    ),
    error = function(e) {
      # Fallback to base R for non-UTF8 or delimiter issues
      tryCatch(
        utils::read.csv(file, na.strings = na_strings,
                        fileEncoding = "latin1", ...),
        error = function(e2) {
          rlang::abort(sprintf("Failed to read file '%s': %s", file,
                               conditionMessage(e2)))
        }
      )
    }
  )

  if (nrow(df) == 0L) {
    rlang::warn(sprintf("File '%s' contains no data rows.", file))
  }

  .attach_metadata(df)
}

# Internal: import Excel files
.import_excel <- function(file, sheet, na_strings, ...) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    rlang::abort(
      "Package 'readxl' is required to import Excel files. ",
      "Install it with: install.packages('readxl')"
    )
  }

  df <- tryCatch(
    readxl::read_excel(file, sheet = sheet, na = na_strings, ...),
    error = function(e) {
      rlang::abort(sprintf("Failed to read Excel file '%s': %s", file,
                           conditionMessage(e)))
    }
  )

  if (nrow(df) == 0L) {
    rlang::warn(sprintf("Excel file '%s' contains no data rows.", file))
  }

  .attach_metadata(df)
}

# Internal: import spatial files (GeoJSON, Shapefile)
.import_spatial <- function(file, ...) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    rlang::abort(
      "Package 'sf' is required to import spatial files. ",
      "Install it with: install.packages('sf')"
    )
  }
  tryCatch(
    sf::st_read(file, quiet = TRUE, ...),
    error = function(e) {
      rlang::abort(sprintf("Failed to read spatial file '%s': %s", file,
                           conditionMessage(e)))
    }
  )
}

# Internal: attach column-type metadata as an attribute
.attach_metadata <- function(df) {
  meta <- data.frame(
    column     = names(df),
    type       = vapply(df, function(col) class(col)[1L], character(1L)),
    n_missing  = vapply(df, function(col) sum(is.na(col)), integer(1L)),
    stringsAsFactors = FALSE
  )
  attr(df, "agristat_metadata") <- meta
  df
}
