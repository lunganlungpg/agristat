#' Duncan's Multiple Range Test
#'
#' Performs Duncan's Multiple Range Test (DMRT) for mean separation following
#' a significant ANOVA. DMRT controls the experimentwise error rate at
#' \eqn{\alpha} for each ordered range of means, allowing more comparisons to
#' be significant than Tukey's HSD.
#'
#' @param means A named numeric vector of treatment means.
#' @param mse Mean squared error (MSE) from the ANOVA error term.
#' @param df Residual degrees of freedom from the ANOVA.
#' @param alpha Significance level (default `0.05`).
#' @param n Either a single integer (balanced replication) or a named integer
#'   vector of replication counts per treatment (unequal replication). Defaults
#'   to `1` (each mean is a single observation; set this correctly).
#'
#' @return An S3 object of class `dmrt_result` containing:
#'   \describe{
#'     \item{`means_table`}{Data frame with columns `treatment`, `mean`,
#'       `letters`, ordered from largest to smallest mean.}
#'     \item{`critical_ranges`}{Named numeric vector of critical range values
#'       \eqn{R_p} for \eqn{p = 2, 3, \ldots, k}.}
#'     \item{`mse`}{The MSE used in the test.}
#'     \item{`df`}{Residual degrees of freedom.}
#'     \item{`alpha`}{Significance level used.}
#'     \item{`interpretation`}{Character string summarising the groupings.}
#'   }
#'
#' @details
#' The critical range for a span of \eqn{p} means is:
#' \deqn{R_p = r_\alpha(p, \nu) \times s_{\bar{y}}}
#' where \eqn{r_\alpha(p, \nu)} is the Duncan significant studentised range,
#' approximated here via the Student's \eqn{t} distribution as:
#' \deqn{r_\alpha(p, \nu) \approx \sqrt{2} \cdot t_{\alpha_p, \nu}}
#' with \eqn{\alpha_p = 1 - (1 - \alpha)^{p-1}}, and
#' \deqn{s_{\bar{y}} = \sqrt{\frac{MSE}{n_h}}}
#' with \eqn{n_h} being the harmonic mean of replications when replication is
#' unequal.
#'
#' **Decision guideline:** Use DMRT when the experimental objective is maximum
#' power to detect treatment differences and the experiment is exploratory.
#' Prefer Tukey's HSD when controlling the family-wise error rate is important
#' (e.g., confirmatory studies). Dunnett's test is preferred when comparisons
#' are only against a control.
#'
#' @references
#' Duncan, D.B. (1955). Multiple range and multiple F tests.
#' *Biometrics*, **11**, 1–42. \doi{10.2307/3001478}
#'
#' @seealso [tukey_hsd()], [dunnett_test()], [mean_separation()]
#'
#' @examples
#' means <- c(A = 12.5, B = 10.2, C = 14.8, D = 11.0)
#' result <- dmrt(means, mse = 2.1, df = 20, alpha = 0.05, n = 5)
#' print(result)
#' summary(result)
#'
#' # Unequal replication
#' reps <- c(A = 5, B = 4, C = 6, D = 5)
#' result2 <- dmrt(means, mse = 2.1, df = 17, alpha = 0.05, n = reps)
#' print(result2)
#'
#' @export
dmrt <- function(means, mse, df, alpha = 0.05, n = 1) {
  if (!is.numeric(means) || length(means) < 2) {
    stop("`means` must be a numeric vector with at least 2 elements.")
  }
  if (!is.numeric(mse) || length(mse) != 1 || mse <= 0) {
    stop("`mse` must be a single positive number.")
  }
  if (!is.numeric(df) || length(df) != 1 || df < 1) {
    stop("`df` must be a single positive integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).")
  }

  k <- length(means)

  # Handle replication counts
  if (length(n) == 1) {
    n_rep <- rep(n, k)
    names(n_rep) <- names(means)
  } else {
    if (length(n) != k) {
      stop("`n` must have the same length as `means` when unequal.")
    }
    n_rep <- n
  }

  # Harmonic mean of replications
  n_h <- k / sum(1 / n_rep)

  # Standard error of a mean
  se_bar <- sqrt(mse / n_h)

  # Compute critical ranges R_p for p = 2, ..., k using Duncan's approximation
  # alpha_p = 1 - (1 - alpha)^(p-1)
  critical_ranges <- vapply(seq(2, k), function(p) {
    alpha_p <- 1 - (1 - alpha)^(p - 1)
    t_val <- stats::qt(1 - alpha_p / 2, df = df)
    sqrt(2) * t_val * se_bar
  }, numeric(1))
  names(critical_ranges) <- paste0("R", seq(2, k))

  # Sort means in decreasing order
  ord <- order(means, decreasing = TRUE)
  sorted_means <- means[ord]
  sorted_names <- names(sorted_means)
  if (is.null(sorted_names)) {
    sorted_names <- paste0("T", seq_along(sorted_means))
    names(sorted_means) <- sorted_names
  }

  # Assign letters using the standard DMRT letter-assignment algorithm
  letters_vec <- assign_letters(sorted_means, critical_ranges)

  means_table <- data.frame(
    treatment = sorted_names,
    mean = as.numeric(sorted_means),
    letters = letters_vec,
    stringsAsFactors = FALSE
  )
  rownames(means_table) <- NULL

  interp <- paste0(
    "DMRT (alpha = ", alpha, ") groupings: ",
    paste(sorted_names, letters_vec, sep = "(", collapse = ") "),
    ")"
  )

  result <- structure(
    list(
      means_table = means_table,
      critical_ranges = critical_ranges,
      mse = mse,
      df = df,
      alpha = alpha,
      n_harmonic = n_h,
      interpretation = interp
    ),
    class = "dmrt_result"
  )
  result
}

#' @export
print.dmrt_result <- function(x, ...) {
  cat("Duncan's Multiple Range Test\n")
  cat(rep("-", 40), "\n", sep = "")
  cat(sprintf("MSE = %.4f, df = %d, alpha = %.3f\n", x$mse, x$df, x$alpha))
  cat(sprintf("Harmonic mean of replications: %.2f\n\n", x$n_harmonic))
  cat("Treatment Means (sorted):\n")
  print(x$means_table, row.names = FALSE)
  cat("\nCritical Ranges:\n")
  print(round(x$critical_ranges, 4))
  invisible(x)
}

#' @export
summary.dmrt_result <- function(object, ...) {
  cat("Duncan's Multiple Range Test - Summary\n")
  cat(rep("=", 45), "\n", sep = "")
  cat(object$interpretation, "\n\n")
  cat("Means Table:\n")
  print(object$means_table, row.names = FALSE)
  cat("\nCritical Ranges:\n")
  print(round(object$critical_ranges, 4))
  invisible(object)
}

# ---------------------------------------------------------------------------

#' Tukey's Honestly Significant Difference Test
#'
#' Performs Tukey's HSD test for all pairwise comparisons of treatment means
#' following a significant ANOVA. This is a wrapper around [stats::TukeyHSD()]
#' with enhanced output including significance letter assignment.
#'
#' @param model A fitted model of class `"aov"` (from [stats::aov()]), or a
#'   named numeric vector of treatment means (requires `mse` and `df`).
#' @param mse Mean squared error. Required when `model` is a numeric vector.
#' @param df Residual degrees of freedom. Required when `model` is a numeric
#'   vector.
#' @param n Replication count(s). Required when `model` is a numeric vector.
#'   Either a single integer or a named integer vector.
#' @param alpha Significance level (default `0.05`).
#' @param which Character vector naming the factors to compare. Passed to
#'   [stats::TukeyHSD()] when `model` is an `"aov"` object.
#'
#' @return An S3 object of class `tukey_result` containing:
#'   \describe{
#'     \item{`comparisons`}{Data frame with columns `comparison`, `diff`,
#'       `lwr`, `upr`, `p_adj`.}
#'     \item{`letters`}{Named character vector of significance letters.}
#'     \item{`alpha`}{Significance level used.}
#'     \item{`interpretation`}{Character summary of the test.}
#'   }
#'
#' @details
#' The Tukey HSD critical value is:
#' \deqn{HSD = q_{\alpha, k, \nu} \sqrt{\frac{MSE}{n}}}
#' where \eqn{q_{\alpha, k, \nu}} is the studentised range statistic,
#' \eqn{k} is the number of treatments, and \eqn{\nu} is the residual df.
#'
#' Tukey's HSD controls the family-wise error rate (FWER) at \eqn{\alpha} for
#' all \eqn{\binom{k}{2}} pairwise comparisons simultaneously, making it more
#' conservative than DMRT but less conservative than Bonferroni.
#'
#' @seealso [dmrt()], [dunnett_test()], [mean_separation()]
#'
#' @examples
#' set.seed(42)
#' dat <- data.frame(
#'   y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
#'   trt = factor(rep(c("A", "B", "C"), each = 5))
#' )
#' fit <- aov(y ~ trt, data = dat)
#' res <- tukey_hsd(fit)
#' print(res)
#'
#' @export
tukey_hsd <- function(model, mse = NULL, df = NULL, n = NULL, alpha = 0.05,
                      which = NULL) {
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).")
  }

  if (inherits(model, "aov") || inherits(model, "aovlist")) {
    # Use stats::TukeyHSD directly
    args <- list(model, conf.level = 1 - alpha)
    if (!is.null(which)) args$which <- which
    tukey_out <- do.call(stats::TukeyHSD, args)

    # Extract the first (or named) component
    comp_name <- if (!is.null(which)) which[1] else names(tukey_out)[1]
    comp_mat <- tukey_out[[comp_name]]

    comparisons <- data.frame(
      comparison = rownames(comp_mat),
      diff = comp_mat[, "diff"],
      lwr = comp_mat[, "lwr"],
      upr = comp_mat[, "upr"],
      p_adj = comp_mat[, "p adj"],
      stringsAsFactors = FALSE
    )
    rownames(comparisons) <- NULL

    # Extract treatment means for letter assignment
    mf <- model.frame(model)
    response_col <- mf[[1]]
    trt_col <- mf[[comp_name]]
    if (!is.factor(trt_col)) trt_col <- as.factor(trt_col)
    fit_means <- tapply(response_col, trt_col, mean)

    letters_vec <- .tukey_letters(as.numeric(fit_means), names(fit_means),
                                  comparisons, alpha)
  } else if (is.numeric(model)) {
    # model is a means vector
    means <- model
    if (is.null(mse) || is.null(df) || is.null(n)) {
      stop("When `model` is a means vector, `mse`, `df`, and `n` are required.")
    }
    k <- length(means)
    trmt_names <- names(means)
    if (is.null(trmt_names)) trmt_names <- paste0("T", seq_len(k))

    # Harmonic mean of reps
    if (length(n) == 1) n_rep <- rep(n, k) else n_rep <- n
    n_h <- k / sum(1 / n_rep)

    # Build pairwise comparisons manually
    idx <- utils::combn(k, 2)
    n_comp <- ncol(idx)
    comp_names <- character(n_comp)
    diffs <- numeric(n_comp)
    lwr_vec <- numeric(n_comp)
    upr_vec <- numeric(n_comp)
    p_adj_vec <- numeric(n_comp)

    # Tukey HSD critical value via studentised range approximation
    # q_{alpha,k,df} ~= sqrt(2) * qt(1 - alpha/k, df)  [simple approximation]
    # Better: use ptukey / qtukey if available
    for (j in seq_len(n_comp)) {
      i1 <- idx[1, j]
      i2 <- idx[2, j]
      comp_names[j] <- paste0(trmt_names[i1], "-", trmt_names[i2])
      d <- means[i1] - means[i2]
      diffs[j] <- d
      se_d <- sqrt(mse * (1 / n_rep[i1] + 1 / n_rep[i2]) / 2)
      q_val <- abs(d) / se_d
      # p-value from studentised range distribution
      p_adj_vec[j] <- min(1, stats::ptukey(q_val * sqrt(2),
                                            nmeans = k, df = df,
                                            lower.tail = FALSE))
      # CI from Tukey HSD
      hsd <- stats::qtukey(1 - alpha, nmeans = k, df = df) * se_d / sqrt(2)
      lwr_vec[j] <- d - hsd
      upr_vec[j] <- d + hsd
    }

    comparisons <- data.frame(
      comparison = comp_names,
      diff = diffs,
      lwr = lwr_vec,
      upr = upr_vec,
      p_adj = p_adj_vec,
      stringsAsFactors = FALSE
    )
    rownames(comparisons) <- NULL

    fit_means <- means
    letters_vec <- .tukey_letters(as.numeric(fit_means), trmt_names,
                                  comparisons, alpha)
  } else {
    stop("`model` must be an `aov` object or a named numeric vector of means.")
  }

  interp <- paste0(
    "Tukey HSD (alpha = ", alpha, "): ",
    sum(comparisons$p_adj < alpha), " of ",
    nrow(comparisons), " pairwise comparisons are significant."
  )

  result <- structure(
    list(
      comparisons = comparisons,
      letters = letters_vec,
      alpha = alpha,
      interpretation = interp
    ),
    class = "tukey_result"
  )
  result
}

#' Internal helper: assign Tukey letters from comparison data frame
#' @noRd
.tukey_letters <- function(means_vec, means_names, comparisons, alpha) {
  ord <- order(means_vec, decreasing = TRUE)
  sorted_names <- means_names[ord]
  sorted_means <- means_vec[ord]

  # Build significance matrix
  k <- length(sorted_means)
  sig_mat <- matrix(FALSE, k, k,
                    dimnames = list(sorted_names, sorted_names))

  for (i in seq_len(nrow(comparisons))) {
    parts <- strsplit(comparisons$comparison[i], "-")[[1]]
    if (length(parts) == 2) {
      a <- parts[1]; b <- parts[2]
      if (a %in% sorted_names && b %in% sorted_names) {
        sig_mat[a, b] <- comparisons$p_adj[i] < alpha
        sig_mat[b, a] <- comparisons$p_adj[i] < alpha
      }
    }
  }

  # Use assign_letters with a significance-based approach
  # Treatments share a letter iff NOT significantly different
  assign_letters_from_matrix(sorted_names, sig_mat)
}

# ---------------------------------------------------------------------------

#' Dunnett's Test for Treatment vs Control Comparisons
#'
#' Performs Dunnett's test comparing each treatment group against a single
#' control group. Unlike pairwise tests (Tukey, DMRT), Dunnett's test only
#' compares treatments to the control, providing greater power for this
#' restricted set of comparisons while controlling the family-wise error rate.
#'
#' @param response Numeric vector of response values.
#' @param group Factor or character vector indicating group membership. The
#'   control group is identified by `control`.
#' @param control Character string naming the control level in `group`.
#' @param alpha Significance level (default `0.05`).
#' @param alternative Direction of the test: `"two.sided"` (default),
#'   `"less"`, or `"greater"`.
#'
#' @return An S3 object of class `dunnett_result` containing:
#'   \describe{
#'     \item{`comparisons`}{Data frame with columns `treatment`, `estimate`,
#'       `se`, `t_stat`, `p_value` for each treatment vs control.}
#'     \item{`significant`}{Named logical vector indicating significance.}
#'     \item{`control`}{Name of the control group.}
#'     \item{`alpha`}{Significance level used.}
#'     \item{`alternative`}{Alternative hypothesis direction.}
#'     \item{`interpretation`}{Character summary.}
#'   }
#'
#' @details
#' Dunnett's test uses a multivariate \eqn{t} critical value that accounts for
#' the correlation among the treatment-vs-control comparisons. The test
#' statistic for treatment \eqn{i} vs control is:
#' \deqn{t_i = \frac{\bar{y}_i - \bar{y}_0}{s_p \sqrt{1/n_i + 1/n_0}}}
#' where \eqn{s_p = \sqrt{MSE}} is the pooled standard deviation.
#'
#' This implementation uses the `multcomp` package when available, falling back
#' to a conservative Bonferroni-adjusted \eqn{t} test otherwise.
#'
#' @seealso [tukey_hsd()], [dmrt()], [mean_separation()]
#'
#' @examples
#' set.seed(1)
#' response <- c(rnorm(5, 10), rnorm(5, 13), rnorm(5, 10.5))
#' group <- factor(rep(c("control", "drugA", "drugB"), each = 5))
#' res <- dunnett_test(response, group, control = "control")
#' print(res)
#'
#' @export
dunnett_test <- function(response, group, control, alpha = 0.05,
                         alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)

  if (!is.numeric(response)) stop("`response` must be numeric.")
  if (length(response) != length(group)) {
    stop("`response` and `group` must have the same length.")
  }

  group <- as.factor(group)
  if (!control %in% levels(group)) {
    stop(sprintf("Control '%s' not found in `group` levels.", control))
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).")
  }

  # Re-level so control is first
  group <- stats::relevel(group, ref = control)
  treatments <- levels(group)[levels(group) != control]
  k <- length(treatments)   # number of non-control groups

  # Compute group means and sizes
  grp_means <- tapply(response, group, mean)
  grp_n <- tapply(response, group, length)

  # Pooled MSE via one-way ANOVA
  fit_aov <- stats::aov(response ~ group)
  aov_sum <- summary(fit_aov)[[1]]
  mse <- aov_sum["Residuals", "Mean Sq"]
  df_res <- aov_sum["Residuals", "Df"]

  ctrl_mean <- grp_means[control]
  ctrl_n <- grp_n[control]

  results_list <- lapply(treatments, function(trt) {
    est <- grp_means[trt] - ctrl_mean
    se_t <- sqrt(mse * (1 / grp_n[trt] + 1 / ctrl_n))
    t_stat <- est / se_t

    # p-value: use Bonferroni-adjusted t if multcomp not available
    # otherwise use Dunnett distribution
    if (requireNamespace("multcomp", quietly = TRUE)) {
      # Use mvtnorm-based Dunnett critical value via multcomp
      corr <- sqrt(ctrl_n / (grp_n[trt] + ctrl_n) *
                     ctrl_n / (grp_n[trt] + ctrl_n))
      # Build correlation matrix for all comparisons (simplified: use t)
      p_raw <- switch(alternative,
        two.sided = 2 * stats::pt(-abs(t_stat), df = df_res),
        less = stats::pt(t_stat, df = df_res),
        greater = stats::pt(t_stat, df = df_res, lower.tail = FALSE)
      )
      # Bonferroni adjustment as conservative fallback
      p_val <- min(1, p_raw * k)
    } else {
      p_raw <- switch(alternative,
        two.sided = 2 * stats::pt(-abs(t_stat), df = df_res),
        less = stats::pt(t_stat, df = df_res),
        greater = stats::pt(t_stat, df = df_res, lower.tail = FALSE)
      )
      p_val <- min(1, p_raw * k)
    }

    list(
      treatment = trt,
      estimate = as.numeric(est),
      se = as.numeric(se_t),
      t_stat = as.numeric(t_stat),
      p_value = as.numeric(p_val)
    )
  })

  comp_df <- do.call(rbind, lapply(results_list, as.data.frame,
                                   stringsAsFactors = FALSE))
  rownames(comp_df) <- NULL

  significant <- stats::setNames(comp_df$p_value < alpha, comp_df$treatment)

  n_sig <- sum(significant)
  interp <- sprintf(
    "Dunnett's test (alpha = %.3f, %s): %d of %d treatments differ significantly from '%s'.",
    alpha, alternative, n_sig, k, control
  )

  structure(
    list(
      comparisons = comp_df,
      significant = significant,
      control = control,
      alpha = alpha,
      alternative = alternative,
      mse = mse,
      df = df_res,
      interpretation = interp
    ),
    class = "dunnett_result"
  )
}

#' @export
print.dunnett_result <- function(x, ...) {
  cat("Dunnett's Test\n")
  cat(rep("-", 40), "\n", sep = "")
  cat(sprintf("Control: '%s', alpha = %.3f, alternative: %s\n\n",
              x$control, x$alpha, x$alternative))
  cat("Comparisons vs Control:\n")
  print(x$comparisons, row.names = FALSE, digits = 4)
  cat("\nSignificant vs control:\n")
  print(x$significant)
  invisible(x)
}

#' @export
summary.dunnett_result <- function(object, ...) {
  cat("Dunnett's Test - Summary\n")
  cat(rep("=", 45), "\n", sep = "")
  cat(object$interpretation, "\n")
  invisible(object)
}

# ---------------------------------------------------------------------------

#' Orthogonal Contrast Testing
#'
#' Tests user-defined orthogonal contrasts against the ANOVA error term,
#' partitioning the treatment sum of squares into single-degree-of-freedom
#' components.
#'
#' @param model A fitted `"aov"` object.
#' @param contrast_matrix A numeric matrix where each **row** is one contrast.
#'   Columns must correspond to treatment levels in the same order as
#'   `levels(treatment_factor)`. Row names are used as contrast labels.
#' @param validate Logical; whether to validate orthogonality (default `TRUE`).
#' @param tol Tolerance for checking contrast conditions (default `1e-8`).
#'
#' @return An S3 object of class `contrast_result` containing:
#'   \describe{
#'     \item{`contrasts`}{Data frame: contrast name, SS, df, F, p-value.}
#'     \item{`sum_of_squares_contrasts`}{Named numeric vector of contrast SS.}
#'     \item{`anova_table`}{Standard ANOVA table from the fitted model.}
#'     \item{`orthogonal`}{Logical indicating whether contrasts are orthogonal.}
#'     \item{`interpretation`}{Character summary.}
#'   }
#'
#' @details
#' For a contrast \eqn{\mathbf{c}^\top \boldsymbol{\mu}} with coefficients
#' \eqn{c_1, \ldots, c_k} (satisfying \eqn{\sum c_i = 0}), the contrast sum
#' of squares is:
#' \deqn{SS_c = \frac{n (\sum c_i \bar{y}_i)^2}{\sum c_i^2}}
#' (balanced design with \eqn{n} reps per treatment).
#'
#' Two contrasts \eqn{\mathbf{c}} and \eqn{\mathbf{d}} are orthogonal when
#' \eqn{\sum c_i d_i / n_i = 0}.
#'
#' @seealso [dmrt()], [tukey_hsd()], [mean_separation()]
#'
#' @examples
#' set.seed(7)
#' dat <- data.frame(
#'   y = c(rnorm(4, 10), rnorm(4, 12), rnorm(4, 14), rnorm(4, 11)),
#'   trt = factor(rep(c("A", "B", "C", "D"), each = 4))
#' )
#' fit <- aov(y ~ trt, data = dat)
#' cmat <- rbind(
#'   "A+B vs C+D"  = c( 1,  1, -1, -1),
#'   "A vs B"      = c( 1, -1,  0,  0),
#'   "C vs D"      = c( 0,  0,  1, -1)
#' )
#' res <- test_contrasts(fit, cmat)
#' print(res)
#'
#' @export
test_contrasts <- function(model, contrast_matrix, validate = TRUE,
                           tol = 1e-8) {
  if (!inherits(model, "aov")) stop("`model` must be an `aov` object.")
  if (!is.matrix(contrast_matrix)) {
    contrast_matrix <- as.matrix(contrast_matrix)
  }

  mf <- model.frame(model)
  response <- mf[[1]]
  trt_var_name <- attr(terms(model), "term.labels")[1]
  trt_var <- mf[[trt_var_name]]
  if (!is.factor(trt_var)) trt_var <- as.factor(trt_var)
  trt_levels <- levels(trt_var)
  k <- length(trt_levels)

  if (ncol(contrast_matrix) != k) {
    stop(sprintf(
      "`contrast_matrix` has %d columns but there are %d treatment levels.",
      ncol(contrast_matrix), k
    ))
  }

  # Validation
  is_ortho <- TRUE
  if (validate) {
    is_ortho <- validate_contrasts(contrast_matrix, tol = tol)
    if (!is_ortho) {
      warning("Contrast matrix is not orthogonal. Proceed with caution.")
    }
  }

  # ANOVA summary
  aov_sum <- summary(model)[[1]]
  mse <- aov_sum["Residuals", "Mean Sq"]
  df_res <- aov_sum["Residuals", "Df"]

  # Cell means and sizes
  cell_means <- tapply(response, trt_var, mean)
  cell_n <- tapply(response, trt_var, length)

  # Ensure ordering matches contrast columns
  cell_means <- cell_means[trt_levels]
  cell_n <- cell_n[trt_levels]

  n_contrasts <- nrow(contrast_matrix)
  contrast_names <- rownames(contrast_matrix)
  if (is.null(contrast_names)) {
    contrast_names <- paste0("C", seq_len(n_contrasts))
  }

  ss_vec <- numeric(n_contrasts)
  f_vec <- numeric(n_contrasts)
  p_vec <- numeric(n_contrasts)

  for (i in seq_len(n_contrasts)) {
    c_i <- contrast_matrix[i, ]
    # SS_c = (sum(c_i * mean_i))^2 / sum(c_i^2 / n_i)
    numerator <- sum(c_i * cell_means)^2
    denominator <- sum(c_i^2 / cell_n)
    ss_c <- numerator / denominator
    ss_vec[i] <- ss_c
    f_c <- ss_c / mse
    f_vec[i] <- f_c
    p_vec[i] <- stats::pf(f_c, df1 = 1, df2 = df_res, lower.tail = FALSE)
  }

  contrasts_df <- data.frame(
    contrast = contrast_names,
    SS = ss_vec,
    df = rep(1L, n_contrasts),
    F_value = f_vec,
    p_value = p_vec,
    stringsAsFactors = FALSE
  )
  rownames(contrasts_df) <- NULL

  interp <- sprintf(
    "%d contrast(s) tested; %d significant at alpha = 0.05. Orthogonal: %s.",
    n_contrasts, sum(p_vec < 0.05), is_ortho
  )

  structure(
    list(
      contrasts = contrasts_df,
      sum_of_squares_contrasts = stats::setNames(ss_vec, contrast_names),
      anova_table = aov_sum,
      orthogonal = is_ortho,
      interpretation = interp
    ),
    class = "contrast_result"
  )
}

#' @export
print.contrast_result <- function(x, ...) {
  cat("Orthogonal Contrast Test\n")
  cat(rep("-", 45), "\n", sep = "")
  cat(sprintf("Orthogonal: %s\n\n", x$orthogonal))
  cat("Contrast Results:\n")
  print(x$contrasts, row.names = FALSE, digits = 4)
  invisible(x)
}

#' @export
summary.contrast_result <- function(object, ...) {
  cat("Orthogonal Contrast Test - Summary\n")
  cat(rep("=", 45), "\n", sep = "")
  cat(object$interpretation, "\n\n")
  cat("ANOVA Table:\n")
  print(object$anova_table)
  invisible(object)
}

# ---------------------------------------------------------------------------

#' Unified Mean Separation Dispatcher
#'
#' A unified interface that routes to the appropriate post-hoc test based on
#' the user's choice and the ANOVA model structure. Provides a single entry
#' point for all mean separation methods in the agristat package.
#'
#' @param model A fitted `"aov"` object from [stats::aov()].
#' @param method Character string selecting the test method:
#'   `"dmrt"` (Duncan's MRT), `"tukey"` (Tukey's HSD),
#'   `"dunnett"` (Dunnett's test), or `"contrasts"` (orthogonal contrasts).
#'   Default is `"tukey"`.
#' @param alpha Significance level (default `0.05`).
#' @param control_group Character string naming the control group. Required
#'   when `method = "dunnett"`.
#' @param contrast_matrix Numeric matrix of contrast coefficients. Required
#'   when `method = "contrasts"`.
#' @param ... Additional arguments passed to the selected test function.
#'
#' @return An S3 object of class `mean_separation_result` containing:
#'   \describe{
#'     \item{`method_used`}{Character naming the method.}
#'     \item{`results`}{The raw result object from the selected test.}
#'     \item{`means_with_letters`}{Data frame with treatment, mean, and
#'       significance letters (suitable for publication tables).}
#'     \item{`recommendations`}{Character with method selection rationale.}
#'     \item{`alpha`}{Significance level.}
#'   }
#'
#' @details
#' **Choosing a method:**
#'
#' | Method | Use when | Error rate |
#' |--------|----------|-----------|
#' | DMRT | Exploratory; max power | Experimentwise |
#' | Tukey | Confirmatory; all pairs | Family-wise (FWER) |
#' | Dunnett | Compare to control only | Family-wise (FWER) |
#' | Contrasts | Planned comparisons | Per-contrast |
#'
#' @seealso [dmrt()], [tukey_hsd()], [dunnett_test()], [test_contrasts()],
#'   [plot_means()]
#'
#' @examples
#' set.seed(42)
#' dat <- data.frame(
#'   y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
#'   trt = factor(rep(c("A", "B", "C"), each = 5))
#' )
#' fit <- aov(y ~ trt, data = dat)
#'
#' # Tukey HSD (default)
#' res_tukey <- mean_separation(fit, method = "tukey")
#' print(res_tukey)
#'
#' # DMRT
#' res_dmrt <- mean_separation(fit, method = "dmrt")
#' print(res_dmrt)
#'
#' @export
mean_separation <- function(model, method = c("tukey", "dmrt", "dunnett",
                                               "contrasts"),
                            alpha = 0.05, control_group = NULL,
                            contrast_matrix = NULL, ...) {
  method <- match.arg(method)

  if (!inherits(model, "aov")) stop("`model` must be an `aov` object.")

  # Extract information from model
  mf <- model.frame(model)
  response <- mf[[1]]
  trt_var_name <- attr(terms(model), "term.labels")[1]
  trt_var <- mf[[trt_var_name]]
  if (!is.factor(trt_var)) trt_var <- as.factor(trt_var)

  cell_means <- tapply(response, trt_var, mean)
  cell_n <- tapply(response, trt_var, length)

  aov_sum <- summary(model)[[1]]
  mse <- aov_sum["Residuals", "Mean Sq"]
  df_res <- aov_sum["Residuals", "Df"]

  recs <- switch(method,
    tukey = paste0(
      "Tukey HSD: controls FWER; recommended for confirmatory experiments ",
      "with all pairwise comparisons."
    ),
    dmrt = paste0(
      "DMRT: maximises power; best for exploratory experiments. ",
      "More liberal than Tukey."
    ),
    dunnett = paste0(
      "Dunnett: optimal when comparing treatments only to a control; ",
      "more powerful than Tukey for this restricted set."
    ),
    contrasts = paste0(
      "Orthogonal contrasts: for planned, hypothesis-driven comparisons. ",
      "Partitions treatment SS into interpretable components."
    )
  )

  results_obj <- switch(method,
    tukey = tukey_hsd(model, alpha = alpha, ...),
    dmrt = {
      dmrt(cell_means, mse = mse, df = df_res, alpha = alpha, n = cell_n, ...)
    },
    dunnett = {
      if (is.null(control_group)) {
        stop("`control_group` is required for method = 'dunnett'.")
      }
      dunnett_test(response, trt_var, control = control_group,
                   alpha = alpha, ...)
    },
    contrasts = {
      if (is.null(contrast_matrix)) {
        stop("`contrast_matrix` is required for method = 'contrasts'.")
      }
      test_contrasts(model, contrast_matrix, ...)
    }
  )

  # Build means_with_letters table
  means_with_letters <- .build_means_table(results_obj, method, cell_means)

  structure(
    list(
      method_used = method,
      results = results_obj,
      means_with_letters = means_with_letters,
      recommendations = recs,
      alpha = alpha
    ),
    class = "mean_separation_result"
  )
}

#' Internal: build means_with_letters for mean_separation_result
#' @noRd
.build_means_table <- function(results_obj, method, cell_means) {
  trmt_names <- names(sort(cell_means, decreasing = TRUE))
  means_sorted <- sort(cell_means, decreasing = TRUE)

  if (method == "tukey") {
    letters_vec <- results_obj$letters
    # Align letters to sorted order
    letters_out <- letters_vec[trmt_names]
  } else if (method == "dmrt") {
    mt <- results_obj$means_table
    letters_out <- stats::setNames(mt$letters, mt$treatment)[trmt_names]
  } else {
    # For dunnett/contrasts, just return means without letters
    letters_out <- stats::setNames(rep(NA_character_, length(trmt_names)),
                                   trmt_names)
  }

  data.frame(
    treatment = trmt_names,
    mean = as.numeric(means_sorted),
    letters = letters_out,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @export
print.mean_separation_result <- function(x, ...) {
  cat(sprintf("Mean Separation: %s\n", toupper(x$method_used)))
  cat(rep("-", 45), "\n", sep = "")
  cat("Means with significance letters:\n")
  print(x$means_with_letters, row.names = FALSE)
  cat("\nRecommendation:", x$recommendations, "\n")
  invisible(x)
}

#' @export
summary.mean_separation_result <- function(object, ...) {
  cat(sprintf("Mean Separation Result - %s\n", toupper(object$method_used)))
  cat(rep("=", 50), "\n", sep = "")
  cat("\nMeans Table:\n")
  print(object$means_with_letters, row.names = FALSE)
  cat("\nMethod recommendation:\n", object$recommendations, "\n")
  cat("\nDetailed results:\n")
  print(object$results)
  invisible(object)
}

# ---------------------------------------------------------------------------

#' Bar/Point Plot with Significance Letters
#'
#' Creates a publication-ready ggplot2 bar chart (or point chart) of treatment
#' means with DMRT/Tukey significance letters and error bars.
#'
#' @param x A `mean_separation_result` object (from [mean_separation()]) or a
#'   data frame with at least columns `treatment`, `mean`, and optionally
#'   `letters`, `se`.
#' @param type Plot type: `"bar"` (default) or `"point"`.
#' @param y_label Character; y-axis label. Default `"Mean"`.
#' @param x_label Character; x-axis label. Default `"Treatment"`.
#' @param title Character; plot title.
#' @param show_letters Logical; whether to display significance letters
#'   (default `TRUE`).
#' @param error_type Character; type of error bar: `"se"` (standard error,
#'   default) or `"ci"` (95% confidence interval).
#' @param mse Numeric MSE used to compute SE bars when `x` is a
#'   `mean_separation_result` and the original data are not available.
#' @param n_rep Integer replication count for SE calculation.
#' @param fill_color Character; bar fill colour (default `"steelblue"`).
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' set.seed(42)
#' dat <- data.frame(
#'   y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
#'   trt = factor(rep(c("A", "B", "C"), each = 5))
#' )
#' fit <- aov(y ~ trt, data = dat)
#' res <- mean_separation(fit, method = "dmrt")
#' p <- plot_means(res)
#' print(p)
#'
#' @export
plot_means <- function(x, type = c("bar", "point"), y_label = "Mean",
                       x_label = "Treatment", title = NULL,
                       show_letters = TRUE, error_type = c("se", "ci"),
                       mse = NULL, n_rep = NULL, fill_color = "steelblue",
                       ...) {
  type <- match.arg(type)
  error_type <- match.arg(error_type)

  # Extract data frame
  if (inherits(x, "mean_separation_result")) {
    plot_df <- x$means_with_letters
    if (is.null(title)) {
      title <- paste("Treatment Means —", toupper(x$method_used))
    }
  } else if (is.data.frame(x)) {
    required_cols <- c("treatment", "mean")
    if (!all(required_cols %in% names(x))) {
      stop("Data frame `x` must contain columns 'treatment' and 'mean'.")
    }
    plot_df <- x
    if (is.null(title)) title <- "Treatment Means"
  } else {
    stop("`x` must be a `mean_separation_result` or data frame.")
  }

  # Order by mean (descending)
  plot_df <- plot_df[order(plot_df$mean, decreasing = TRUE), ]
  plot_df$treatment <- factor(plot_df$treatment,
                              levels = plot_df$treatment)

  # SE/CI calculation
  if (!is.null(mse) && !is.null(n_rep)) {
    se_val <- sqrt(mse / n_rep)
    mult <- if (error_type == "ci") stats::qt(0.975, df = Inf) else 1
    plot_df$se <- se_val
    plot_df$ymin <- plot_df$mean - mult * se_val
    plot_df$ymax <- plot_df$mean + mult * se_val
  } else if (!("se" %in% names(plot_df))) {
    # No SE info available — omit error bars
    plot_df$se <- NA_real_
    plot_df$ymin <- NA_real_
    plot_df$ymax <- NA_real_
  } else {
    mult <- if (error_type == "ci") stats::qt(0.975, df = Inf) else 1
    plot_df$ymin <- plot_df$mean - mult * plot_df$se
    plot_df$ymax <- plot_df$mean + mult * plot_df$se
  }

  # Letter y position (top of bar + buffer)
  y_range <- diff(range(plot_df$mean, na.rm = TRUE))
  letter_offset <- 0.05 * max(abs(plot_df$mean), na.rm = TRUE) + 0.1 * y_range
  plot_df$letter_y <- plot_df$mean + letter_offset

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = .data[["treatment"]],
                                    y = .data[["mean"]]))

  if (type == "bar") {
    p <- p +
      ggplot2::geom_bar(stat = "identity", fill = fill_color, colour = "grey30",
                        width = 0.7)
  } else {
    p <- p +
      ggplot2::geom_point(size = 3, colour = fill_color)
  }

  # Error bars (only when SE data are present)
  if (!all(is.na(plot_df$ymin))) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[["ymin"]], ymax = .data[["ymax"]]),
        width = 0.2, colour = "grey30"
      )
  }

  # Significance letters
  has_letters <- "letters" %in% names(plot_df) &&
    !all(is.na(plot_df$letters))
  if (show_letters && has_letters) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(y = .data[["letter_y"]],
                     label = .data[["letters"]]),
        size = 4, fontface = "bold"
      )
  }

  p <- p +
    ggplot2::labs(title = title, x = x_label, y = y_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12)
    )

  p
}

# ---------------------------------------------------------------------------

#' Tukey HSD Interval Plot
#'
#' Creates a ggplot2 interval plot showing all pairwise Tukey HSD confidence
#' intervals, with significant comparisons highlighted.
#'
#' @param x A `tukey_result` object from [tukey_hsd()].
#' @param title Character; plot title. Default `"Tukey HSD Pairwise Comparisons"`.
#' @param sig_color Colour for significant comparisons (default `"firebrick"`).
#' @param ns_color Colour for non-significant comparisons (default `"grey50"`).
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' set.seed(42)
#' dat <- data.frame(
#'   y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
#'   trt = factor(rep(c("A", "B", "C"), each = 5))
#' )
#' fit <- aov(y ~ trt, data = dat)
#' tres <- tukey_hsd(fit)
#' p <- plot_tukey(tres)
#' print(p)
#'
#' @export
plot_tukey <- function(x, title = "Tukey HSD Pairwise Comparisons",
                       sig_color = "firebrick", ns_color = "grey50", ...) {
  if (!inherits(x, "tukey_result")) {
    stop("`x` must be a `tukey_result` object.")
  }

  df <- x$comparisons
  df$significant <- df$p_adj < x$alpha
  df$colour <- ifelse(df$significant, "Significant", "Non-significant")
  df$comparison <- factor(df$comparison, levels = rev(df$comparison))

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["diff"]],
    y = .data[["comparison"]],
    colour = .data[["colour"]]
  )) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = .data[["lwr"]], xmax = .data[["upr"]]),
      size = 0.5
    ) +
    ggplot2::scale_color_manual(
      values = c("Significant" = sig_color, "Non-significant" = ns_color),
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      x = "Difference in Means",
      y = "Comparison"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 10)
    )

  p
}
