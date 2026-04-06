#' Assign Significance Letters to Ordered Means
#'
#' Implements the letter-assignment algorithm for mean separation tests (DMRT,
#' Tukey, etc.). Means that are not significantly different share a letter.
#' The algorithm scans through ordered means and assigns the next available
#' letter whenever a significant difference is found.
#'
#' @param sorted_means Named numeric vector of means **sorted in decreasing
#'   order**. Names are the treatment labels.
#' @param critical_ranges Named numeric vector of critical range values
#'   \eqn{R_p} for \eqn{p = 2, 3, \ldots, k}. Names must be `"R2"`, `"R3"`,
#'   …, `"Rk"`. Produced by [dmrt()].
#'
#' @return A character vector of significance letters, same length and names as
#'   `sorted_means`.
#'
#' @details
#' The algorithm follows these steps:
#' 1. Begin with letter `"a"`.
#' 2. For the highest mean, try to extend the current letter group to include
#'    subsequent means. A mean joins the current group if the range
#'    (highest − current) \eqn{<} \eqn{R_p} for the appropriate span \eqn{p}.
#' 3. When a significant difference is found, start a new letter.
#' 4. A treatment may carry multiple letters if it is not significantly
#'    different from treatments in two separate groups.
#'
#' @examples
#' means_sorted <- c(A = 14.8, C = 12.5, D = 11.0, B = 10.2)
#' cr <- c(R2 = 2.1, R3 = 2.3, R4 = 2.4)
#' assign_letters(means_sorted, cr)
#'
#' @export
assign_letters <- function(sorted_means, critical_ranges) {
  k <- length(sorted_means)
  if (k == 0) return(character(0))
  if (k == 1) return(stats::setNames("a", names(sorted_means)))

  trmt_names <- names(sorted_means)
  if (is.null(trmt_names)) trmt_names <- paste0("T", seq_len(k))

  # Each treatment accumulates a set of letter indices
  letter_sets <- vector("list", k)
  for (i in seq_len(k)) letter_sets[[i]] <- integer(0)

  current_letter <- 1L

  # Standard DMRT letter assignment:
  # For each starting position i, scan forward; assign the same letter until
  # the range exceeds the critical range R_p (where p = j - i + 1).
  i <- 1L
  while (i <= k) {
    # Start a new letter group from position i
    letter_sets[[i]] <- c(letter_sets[[i]], current_letter)

    j <- i + 1L
    while (j <= k) {
      p <- j - i + 1L  # span
      r_name <- paste0("R", p)
      r_p <- if (r_name %in% names(critical_ranges)) {
        critical_ranges[[r_name]]
      } else {
        critical_ranges[[length(critical_ranges)]]  # use largest available
      }

      range_ij <- sorted_means[i] - sorted_means[j]

      if (range_ij < r_p) {
        # Not significantly different — include j in current letter group
        letter_sets[[j]] <- c(letter_sets[[j]], current_letter)
        j <- j + 1L
      } else {
        # Significantly different — stop extending this group
        break
      }
    }

    # Move to the first position not yet in any letter group
    # (i.e., find the next i that needs a new letter)
    next_i <- i + 1L
    while (next_i <= k && current_letter %in% letter_sets[[next_i]]) {
      next_i <- next_i + 1L
    }
    if (next_i > k) break

    current_letter <- current_letter + 1L
    i <- next_i
  }

  # Convert letter indices to characters
  all_letters <- c(letters, paste0(letters, letters))  # a..z, aa..zz
  result <- vapply(letter_sets, function(idx) {
    if (length(idx) == 0) return("?")
    paste(all_letters[sort(unique(idx))], collapse = "")
  }, character(1))

  stats::setNames(result, trmt_names)
}

# ---------------------------------------------------------------------------

#' Assign Letters from a Significance Matrix
#'
#' Alternative letter-assignment entry point used internally by Tukey HSD.
#' Given a matrix where `sig_mat[i, j] = TRUE` means treatments `i` and `j`
#' are significantly different, produces a letter vector.
#'
#' @param sorted_names Character vector of treatment names (sorted by mean,
#'   decreasing).
#' @param sig_mat Logical matrix (k x k) of significance indicators, with row
#'   and column names equal to `sorted_names`.
#'
#' @return A named character vector of significance letters.
#'
#' @keywords internal
#' @export
assign_letters_from_matrix <- function(sorted_names, sig_mat) {
  k <- length(sorted_names)
  if (k == 0) return(character(0))
  if (k == 1) return(stats::setNames("a", sorted_names))

  letter_sets <- vector("list", k)
  for (i in seq_len(k)) letter_sets[[i]] <- integer(0)

  current_letter <- 1L
  i <- 1L

  while (i <= k) {
    letter_sets[[i]] <- c(letter_sets[[i]], current_letter)

    j <- i + 1L
    while (j <= k) {
      name_i <- sorted_names[i]
      name_j <- sorted_names[j]
      # Check if any treatment in the current group is significantly different
      # from j (conservative: extend group only if ALL are non-sig)
      group_members <- sorted_names[vapply(letter_sets[seq_len(j - 1)],
                                           function(s) current_letter %in% s,
                                           logical(1))]
      any_sig <- any(sig_mat[group_members, name_j], na.rm = TRUE)

      if (!any_sig) {
        letter_sets[[j]] <- c(letter_sets[[j]], current_letter)
        j <- j + 1L
      } else {
        break
      }
    }

    next_i <- i + 1L
    while (next_i <= k && current_letter %in% letter_sets[[next_i]]) {
      next_i <- next_i + 1L
    }
    if (next_i > k) break

    current_letter <- current_letter + 1L
    i <- next_i
  }

  all_letters <- c(letters, paste0(letters, letters))
  result <- vapply(letter_sets, function(idx) {
    if (length(idx) == 0) return("?")
    paste(all_letters[sort(unique(idx))], collapse = "")
  }, character(1))

  stats::setNames(result, sorted_names)
}

# ---------------------------------------------------------------------------

#' Validate Orthogonality of a Contrast Matrix
#'
#' Checks whether the rows of a contrast matrix satisfy the two conditions for
#' orthogonality: (1) each row sums to zero, and (2) the inner product of any
#' two distinct rows is zero.
#'
#' @param contrast_matrix A numeric matrix with one contrast per row.
#' @param tol Numerical tolerance for zero-checking (default `1e-8`).
#'
#' @return A single logical value: `TRUE` if all conditions are met, `FALSE`
#'   otherwise. Diagnostic messages are printed via `message()`.
#'
#' @examples
#' cmat <- rbind(
#'   c( 1,  1, -1, -1),
#'   c( 1, -1,  0,  0),
#'   c( 0,  0,  1, -1)
#' )
#' validate_contrasts(cmat)
#'
#' # Non-orthogonal example
#' bad <- rbind(c(1, -1, 0), c(1, 0, -1))
#' validate_contrasts(bad)
#'
#' @export
validate_contrasts <- function(contrast_matrix, tol = 1e-8) {
  if (!is.matrix(contrast_matrix)) {
    contrast_matrix <- as.matrix(contrast_matrix)
  }
  nr <- nrow(contrast_matrix)
  all_ok <- TRUE

  # Check row sums = 0
  row_sums <- rowSums(contrast_matrix)
  bad_rows <- which(abs(row_sums) > tol)
  if (length(bad_rows) > 0) {
    message(sprintf(
      "Row(s) %s do not sum to zero (sums: %s).",
      paste(bad_rows, collapse = ", "),
      paste(round(row_sums[bad_rows], 6), collapse = ", ")
    ))
    all_ok <- FALSE
  }

  # Check pairwise inner products = 0
  if (nr > 1) {
    for (i in seq_len(nr - 1)) {
      for (j in seq(i + 1, nr)) {
        ip <- sum(contrast_matrix[i, ] * contrast_matrix[j, ])
        if (abs(ip) > tol) {
          message(sprintf(
            "Rows %d and %d are not orthogonal (inner product = %.6f).", i, j, ip
          ))
          all_ok <- FALSE
        }
      }
    }
  }

  all_ok
}

# ---------------------------------------------------------------------------

#' Format Contrast Output for Display
#'
#' Converts a `contrast_result` object into a neatly formatted character
#' string suitable for inclusion in reports or vignettes.
#'
#' @param contrast_result A `contrast_result` object from [test_contrasts()].
#' @param digits Integer; number of decimal places (default `4`).
#' @param alpha Significance level for marking significance (default `0.05`).
#'
#' @return A character string (invisibly) with a formatted table. The string
#'   is also printed to the console.
#'
#' @examples
#' set.seed(7)
#' dat <- data.frame(
#'   y = c(rnorm(4, 10), rnorm(4, 12), rnorm(4, 14), rnorm(4, 11)),
#'   trt = factor(rep(c("A", "B", "C", "D"), each = 4))
#' )
#' fit <- aov(y ~ trt, data = dat)
#' cmat <- rbind(
#'   "A+B vs C+D" = c( 1,  1, -1, -1),
#'   "A vs B"     = c( 1, -1,  0,  0),
#'   "C vs D"     = c( 0,  0,  1, -1)
#' )
#' res <- test_contrasts(fit, cmat)
#' format_contrast_output(res)
#'
#' @export
format_contrast_output <- function(contrast_result, digits = 4, alpha = 0.05) {
  if (!inherits(contrast_result, "contrast_result")) {
    stop("`contrast_result` must be a `contrast_result` object.")
  }

  df <- contrast_result$contrasts
  df$sig <- ifelse(df$p_value < alpha, "*", " ")

  header <- sprintf("%-20s %10s %4s %10s %10s %4s\n",
                    "Contrast", "SS", "df", "F", "p-value", "Sig")
  sep_line <- paste(rep("-", nchar(header) - 1), collapse = "")

  rows <- apply(df, 1, function(r) {
    sprintf("%-20s %10.4f %4s %10.4f %10.4f %4s",
            r["contrast"],
            as.numeric(r["SS"]),
            r["df"],
            as.numeric(r["F_value"]),
            as.numeric(r["p_value"]),
            r["sig"])
  })

  out <- paste0(
    "\nOrthogonal Contrasts\n",
    sep_line, "\n",
    header,
    sep_line, "\n",
    paste(rows, collapse = "\n"), "\n",
    sep_line, "\n",
    sprintf("Orthogonal: %s\n", contrast_result$orthogonal)
  )

  cat(out)
  invisible(out)
}

# ---------------------------------------------------------------------------
# Internal utility: p-value formatting
# ---------------------------------------------------------------------------

#' Format p-value for display
#' @param p Numeric p-value.
#' @param digits Integer digits.
#' @noRd
.fmt_pval <- function(p, digits = 4) {
  threshold <- 10^(-digits)
  ifelse(p < threshold, paste0("< ", threshold), round(p, digits))
}
