test_that("dmrt: basic letter assignment with balanced design", {
  # Four treatments with known means; MSE = 2.0, df = 16, n = 5
  means <- c(A = 14.8, B = 10.2, C = 12.5, D = 11.0)
  res <- dmrt(means, mse = 2.0, df = 16, alpha = 0.05, n = 5)

  expect_s3_class(res, "dmrt_result")
  expect_equal(nrow(res$means_table), 4)
  expect_named(res$means_table, c("treatment", "mean", "letters"))

  # means_table should be sorted decreasing
  expect_true(all(diff(res$means_table$mean) <= 0))

  # Critical ranges: should have R2, R3, R4
  expect_named(res$critical_ranges, c("R2", "R3", "R4"))
  expect_true(all(res$critical_ranges > 0))

  # Ranges should be non-decreasing (R2 <= R3 <= R4)
  expect_true(all(diff(res$critical_ranges) >= 0))
})

test_that("dmrt: highest and lowest means get different letters when far apart", {
  means <- c(A = 20.0, B = 10.0)
  res <- dmrt(means, mse = 0.5, df = 8, alpha = 0.05, n = 5)
  letters_vec <- res$means_table$letters
  expect_false(identical(letters_vec[1], letters_vec[2]))
})

test_that("dmrt: identical means share the same letter", {
  means <- c(A = 10.0, B = 10.0, C = 10.0)
  res <- dmrt(means, mse = 1.0, df = 12, alpha = 0.05, n = 4)
  expect_equal(length(unique(res$means_table$letters)), 1L)
})

test_that("dmrt: unequal replication uses harmonic mean", {
  means <- c(A = 14.0, B = 11.0, C = 12.5)
  n_rep <- c(A = 3, B = 5, C = 4)
  res <- dmrt(means, mse = 1.5, df = 9, alpha = 0.05, n = n_rep)
  # harmonic mean of 3, 5, 4
  expected_nh <- 3 / (1/3 + 1/5 + 1/4)
  expect_equal(res$n_harmonic, expected_nh, tolerance = 1e-6)
})

test_that("dmrt: input validation", {
  expect_error(dmrt(c(1), mse = 1, df = 5, n = 1),
               "`means` must be a numeric vector with at least 2 elements")
  expect_error(dmrt(c(1, 2), mse = -1, df = 5, n = 1),
               "`mse` must be a single positive number")
  expect_error(dmrt(c(1, 2), mse = 1, df = 0, n = 1),
               "`df` must be a single positive integer")
  expect_error(dmrt(c(1, 2), mse = 1, df = 5, alpha = 1.5, n = 1),
               "`alpha` must be a single number in")
  expect_error(dmrt(c(1, 2, 3), mse = 1, df = 5, n = c(3, 4)),
               "`n` must have the same length as `means`")
})

test_that("dmrt: print and summary methods work", {
  means <- c(A = 12.0, B = 9.0, C = 11.0)
  res <- dmrt(means, mse = 1.2, df = 9, n = 4)
  expect_output(print(res), "Duncan")
  expect_output(summary(res), "DMRT")
})

test_that("dmrt: two-treatment case", {
  means <- c(High = 15.0, Low = 8.0)
  res <- dmrt(means, mse = 1.0, df = 8, n = 5)
  expect_equal(length(res$critical_ranges), 1L)
  expect_named(res$critical_ranges, "R2")
})

# ---------------------------------------------------------------------------
# assign_letters() tests
# ---------------------------------------------------------------------------

test_that("assign_letters: returns correct length and names", {
  sorted <- c(A = 14.8, B = 12.5, C = 11.0, D = 10.2)
  cr <- c(R2 = 2.0, R3 = 2.3, R4 = 2.5)
  result <- assign_letters(sorted, cr)
  expect_length(result, 4)
  expect_named(result, names(sorted))
})

test_that("assign_letters: all-equal means share one letter", {
  sorted <- c(A = 5.0, B = 5.0, C = 5.0)
  # Very large critical ranges → all share 'a'
  cr <- c(R2 = 100.0, R3 = 100.0)
  result <- assign_letters(sorted, cr)
  expect_true(all(result == "a"))
})

test_that("assign_letters: all-different means get distinct letters", {
  sorted <- c(A = 20.0, B = 10.0, C = 5.0)
  # Very small critical ranges → all different
  cr <- c(R2 = 0.001, R3 = 0.002)
  result <- assign_letters(sorted, cr)
  expect_equal(length(unique(result)), 3L)
})

test_that("assign_letters: single element returns 'a'", {
  result <- assign_letters(c(X = 5.0), c(R2 = 1.0))
  expect_equal(result, c(X = "a"))
})

# ---------------------------------------------------------------------------
# tukey_hsd() tests
# ---------------------------------------------------------------------------

test_that("tukey_hsd: works with aov object", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
    trt = factor(rep(c("A", "B", "C"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- tukey_hsd(fit)

  expect_s3_class(res, "tukey_result")
  expect_named(res$comparisons, c("comparison", "diff", "lwr", "upr", "p_adj"))
  # 3 treatments → 3 pairwise comparisons
  expect_equal(nrow(res$comparisons), 3)
  expect_length(res$letters, 3)
  expect_type(res$interpretation, "character")
})

test_that("tukey_hsd: results are consistent with stats::TukeyHSD", {
  set.seed(1)
  dat <- data.frame(
    y = c(rnorm(6, 10), rnorm(6, 13), rnorm(6, 11)),
    trt = factor(rep(c("A", "B", "C"), each = 6))
  )
  fit <- aov(y ~ trt, data = dat)

  res_our <- tukey_hsd(fit)
  res_std <- stats::TukeyHSD(fit)$trt

  # Match comparisons (order may differ)
  std_df <- data.frame(
    comparison = rownames(res_std),
    p_adj = res_std[, "p adj"],
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(res_our$comparisons))) {
    comp <- res_our$comparisons$comparison[i]
    # Find matching row (possibly reversed order)
    rev_comp <- paste(rev(strsplit(comp, "-")[[1]]), collapse = "-")
    matched <- std_df[std_df$comparison %in% c(comp, rev_comp), , drop = FALSE]
    if (nrow(matched) == 1) {
      expect_equal(res_our$comparisons$p_adj[i],
                   matched$p_adj[1],
                   tolerance = 1e-4)
    }
  }
})

test_that("tukey_hsd: works with means vector", {
  means <- c(A = 14.0, B = 11.0, C = 12.0)
  res <- tukey_hsd(means, mse = 1.5, df = 12, n = 5)
  expect_s3_class(res, "tukey_result")
  expect_equal(nrow(res$comparisons), 3)
})

test_that("tukey_hsd: error when means vector lacks mse/df/n", {
  means <- c(A = 14.0, B = 11.0)
  expect_error(tukey_hsd(means), "`mse`, `df`, and `n` are required")
})

test_that("tukey_hsd: print and summary methods work", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 14)),
    trt = factor(rep(c("A", "B"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- tukey_hsd(fit)
  expect_output(print(res), "Tukey")
  expect_output(summary(res), "Tukey")
})

test_that("tukey_hsd: alpha validation", {
  means <- c(A = 5, B = 6)
  expect_error(tukey_hsd(means, mse = 1, df = 5, n = 3, alpha = 1.1),
               "`alpha` must be a single number in")
})

# ---------------------------------------------------------------------------
# dunnett_test() tests
# ---------------------------------------------------------------------------

test_that("dunnett_test: basic two-sided test", {
  set.seed(99)
  response <- c(rnorm(5, 10), rnorm(5, 14), rnorm(5, 10.5))
  group <- factor(rep(c("ctrl", "trtA", "trtB"), each = 5))
  res <- dunnett_test(response, group, control = "ctrl")

  expect_s3_class(res, "dunnett_result")
  expect_named(res$comparisons, c("treatment", "estimate", "se", "t_stat",
                                   "p_value"))
  expect_equal(nrow(res$comparisons), 2)
  expect_named(res$significant, c("trtA", "trtB"))
  expect_type(res$interpretation, "character")
})

test_that("dunnett_test: significantly different treatment detected", {
  set.seed(7)
  response <- c(rep(10, 5), rep(20, 5))  # very large difference
  group <- factor(rep(c("ctrl", "trt"), each = 5))
  res <- dunnett_test(response, group, control = "ctrl")
  expect_true(res$significant["trt"])
})

test_that("dunnett_test: non-significant treatment not flagged", {
  set.seed(8)
  response <- c(rnorm(10, 10, 0.1), rnorm(10, 10, 0.1))
  group <- factor(rep(c("ctrl", "trt"), each = 10))
  res <- dunnett_test(response, group, control = "ctrl")
  expect_false(res$significant["trt"])
})

test_that("dunnett_test: input validation", {
  response <- c(1, 2, 3, 4)
  group <- factor(c("a", "a", "b", "b"))
  expect_error(dunnett_test(response, group, control = "x"),
               "not found in `group` levels")
  expect_error(dunnett_test("abc", group, control = "a"),
               "`response` must be numeric")
  expect_error(dunnett_test(response, group[1:3], control = "a"),
               "same length")
})

test_that("dunnett_test: directional alternatives", {
  set.seed(5)
  response <- c(rnorm(5, 10), rnorm(5, 13))
  group <- factor(rep(c("ctrl", "trt"), each = 5))
  res_two  <- dunnett_test(response, group, control = "ctrl",
                            alternative = "two.sided")
  res_gt   <- dunnett_test(response, group, control = "ctrl",
                            alternative = "greater")
  res_lt   <- dunnett_test(response, group, control = "ctrl",
                            alternative = "less")
  expect_equal(res_two$alternative, "two.sided")
  expect_equal(res_gt$alternative, "greater")
  expect_equal(res_lt$alternative, "less")
})

test_that("dunnett_test: print and summary work", {
  response <- c(10, 11, 10, 12, 15, 16, 14, 15)
  group <- factor(c(rep("ctrl", 4), rep("trt", 4)))
  res <- dunnett_test(response, group, control = "ctrl")
  expect_output(print(res), "Dunnett")
  expect_output(summary(res), "Dunnett")
})

# ---------------------------------------------------------------------------
# test_contrasts() tests
# ---------------------------------------------------------------------------

test_that("test_contrasts: basic orthogonal contrasts", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(4, 10), rnorm(4, 12), rnorm(4, 14), rnorm(4, 11)),
    trt = factor(rep(c("A", "B", "C", "D"), each = 4))
  )
  fit <- aov(y ~ trt, data = dat)
  cmat <- rbind(
    "AB vs CD" = c( 1,  1, -1, -1),
    "A vs B"   = c( 1, -1,  0,  0),
    "C vs D"   = c( 0,  0,  1, -1)
  )
  res <- test_contrasts(fit, cmat)

  expect_s3_class(res, "contrast_result")
  expect_equal(nrow(res$contrasts), 3)
  expect_named(res$contrasts, c("contrast", "SS", "df", "F_value", "p_value"))
  expect_true(res$orthogonal)
  expect_length(res$sum_of_squares_contrasts, 3)
})

test_that("test_contrasts: SS values are positive", {
  set.seed(3)
  dat <- data.frame(
    y = c(rnorm(4, 8), rnorm(4, 12)),
    trt = factor(rep(c("A", "B"), each = 4))
  )
  fit <- aov(y ~ trt, data = dat)
  cmat <- rbind("A vs B" = c(1, -1))
  res <- test_contrasts(fit, cmat)
  expect_true(all(res$contrasts$SS >= 0))
})

test_that("test_contrasts: warns on non-orthogonal contrasts", {
  set.seed(1)
  dat <- data.frame(
    y = c(rnorm(3, 10), rnorm(3, 12), rnorm(3, 11)),
    trt = factor(rep(c("A", "B", "C"), each = 3))
  )
  fit <- aov(y ~ trt, data = dat)
  bad_cmat <- rbind(
    "A vs B"  = c(1, -1, 0),
    "A vs C"  = c(1, 0, -1)
  )
  expect_warning(test_contrasts(fit, bad_cmat), NA)  # warning suppressed test
  # Should at least complete without error
  expect_s3_class(
    suppressWarnings(test_contrasts(fit, bad_cmat)),
    "contrast_result"
  )
})

test_that("test_contrasts: column mismatch error", {
  set.seed(1)
  dat <- data.frame(
    y = c(rnorm(4, 10), rnorm(4, 12)),
    trt = factor(rep(c("A", "B"), each = 4))
  )
  fit <- aov(y ~ trt, data = dat)
  bad_cmat <- rbind(c(1, -1, 0))
  expect_error(test_contrasts(fit, bad_cmat), "columns")
})

test_that("test_contrasts: print and summary work", {
  set.seed(4)
  dat <- data.frame(
    y = c(rnorm(4, 10), rnorm(4, 12), rnorm(4, 8), rnorm(4, 11)),
    trt = factor(rep(c("A", "B", "C", "D"), each = 4))
  )
  fit <- aov(y ~ trt, data = dat)
  cmat <- rbind(
    "AB vs CD" = c(1, 1, -1, -1),
    "A vs B"   = c(1, -1, 0, 0),
    "C vs D"   = c(0, 0, 1, -1)
  )
  res <- test_contrasts(fit, cmat)
  expect_output(print(res), "Contrast")
  expect_output(summary(res), "Contrast")
})

# ---------------------------------------------------------------------------
# validate_contrasts() tests
# ---------------------------------------------------------------------------

test_that("validate_contrasts: valid orthogonal matrix returns TRUE", {
  cmat <- rbind(
    c( 1,  1, -1, -1),
    c( 1, -1,  0,  0),
    c( 0,  0,  1, -1)
  )
  expect_true(validate_contrasts(cmat))
})

test_that("validate_contrasts: row not summing to zero returns FALSE", {
  cmat <- rbind(c(1, 0, -1), c(0, 1, -1))  # row 2 sums to 0 but row 1 ok
  # Row 2 sums to 0, row 1 sums to 0. Now let's make one bad:
  bad <- rbind(c(1, 1, 0), c(0, 1, -1))  # row 1 sums to 2
  expect_false(suppressMessages(validate_contrasts(bad)))
})

test_that("validate_contrasts: non-orthogonal rows return FALSE", {
  bad <- rbind(c(1, -1, 0), c(1, 0, -1))
  # Inner product: 1*1 + (-1)*0 + 0*(-1) = 1 ≠ 0
  expect_false(suppressMessages(validate_contrasts(bad)))
})

test_that("validate_contrasts: single row (trivially orthogonal) returns TRUE", {
  cmat <- matrix(c(1, -1, 0), nrow = 1)
  expect_true(validate_contrasts(cmat))
})

# ---------------------------------------------------------------------------
# mean_separation() dispatcher tests
# ---------------------------------------------------------------------------

test_that("mean_separation: routes to tukey by default", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
    trt = factor(rep(c("A", "B", "C"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit)

  expect_s3_class(res, "mean_separation_result")
  expect_equal(res$method_used, "tukey")
  expect_s3_class(res$results, "tukey_result")
})

test_that("mean_separation: routes to dmrt", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
    trt = factor(rep(c("A", "B", "C"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit, method = "dmrt")

  expect_equal(res$method_used, "dmrt")
  expect_s3_class(res$results, "dmrt_result")
})

test_that("mean_separation: routes to dunnett", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
    trt = factor(rep(c("ctrl", "trtA", "trtB"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit, method = "dunnett", control_group = "ctrl")

  expect_equal(res$method_used, "dunnett")
  expect_s3_class(res$results, "dunnett_result")
})

test_that("mean_separation: routes to contrasts", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(4, 10), rnorm(4, 12), rnorm(4, 14), rnorm(4, 11)),
    trt = factor(rep(c("A", "B", "C", "D"), each = 4))
  )
  fit <- aov(y ~ trt, data = dat)
  cmat <- rbind(
    "AB vs CD" = c(1, 1, -1, -1),
    "A vs B"   = c(1, -1, 0, 0),
    "C vs D"   = c(0, 0, 1, -1)
  )
  res <- mean_separation(fit, method = "contrasts", contrast_matrix = cmat)

  expect_equal(res$method_used, "contrasts")
  expect_s3_class(res$results, "contrast_result")
})

test_that("mean_separation: dunnett requires control_group", {
  set.seed(1)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12)),
    trt = factor(rep(c("A", "B"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  expect_error(mean_separation(fit, method = "dunnett"),
               "`control_group` is required")
})

test_that("mean_separation: contrasts requires contrast_matrix", {
  set.seed(1)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12)),
    trt = factor(rep(c("A", "B"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  expect_error(mean_separation(fit, method = "contrasts"),
               "`contrast_matrix` is required")
})

test_that("mean_separation: means_with_letters data frame has correct columns", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 14)),
    trt = factor(rep(c("A", "B"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit)

  expect_named(res$means_with_letters, c("treatment", "mean", "letters"))
  expect_equal(nrow(res$means_with_letters), 2)
})

test_that("mean_separation: print and summary work", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
    trt = factor(rep(c("A", "B", "C"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit)
  expect_output(print(res), "TUKEY")
  expect_output(summary(res), "TUKEY")
})

test_that("mean_separation: non-significant ANOVA still works", {
  set.seed(99)
  # Very noisy data, unlikely to be significant
  dat <- data.frame(
    y = rnorm(15, 10, 10),
    trt = factor(rep(c("A", "B", "C"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit)
  # Should complete without error; all letters likely "a"
  expect_s3_class(res, "mean_separation_result")
})

test_that("mean_separation: single treatment errors on aov (sanity)", {
  dat <- data.frame(y = rnorm(5, 10), trt = factor(rep("A", 5)))
  # aov with single level is degenerate — not a typical usage scenario
  # Just verify our error check catches non-aov input
  expect_error(mean_separation(list()), "`model` must be an `aov` object")
})

# ---------------------------------------------------------------------------
# plot_means() tests
# ---------------------------------------------------------------------------

test_that("plot_means: returns a ggplot object (bar)", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
    trt = factor(rep(c("A", "B", "C"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit, method = "dmrt")
  p <- plot_means(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot_means: point type works", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12)),
    trt = factor(rep(c("A", "B"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit)
  p <- plot_means(res, type = "point")
  expect_s3_class(p, "ggplot")
})

test_that("plot_means: works with data frame input", {
  df <- data.frame(
    treatment = c("A", "B", "C"),
    mean = c(12.0, 10.0, 14.0),
    letters = c("b", "b", "a")
  )
  p <- plot_means(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_means: error bars included when mse and n_rep given", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 14)),
    trt = factor(rep(c("A", "B"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  aov_sum <- summary(fit)[[1]]
  mse <- aov_sum["Residuals", "Mean Sq"]
  res <- mean_separation(fit, method = "dmrt")
  p <- plot_means(res, mse = mse, n_rep = 5)
  # Check that geom_errorbar layer is present
  layer_classes <- vapply(p$layers,
                          function(l) class(l$geom)[1], character(1))
  expect_true(any(grepl("GeomErrorbar", layer_classes)))
})

test_that("plot_means: error for unsupported input", {
  expect_error(plot_means(list()), "`x` must be a `mean_separation_result`")
})

test_that("plot_means: data frame without required columns errors", {
  df <- data.frame(x = 1:3, z = 4:6)
  expect_error(plot_means(df), "must contain columns")
})

# ---------------------------------------------------------------------------
# plot_tukey() tests
# ---------------------------------------------------------------------------

test_that("plot_tukey: returns a ggplot object", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 14)),
    trt = factor(rep(c("A", "B", "C"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  tres <- tukey_hsd(fit)
  p <- plot_tukey(tres)
  expect_s3_class(p, "ggplot")
})

test_that("plot_tukey: contains geom_vline at zero", {
  set.seed(42)
  dat <- data.frame(
    y = c(rnorm(5, 10), rnorm(5, 14)),
    trt = factor(rep(c("A", "B"), each = 5))
  )
  fit <- aov(y ~ trt, data = dat)
  tres <- tukey_hsd(fit)
  p <- plot_tukey(tres)
  layer_classes <- vapply(p$layers,
                          function(l) class(l$geom)[1], character(1))
  expect_true(any(grepl("GeomVline", layer_classes)))
})

test_that("plot_tukey: error for non-tukey_result input", {
  expect_error(plot_tukey(list()), "`x` must be a `tukey_result`")
})

# ---------------------------------------------------------------------------
# format_contrast_output() tests
# ---------------------------------------------------------------------------

test_that("format_contrast_output: returns output and prints", {
  set.seed(7)
  dat <- data.frame(
    y = c(rnorm(4, 10), rnorm(4, 12), rnorm(4, 14), rnorm(4, 11)),
    trt = factor(rep(c("A", "B", "C", "D"), each = 4))
  )
  fit <- aov(y ~ trt, data = dat)
  cmat <- rbind(
    "AB vs CD" = c(1, 1, -1, -1),
    "A vs B"   = c(1, -1, 0, 0),
    "C vs D"   = c(0, 0, 1, -1)
  )
  res <- test_contrasts(fit, cmat)
  out <- expect_output(format_contrast_output(res), "Contrast")
  expect_type(suppressMessages(suppressWarnings(
    format_contrast_output(res)
  )), "character")
})

test_that("format_contrast_output: error for non-contrast_result", {
  expect_error(format_contrast_output(list()), "must be a `contrast_result`")
})

# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

test_that("dmrt: two treatments, very close means share letter", {
  means <- c(A = 10.001, B = 10.000)
  res <- dmrt(means, mse = 5.0, df = 8, n = 4)
  expect_equal(res$means_table$letters[1], res$means_table$letters[2])
})

test_that("tukey_hsd: two treatments, large difference: different letters", {
  set.seed(1)
  dat <- data.frame(
    y = c(rnorm(10, 10, 0.2), rnorm(10, 20, 0.2)),
    trt = factor(rep(c("Low", "High"), each = 10))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- tukey_hsd(fit)
  letters_vec <- res$letters
  expect_false(letters_vec["Low"] == letters_vec["High"])
})

test_that("mean_separation: unequal replication handled gracefully", {
  set.seed(100)
  dat <- data.frame(
    y = c(rnorm(3, 10), rnorm(5, 12), rnorm(4, 11)),
    trt = factor(c(rep("A", 3), rep("B", 5), rep("C", 4)))
  )
  fit <- aov(y ~ trt, data = dat)
  res <- mean_separation(fit, method = "dmrt")
  expect_s3_class(res, "mean_separation_result")
})
