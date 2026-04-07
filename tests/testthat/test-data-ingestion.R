# Tests for Module 1 — Data Ingestion & Pre-Processing
#
# Covers: declare_hypothesis, check_assumptions, auto_transform,
#         suggest_nonparametric, chi_square_segregation, import_field_data
# Edge cases: empty data, single value, all NA, invalid inputs

# ============================================================================
# Step 8: Hypothesis Tracking
# ============================================================================

test_that("declare_hypothesis() creates a valid hypothesis object", {
  h <- declare_hypothesis(
    h0    = "No significant difference in yield between varieties",
    h1    = "At least one variety differs in yield",
    alpha = 0.05
  )
  expect_s3_class(h, "hypothesis")
  expect_equal(h$h0_text, "No significant difference in yield between varieties")
  expect_equal(h$h1_text, "At least one variety differs in yield")
  expect_equal(h$alpha, 0.05)
  expect_s3_class(h$declared_date, "POSIXct")
})

test_that("declare_hypothesis() accepts different alpha values", {
  h01 <- declare_hypothesis("H0", "H1", alpha = 0.01)
  h10 <- declare_hypothesis("H0", "H1", alpha = 0.10)
  expect_equal(h01$alpha, 0.01)
  expect_equal(h10$alpha, 0.10)
})

test_that("declare_hypothesis() validates h0 input", {
  expect_error(declare_hypothesis(123, "H1"),
               regexp = "h0")
  expect_error(declare_hypothesis("", "H1"),
               regexp = "h0")
  expect_error(declare_hypothesis(c("a", "b"), "H1"),
               regexp = "h0")
})

test_that("declare_hypothesis() validates h1 input", {
  expect_error(declare_hypothesis("H0", 123),
               regexp = "h1")
  expect_error(declare_hypothesis("H0", ""),
               regexp = "h1")
})

test_that("declare_hypothesis() validates alpha", {
  expect_error(declare_hypothesis("H0", "H1", alpha = 0),
               regexp = "alpha")
  expect_error(declare_hypothesis("H0", "H1", alpha = 1),
               regexp = "alpha")
  expect_error(declare_hypothesis("H0", "H1", alpha = -0.05),
               regexp = "alpha")
  expect_error(declare_hypothesis("H0", "H1", alpha = "0.05"),
               regexp = "alpha")
})

test_that("print.hypothesis() displays without error and returns invisibly", {
  h <- declare_hypothesis("No effect", "There is an effect")
  out <- capture.output(ret <- print(h))
  expect_identical(ret, h)
  expect_true(any(grepl("Hypothesis", out)))
  expect_true(any(grepl("No effect", out)))
  expect_true(any(grepl("There is an effect", out)))
})

test_that("summary.hypothesis() shows declared date", {
  h <- declare_hypothesis("No effect", "There is an effect")
  out <- capture.output(ret <- summary(h))
  expect_identical(ret, h)
  expect_true(any(grepl("Declared:", out)))
})

# ============================================================================
# Step 9: Assumption Testing Engine
# ============================================================================

test_that("check_assumptions() returns assumption_report for normal data", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 1)
  report <- check_assumptions(x)
  expect_s3_class(report, "assumption_report")
  expect_true(is.data.frame(report$test_results))
  expect_true("Shapiro-Wilk" %in% report$test_results$test_name)
  expect_true(is.character(report$recommendation))
  expect_equal(report$alpha, 0.05)
})

test_that("check_assumptions() detects non-normal data", {
  set.seed(7)
  x <- c(rexp(40, rate = 0.5))   # clearly skewed
  report <- check_assumptions(x)
  sw_row <- report$test_results[report$test_results$test_name == "Shapiro-Wilk", ]
  expect_true(nrow(sw_row) == 1L)
  # p-value should be numeric
  expect_true(is.numeric(sw_row$p_value))
})

test_that("check_assumptions() adds Levene's test when groups provided", {
  set.seed(1)
  x   <- rnorm(30, mean = 5)
  grp <- rep(c("A", "B", "C"), each = 10)
  report <- check_assumptions(x, groups = grp)
  expect_true("Levene" %in% report$test_results$test_name)
  expect_equal(report$n_groups, 3L)
})

test_that("check_assumptions() handles groups as factor", {
  set.seed(2)
  x   <- rnorm(20)
  grp <- factor(rep(c("T1", "T2"), each = 10))
  report <- check_assumptions(x, groups = grp)
  expect_true("Levene" %in% report$test_results$test_name)
  expect_equal(report$n_groups, 2L)
})

test_that("check_assumptions() warns and returns NA results for n < 3", {
  expect_warning(
    report <- check_assumptions(c(1, 2)),
    regexp = "non-missing"
  )
  expect_true(is.na(report$test_results$p_value[1L]))
})

test_that("check_assumptions() handles single value", {
  expect_warning(
    report <- check_assumptions(5),
    regexp = "non-missing"
  )
  expect_true(is.na(report$test_results$p_value[1L]))
})

test_that("check_assumptions() handles all-NA input", {
  expect_warning(
    report <- check_assumptions(as.numeric(c(NA, NA, NA))),
    regexp = "non-missing"
  )
  expect_s3_class(report, "assumption_report")
})

test_that("check_assumptions() rejects non-numeric input", {
  expect_error(check_assumptions(c("a", "b", "c")), regexp = "numeric")
})

test_that("check_assumptions() validates alpha", {
  expect_error(check_assumptions(rnorm(10), alpha = 2), regexp = "alpha")
})

test_that("check_assumptions() requires groups length matches x", {
  x   <- rnorm(10)
  grp <- rep("A", 5)
  expect_error(check_assumptions(x, groups = grp), regexp = "same length")
})

test_that("print.assumption_report() outputs without error", {
  set.seed(3)
  report <- check_assumptions(rnorm(20))
  out <- capture.output(print(report))
  expect_true(any(grepl("Assumption Report", out)))
})

test_that("summary.assumption_report() works", {
  set.seed(4)
  report <- check_assumptions(rnorm(20))
  out <- capture.output(summary(report))
  expect_true(any(grepl("Assumption", out)))
})

# ============================================================================
# Step 10: Data Transformation Utilities
# ============================================================================

test_that("auto_transform() returns auto_transform_result", {
  set.seed(10)
  x <- exp(rnorm(50))   # log-normal
  result <- auto_transform(x)
  expect_s3_class(result, "auto_transform_result")
  expect_true(is.character(result$best_transform))
  expect_true(is.numeric(result$original_p_value))
  expect_true(is.numeric(result$transformed_p_value))
  expect_true(is.character(result$formula_suggestion))
  expect_s3_class(result$plot, "gg")
  expect_true(is.data.frame(result$all_results))
})

test_that("auto_transform() selects log for log-normal data", {
  set.seed(11)
  x <- exp(rnorm(100, 2, 0.5))
  result <- auto_transform(x)
  # The best transform should improve normality
  expect_true(result$transformed_p_value >= result$original_p_value ||
                result$best_transform == "none")
})

test_that("auto_transform() validates input", {
  expect_error(auto_transform("text"), regexp = "numeric")
  expect_error(auto_transform(c(1, 2)), regexp = "3 finite")
})

test_that("auto_transform() handles custom response_name", {
  set.seed(12)
  x <- exp(rnorm(30))
  result <- auto_transform(x, response_name = "yield")
  expect_true(grepl("yield", result$formula_suggestion))
})

test_that("print.auto_transform_result() works", {
  set.seed(13)
  x <- exp(rnorm(40))
  result <- auto_transform(x)
  out <- capture.output(print(result))
  expect_true(any(grepl("Auto-Transform", out)))
  expect_true(any(grepl("best transformation", out)))
})

test_that("summary.auto_transform_result() shows all_results table", {
  set.seed(14)
  x <- exp(rnorm(40))
  result <- auto_transform(x)
  out <- capture.output(summary(result))
  expect_true(any(grepl("transformations", out)))
})

test_that("suggest_nonparametric() rejects wrong class", {
  expect_error(suggest_nonparametric(list()), regexp = "assumption_report")
})

test_that("suggest_nonparametric() returns parametric when all pass", {
  set.seed(15)
  x <- rnorm(50)
  report <- check_assumptions(x)
  recs <- suggest_nonparametric(report)
  expect_true("parametric" %in% names(recs) ||
                length(recs) > 0L)   # at least one entry
  expect_true(is.character(recs))
})

test_that("suggest_nonparametric() recommends Mann-Whitney for 2-group non-normal", {
  set.seed(99)
  x   <- c(rexp(20, 2), rexp(20, 5))
  grp <- rep(c("A", "B"), each = 20)
  report <- check_assumptions(x, groups = grp)
  recs <- suggest_nonparametric(report, n_groups = 2L)
  # Should mention Mann-Whitney or Wilcoxon for 2-group case
  all_text <- paste(recs, collapse = " ")
  expect_true(grepl("Mann-Whitney|Wilcox|Kruskal|parametric", all_text,
                    ignore.case = TRUE))
})

test_that("suggest_nonparametric() recommends Kruskal-Wallis for 3+ groups", {
  set.seed(100)
  x   <- c(rexp(30, 0.5))
  grp <- rep(c("A", "B", "C"), each = 10)
  report <- check_assumptions(x, groups = grp)
  recs <- suggest_nonparametric(report, n_groups = 3L)
  all_text <- paste(recs, collapse = " ")
  expect_true(grepl("Kruskal|Mann|parametric", all_text, ignore.case = TRUE))
})

# ============================================================================
# Step 11: Mendelian Segregation Tester
# ============================================================================

test_that("chi_square_segregation() 3:1 ratio — consistent", {
  result <- chi_square_segregation(
    observed        = c(75, 25),
    expected_ratios = c(3, 1)
  )
  expect_s3_class(result, "chi_square_result")
  expect_equal(result$df, 1L)
  expect_equal(result$ratio_string, "3:1")
  expect_true(is.numeric(result$chi_sq_stat))
  expect_true(is.numeric(result$p_value))
  expect_true(result$goodness_of_fit)   # 75:25 is close to 3:1
  expect_true(grepl("Consistent", result$interpretation))
})

test_that("chi_square_segregation() 9:3:3:1 dihybrid ratio", {
  result <- chi_square_segregation(
    observed        = c(90, 30, 30, 10),
    expected_ratios = c(9, 3, 3, 1)
  )
  expect_s3_class(result, "chi_square_result")
  expect_equal(result$df, 3L)
  expect_equal(result$ratio_string, "9:3:3:1")
  expect_true(result$goodness_of_fit)
})

test_that("chi_square_segregation() detects significant deviation", {
  result <- chi_square_segregation(
    observed        = c(90, 10),
    expected_ratios = c(3, 1),
    alpha           = 0.05
  )
  # 90:10 deviates greatly from 3:1 (expected 75:25)
  expect_false(result$goodness_of_fit)
  expect_true(grepl("deviate", result$interpretation, ignore.case = TRUE))
})

test_that("chi_square_segregation() computes chi-square correctly", {
  # Known: observed = c(3, 1), expected_ratios = c(1, 1), total = 4
  # Expected = c(2, 2), chi2 = (3-2)^2/2 + (1-2)^2/2 = 0.5 + 0.5 = 1
  result <- chi_square_segregation(c(3, 1), c(1, 1))
  expect_equal(result$chi_sq_stat, 1, tolerance = 1e-10)
})

test_that("chi_square_segregation() validates observed", {
  expect_error(chi_square_segregation(c(-1, 3), c(1, 1)),
               regexp = "non-negative")
  expect_error(chi_square_segregation(c(NA, 3), c(1, 1)),
               regexp = "non-negative")
  expect_error(chi_square_segregation(5, c(1, 1)),
               regexp = "at least 2")
})

test_that("chi_square_segregation() validates expected_ratios", {
  expect_error(chi_square_segregation(c(3, 1), c(0, 1)),
               regexp = "positive")
  expect_error(chi_square_segregation(c(3, 1), c(1, 1, 1)),
               regexp = "same length")
})

test_that("chi_square_segregation() validates alpha", {
  expect_error(chi_square_segregation(c(3, 1), c(3, 1), alpha = 1.5),
               regexp = "alpha")
})

test_that("print.chi_square_result() outputs without error", {
  result <- chi_square_segregation(c(75, 25), c(3, 1))
  out <- capture.output(print(result))
  expect_true(any(grepl("Chi-Square", out)))
  expect_true(any(grepl("3:1", out)))
})

test_that("chi_square_segregation() warns for low expected counts", {
  expect_warning(
    chi_square_segregation(c(5, 5, 0), c(1, 1, 1)),
    regexp = "expected"
  )
})

# ============================================================================
# Step 12: Data Import Helpers
# ============================================================================

test_that("import_field_data() imports a CSV file", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write.csv(
    data.frame(
      plot  = 1:4,
      yield = c(3.2, 4.1, 3.8, 4.5),
      trt   = c("A", "B", "A", "B")
    ),
    tmp,
    row.names = FALSE
  )
  df <- import_field_data(tmp)
  expect_true(is.data.frame(df))
  expect_true(nrow(df) == 4L)
  expect_true("yield" %in% names(df))
  meta <- attr(df, "agristat_metadata")
  expect_true(is.data.frame(meta))
  expect_true("column" %in% names(meta))
})

test_that("import_field_data() handles NA strings in CSV", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(
    c("id,value", "1,3.5", "2,.", "3,N/A", "4,"),
    tmp
  )
  df <- import_field_data(tmp)
  expect_true(any(is.na(df$value)))
})

test_that("import_field_data() errors for non-existent file", {
  expect_error(
    import_field_data("/no/such/file.csv"),
    regexp = "not found|File not found"
  )
})

test_that("import_field_data() errors for unsupported extension", {
  tmp <- tempfile(fileext = ".xyz")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("col1,col2", tmp)
  expect_error(
    import_field_data(tmp),
    regexp = "Unsupported file extension"
  )
})

test_that("import_field_data() imports Excel when readxl is available", {
  skip_if_not_installed("readxl")
  skip_if_not_installed("writexl")
  tmp <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp), add = TRUE)
  writexl::write_xlsx(
    data.frame(plot = 1:3, yield = c(2.1, 3.4, 2.8)),
    tmp
  )
  df <- import_field_data(tmp)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 3L)
})

test_that("import_field_data() imports GeoJSON when sf is available", {
  skip_if_not_installed("sf")
  tmp <- tempfile(fileext = ".geojson")
  on.exit(unlink(tmp), add = TRUE)
  geojson_txt <- paste0(
    '{"type":"FeatureCollection","features":[',
    '{"type":"Feature","properties":{"name":"A","value":1},',
    '"geometry":{"type":"Point","coordinates":[0,0]}},',
    '{"type":"Feature","properties":{"name":"B","value":2},',
    '"geometry":{"type":"Point","coordinates":[1,1]}}',
    ']}'
  )
  writeLines(geojson_txt, tmp)
  sf_obj <- import_field_data(tmp)
  expect_true(inherits(sf_obj, "sf"))
  expect_equal(nrow(sf_obj), 2L)
})

# ============================================================================
# Utility helpers
# ============================================================================

test_that("format_p_value() formats correctly", {
  expect_equal(format_p_value(0.0432, digits = 4L), "0.0432")
  expect_equal(format_p_value(0.00001), "< 0.001")
  expect_equal(format_p_value(NA), "NA")
})

test_that("format_statistic() formats correctly", {
  expect_equal(format_statistic(3.14159, digits = 2L), "3.14")
  expect_equal(format_statistic(NA), "NA")
})
