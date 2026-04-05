test_that("agristat package loads", {
  expect_true(requireNamespace("agristat", quietly = TRUE))
})

# Placeholder tests for Module 1 — Data Ingestion & Pre-Processing.
# Full tests will be added in Phase 2 when functions are implemented.

# test_that("declare_hypothesis() creates a hypothesis object", {
#   h <- declare_hypothesis(
#     H0 = "No significant difference between treatments",
#     H1 = "At least one treatment differs",
#     alpha = 0.05
#   )
#   expect_s3_class(h, "hypothesis")
#   expect_equal(h$alpha, 0.05)
# })

# test_that("check_assumptions() returns assumption results", {
#   df <- data.frame(
#     treatment = rep(c("A", "B", "C"), each = 10),
#     value     = rnorm(30, mean = 5)
#   )
#   result <- check_assumptions(df, response = "value", group = "treatment")
#   expect_s3_class(result, "assumption_result")
#   expect_true(!is.null(result$normality))
#   expect_true(!is.null(result$homogeneity))
# })

# test_that("chi_square_segregation() works for 3:1 ratio", {
#   result <- chi_square_segregation(
#     observed = c(75, 25),
#     expected_ratio = c(3, 1)
#   )
#   expect_s3_class(result, "htest")
# })
