test_that("agristat package loads for post-hoc tests", {
  expect_true(requireNamespace("agristat", quietly = TRUE))
})

# Placeholder tests for Module 4 — Mean Separation & Post-Hoc Tests.
# Full tests will be added in Phase 5 when functions are implemented.

# test_that("tukey_test() returns a valid multiple comparison result", {
#   df <- data.frame(
#     treatment = factor(rep(c("A", "B", "C"), each = 10)),
#     value     = c(rnorm(10, 5), rnorm(10, 7), rnorm(10, 6))
#   )
#   model  <- aov(value ~ treatment, data = df)
#   result <- tukey_test(model)
#   expect_s3_class(result, "posthoc_result")
#   expect_true(!is.null(result$letters))
# })

# test_that("format_mean_table() returns a data frame", {
#   df <- data.frame(
#     treatment = factor(rep(c("A", "B"), each = 5)),
#     value     = rnorm(10, 5)
#   )
#   model  <- aov(value ~ treatment, data = df)
#   result <- tukey_test(model)
#   tbl    <- format_mean_table(result)
#   expect_s3_class(tbl, "data.frame")
# })
