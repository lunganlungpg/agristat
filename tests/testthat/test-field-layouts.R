test_that("agristat package loads for field layout tests", {
  expect_true(requireNamespace("agristat", quietly = TRUE))
})

# Placeholder tests for Module 2 — Field Layout Engine.
# Full tests will be added in Phase 3 when functions are implemented.

# test_that("design_crd() creates a valid CRD layout", {
#   layout <- design_crd(treatments = c("A", "B", "C"), reps = 4)
#   expect_s3_class(layout, "field_design")
#   expect_equal(nrow(layout$plan), 12)
# })

# test_that("design_rcbd() creates a valid RCBD layout", {
#   layout <- design_rcbd(treatments = c("A", "B", "C", "D"), blocks = 3)
#   expect_s3_class(layout, "field_design")
#   expect_equal(layout$n_blocks, 3)
# })

# test_that("plot_field_layout() returns a ggplot object", {
#   layout <- design_crd(treatments = c("A", "B"), reps = 3)
#   p <- plot_field_layout(layout)
#   expect_s3_class(p, "ggplot")
# })
