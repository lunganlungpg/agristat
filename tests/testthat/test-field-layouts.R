# Tests for Module 2 — Field Layout Engine
#
# Full tests will be added in Phase 3 when functions are implemented.

test_that("agristat package loads successfully", {
  expect_true(requireNamespace("agristat", quietly = TRUE))
})
