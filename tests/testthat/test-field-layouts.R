# Tests for Module 2 — Field Layout Engine
#
# Full tests will be added in Phase 3 when functions are implemented.

test_that("agristat package loads successfully", {
  expect_true(requireNamespace("agristat", quietly = TRUE))
})

# ============================================================================
# Step 16 Tests: CRD
# ============================================================================

test_that("design_crd() returns a field_design object", {
  d <- design_crd(c("A", "B", "C"), replicates = 3L, seed = 1L)
  expect_s3_class(d, "field_design")
  expect_s3_class(d, "crd_design")
})

test_that("design_crd() layout has correct dimensions", {
  d <- design_crd(c("A", "B", "C"), replicates = 4L, seed = 42L)
  expect_equal(nrow(d$layout), 12L)
  expect_named(d$layout, c("plot_id", "treatment", "block", "row", "column"))
})

test_that("design_crd() block column is NA", {
  d <- design_crd(c("T1", "T2"), replicates = 3L, seed = 7L)
  expect_true(all(is.na(d$layout$block)))
})

test_that("design_crd() each treatment appears replicates times", {
  d <- design_crd(c("A", "B", "C"), replicates = 5L, seed = 1L)
  tbl <- table(d$layout$treatment)
  expect_true(all(tbl == 5L))
})

test_that("design_crd() is reproducible with seed", {
  d1 <- design_crd(c("A", "B", "C"), replicates = 4L, seed = 99L)
  d2 <- design_crd(c("A", "B", "C"), replicates = 4L, seed = 99L)
  expect_identical(d1$layout$treatment, d2$layout$treatment)
})

test_that("design_crd() with different seeds produces different layouts", {
  d1 <- design_crd(c("A", "B", "C"), replicates = 6L, seed = 1L)
  d2 <- design_crd(c("A", "B", "C"), replicates = 6L, seed = 2L)
  # Very unlikely to be identical by chance with 18 plots
  expect_false(identical(d1$layout$treatment, d2$layout$treatment))
})

test_that("design_crd() errors on single treatment", {
  expect_error(design_crd(c("A"), replicates = 3L))
})

test_that("design_crd() errors on non-positive replicates", {
  expect_error(design_crd(c("A", "B"), replicates = 0L))
})

test_that("design_crd() plot_id is sequential from 1", {
  d <- design_crd(c("A", "B"), replicates = 3L, seed = 1L)
  expect_equal(d$layout$plot_id, 1L:6L)
})

test_that("design_crd() stores metadata correctly", {
  d <- design_crd(c("A", "B", "C"), replicates = 4L, seed = 1L)
  expect_equal(d$replications, 4L)
  expect_equal(d$randomization_seed, 1L)
  expect_equal(d$design_type, "CRD")
})

# ============================================================================
# Step 16 Tests: RCBD
# ============================================================================

test_that("design_rcbd() returns a field_design object", {
  d <- design_rcbd(c("A", "B", "C"), blocks = 3L, seed = 1L)
  expect_s3_class(d, "field_design")
  expect_s3_class(d, "rcbd_design")
})

test_that("design_rcbd() block structure is balanced", {
  d <- design_rcbd(c("A", "B", "C", "D"), blocks = 4L, seed = 7L)
  # Each block must contain each treatment exactly once
  for (blk in unique(d$layout$block)) {
    trts_in_block <- d$layout$treatment[d$layout$block == blk]
    expect_equal(sort(trts_in_block), sort(c("A", "B", "C", "D")))
  }
})

test_that("design_rcbd() layout has correct number of plots", {
  d <- design_rcbd(c("V1", "V2", "V3"), blocks = 5L, seed = 2L)
  expect_equal(nrow(d$layout), 15L)
})

test_that("design_rcbd() each treatment appears in each block", {
  d <- design_rcbd(c("A", "B"), blocks = 6L, seed = 3L)
  tbl <- table(d$layout$treatment, d$layout$block)
  expect_true(all(tbl == 1L))
})

test_that("design_rcbd() is reproducible with seed", {
  d1 <- design_rcbd(c("A", "B", "C"), blocks = 3L, seed = 42L)
  d2 <- design_rcbd(c("A", "B", "C"), blocks = 3L, seed = 42L)
  expect_identical(d1$layout$treatment, d2$layout$treatment)
})

test_that("design_rcbd() block column has correct values", {
  d <- design_rcbd(c("A", "B"), blocks = 3L, seed = 1L)
  expect_equal(sort(unique(d$layout$block)), c("1", "2", "3"))
})

test_that("design_rcbd() errors on zero blocks", {
  expect_error(design_rcbd(c("A", "B"), blocks = 0L))
})

# ============================================================================
# Step 17 Tests: Factorial
# ============================================================================

test_that("design_factorial() returns a field_design object", {
  d <- design_factorial(
    factor_names  = c("N", "V"),
    factor_levels = list(N = c("Low", "High"), V = c("V1", "V2", "V3")),
    replicates    = 2L, seed = 1L
  )
  expect_s3_class(d, "field_design")
  expect_s3_class(d, "factorial_design")
})

test_that("design_factorial() generates complete grid", {
  d <- design_factorial(
    factor_names  = c("A", "B"),
    factor_levels = list(A = c("a1", "a2"), B = c("b1", "b2", "b3")),
    replicates    = 3L, seed = 1L
  )
  # 2 × 3 × 3 replicates = 18 plots
  expect_equal(nrow(d$layout), 18L)
})

test_that("design_factorial() has correct column names", {
  d <- design_factorial(
    factor_names  = c("Nitrogen", "Variety"),
    factor_levels = list(Nitrogen = c("Low", "High"),
                         Variety  = c("V1", "V2")),
    replicates = 2L, seed = 1L
  )
  expect_true(all(c("Nitrogen", "Variety", "treatment",
                    "plot_id", "row", "column") %in% names(d$layout)))
})

test_that("design_factorial() each cell has correct number of replicates", {
  d <- design_factorial(
    factor_names  = c("A", "B"),
    factor_levels = list(A = c("a1", "a2"), B = c("b1", "b2")),
    replicates    = 4L, seed = 5L
  )
  # 4 cells × 4 reps = 16
  expect_equal(nrow(d$layout), 16L)
  tbl <- table(d$layout$treatment)
  expect_true(all(tbl == 4L))
})

test_that("design_factorial() contrast_coding includes interactions", {
  d <- design_factorial(
    factor_names  = c("A", "B"),
    factor_levels = list(A = c("a1", "a2"), B = c("b1", "b2")),
    replicates    = 2L, seed = 1L
  )
  cc <- d$metadata$contrast_coding
  expect_true(any(grepl("\\*", cc)))  # interaction terms present
})

test_that("design_factorial() 3-way design has correct n_plots", {
  d <- design_factorial(
    factor_names  = c("A", "B", "C"),
    factor_levels = list(A = c("a1", "a2"),
                         B = c("b1", "b2"),
                         C = c("c1", "c2", "c3")),
    replicates    = 2L, seed = 1L
  )
  # 2 × 2 × 3 × 2 = 24
  expect_equal(nrow(d$layout), 24L)
})

test_that("design_factorial() errors with mismatched factor_names", {
  expect_error(
    design_factorial(
      factor_names  = c("A", "B"),
      factor_levels = list(A = c("a1", "a2")),  # missing B
      replicates    = 2L
    )
  )
})

# ============================================================================
# Step 18 Tests: Split-Plot
# ============================================================================

test_that("design_split_plot() returns correct class", {
  d <- design_split_plot(c("I+", "I-"), c("V1", "V2", "V3"),
                         blocks = 3L, seed = 1L)
  expect_s3_class(d, "field_design")
  expect_s3_class(d, "split_plot_design")
})

test_that("design_split_plot() has correct number of plots", {
  d <- design_split_plot(c("WP1", "WP2"), c("SP1", "SP2", "SP3"),
                         blocks = 4L, seed = 2L)
  # 2 × 3 × 4 blocks = 24
  expect_equal(nrow(d$layout), 24L)
})

test_that("design_split_plot() error structure has two strata", {
  d <- design_split_plot(c("I+", "I-"), c("V1", "V2"),
                         blocks = 3L, seed = 1L)
  es <- d$error_structure
  expect_true("whole_plot_error" %in% names(es))
  expect_true("sub_plot_error" %in% names(es))
})

test_that("design_split_plot() whole_plot_level and sub_plot_level present", {
  d <- design_split_plot(c("T1", "T2"), c("S1", "S2"), blocks = 2L, seed = 1L)
  expect_true("whole_plot_level" %in% names(d$layout))
  expect_true("sub_plot_level"   %in% names(d$layout))
})

test_that("design_split_plot() errors on < 2 blocks", {
  expect_error(
    design_split_plot(c("A", "B"), c("S1", "S2"), blocks = 1L)
  )
})

test_that("design_split_split_plot() returns correct class", {
  d <- design_split_split_plot(
    c("A", "B"), c("C", "D"), c("E", "F"), blocks = 3L, seed = 1L
  )
  expect_s3_class(d, "split_split_plot_design")
})

test_that("design_split_split_plot() has three error strata", {
  d <- design_split_split_plot(
    c("A", "B"), c("C", "D"), c("E", "F"), blocks = 3L, seed = 1L
  )
  es_names <- names(d$error_structure)
  expect_true("whole_plot_error"     %in% es_names)
  expect_true("sub_plot_error"       %in% es_names)
  expect_true("sub_sub_plot_error"   %in% es_names)
})

test_that("design_strip_plot() returns correct class", {
  d <- design_strip_plot(c("R1", "R2", "R3"), c("C1", "C2"), blocks = 3L,
                         seed = 1L)
  expect_s3_class(d, "strip_plot_design")
})

test_that("design_strip_plot() has three error strata", {
  d <- design_strip_plot(c("R1", "R2"), c("C1", "C2"), blocks = 3L, seed = 1L)
  es_names <- names(d$error_structure)
  expect_true("row_error"          %in% es_names)
  expect_true("col_error"          %in% es_names)
  expect_true("interaction_error"  %in% es_names)
})

# ============================================================================
# Step 19 Tests: Latin Square
# ============================================================================

test_that("design_latin_square() returns correct class", {
  d <- design_latin_square(c("A", "B", "C"), seed = 1L)
  expect_s3_class(d, "field_design")
  expect_s3_class(d, "latin_square_design")
})

test_that("design_latin_square() produces n×n layout", {
  d <- design_latin_square(c("A", "B", "C", "D"), seed = 1L)
  expect_equal(nrow(d$layout), 16L)
  expect_equal(ncol(d$layout), 5L)
})

test_that("design_latin_square() each treatment once per row", {
  d <- design_latin_square(c("A", "B", "C"), seed = 42L)
  for (r in unique(d$layout$row)) {
    trts <- d$layout$treatment[d$layout$row == r]
    expect_equal(sort(trts), sort(c("A", "B", "C")))
  }
})

test_that("design_latin_square() each treatment once per column", {
  d <- design_latin_square(c("A", "B", "C"), seed = 42L)
  for (cc in unique(d$layout$column)) {
    trts <- d$layout$treatment[d$layout$column == cc]
    expect_equal(sort(trts), sort(c("A", "B", "C")))
  }
})

test_that("design_latin_square() orthogonality: unique in both row and col", {
  d <- design_latin_square(c("A", "B", "C", "D"), seed = 7L)
  sq <- matrix(d$layout$treatment[order(d$layout$row, d$layout$column)],
               nrow = 4L, ncol = 4L)
  # Each row: 4 distinct treatments
  expect_true(all(apply(sq, 1L, function(r) length(unique(r))) == 4L))
  # Each col: 4 distinct treatments
  expect_true(all(apply(sq, 2L, function(cc) length(unique(cc))) == 4L))
})

test_that("design_latin_square() is reproducible", {
  d1 <- design_latin_square(c("T1", "T2", "T3"), seed = 1L)
  d2 <- design_latin_square(c("T1", "T2", "T3"), seed = 1L)
  expect_identical(d1$layout$treatment, d2$layout$treatment)
})

test_that("design_latin_square() errors on single treatment", {
  expect_error(design_latin_square(c("A")))
})

# ============================================================================
# Step 20 Tests: ANOVA Engine
# ============================================================================

test_that("analyze_design() returns agristat_aov for CRD", {
  set.seed(1L)
  d <- design_crd(c("A", "B", "C"), replicates = 5L, seed = 1L)
  df <- d$layout
  df$yield <- c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 11))
  m <- analyze_design(d, response = "yield", data = df)
  expect_s3_class(m, "agristat_aov")
})

test_that("analyze_design() CRD ANOVA has treatment term", {
  set.seed(1L)
  d <- design_crd(c("A", "B", "C"), replicates = 5L, seed = 1L)
  df <- d$layout
  df$yield <- c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 11))
  m <- analyze_design(d, response = "yield", data = df)
  sm <- summary(m$model)
  # summary(aov) returns a list; first element has rownames with terms
  term_names <- rownames(sm[[1L]])
  expect_true(any(grepl("treatment", term_names, ignore.case = TRUE)))
})

test_that("analyze_design() RCBD includes block term", {
  set.seed(2L)
  d <- design_rcbd(c("A", "B", "C"), blocks = 4L, seed = 2L)
  df <- d$layout
  df$yield <- rnorm(nrow(df), 10, 2)
  m <- analyze_design(d, response = "yield", data = df)
  expect_s3_class(m, "agristat_aov")
  sm <- summary(m$model)
  term_names <- rownames(sm[[1L]])
  expect_true(any(grepl("block", term_names, ignore.case = TRUE)))
  expect_true(any(grepl("treatment", term_names, ignore.case = TRUE)))
})

test_that("analyze_design() factorial has interaction term", {
  set.seed(3L)
  d <- design_factorial(
    factor_names  = c("N", "V"),
    factor_levels = list(N = c("Low", "High"), V = c("V1", "V2")),
    replicates    = 3L, seed = 3L
  )
  df <- d$layout
  df$yield <- rnorm(nrow(df), 10, 2)
  m <- analyze_design(d, response = "yield", data = df)
  expect_s3_class(m, "agristat_aov")
  sm <- summary(m$model)
  term_names <- rownames(sm[[1L]])
  expect_true(any(grepl("N:V|N\\*V", term_names)))
})

test_that("analyze_design() Latin Square has row + column + treatment terms", {
  set.seed(4L)
  d <- design_latin_square(c("A", "B", "C", "D"), seed = 4L)
  df <- d$layout
  df$yield <- rnorm(nrow(df), 10, 2)
  m <- analyze_design(d, response = "yield", data = df)
  expect_s3_class(m, "agristat_aov")
  sm <- summary(m$model)
  term_names <- rownames(sm[[1L]])
  expect_true(any(grepl("row",       term_names, ignore.case = TRUE)))
  expect_true(any(grepl("column",    term_names, ignore.case = TRUE)))
  expect_true(any(grepl("treatment", term_names, ignore.case = TRUE)))
})

test_that("analyze_design() errors on non-field_design input", {
  expect_error(analyze_design(list(), response = "y", data = data.frame(y = 1)))
})

test_that("analyze_design() errors when response column missing", {
  d <- design_crd(c("A", "B"), replicates = 3L, seed = 1L)
  df <- d$layout
  expect_error(analyze_design(d, response = "yield", data = df))
})

# ============================================================================
# Step 21 Tests: Error Term Helpers
# ============================================================================

test_that("compute_error_terms() returns named numeric vector", {
  d  <- design_split_plot(c("I+", "I-"), c("V1", "V2", "V3"),
                          blocks = 4L, seed = 1L)
  et <- compute_error_terms(d)
  expect_true(is.numeric(et))
  expect_true(!is.null(names(et)))
  expect_true(all(et > 0))
})

test_that("compute_error_terms() errors on non-field_design", {
  expect_error(compute_error_terms(list()))
})

test_that("get_error_structure() returns list", {
  d <- design_rcbd(c("A", "B", "C"), blocks = 3L, seed = 1L)
  es <- get_error_structure(d)
  expect_type(es, "list")
})

test_that("get_error_structure() split-plot has two strata", {
  d  <- design_split_plot(c("I+", "I-"), c("V1", "V2"), blocks = 3L, seed = 1L)
  es <- get_error_structure(d)
  expect_true("whole_plot_error" %in% names(es))
  expect_true("sub_plot_error"   %in% names(es))
})

test_that("extract_variance_components() errors on non-lmerMod", {
  d <- design_crd(c("A", "B"), replicates = 3L, seed = 1L)
  df <- d$layout
  df$yield <- rnorm(6)
  m <- analyze_design(d, response = "yield", data = df)
  # CRD uses aov, not lmerMod — should warn and return NA
  expect_warning(vc <- extract_variance_components(m$model))
  expect_true(is.na(vc))
})

# ============================================================================
# Step 22 Tests: Visualization
# ============================================================================

test_that("plot_field_layout() returns a ggplot object for CRD", {
  d <- design_crd(c("A", "B", "C"), replicates = 4L, seed = 1L)
  p <- plot_field_layout(d)
  expect_s3_class(p, "gg")
})

test_that("plot_field_layout() returns a ggplot object for RCBD", {
  d <- design_rcbd(c("V1", "V2", "V3", "V4"), blocks = 3L, seed = 1L)
  p <- plot_field_layout(d)
  expect_s3_class(p, "gg")
})

test_that("plot_field_layout() returns a ggplot object for Latin Square", {
  d <- design_latin_square(c("A", "B", "C"), seed = 1L)
  p <- plot_field_layout(d)
  expect_s3_class(p, "gg")
})

test_that("plot_field_layout() errors on non-field_design", {
  expect_error(plot_field_layout(list()))
})

test_that("plot_field_layout() accepts custom title", {
  d <- design_crd(c("A", "B"), replicates = 3L, seed = 1L)
  p <- plot_field_layout(d, title = "My Layout")
  expect_equal(p$labels$title, "My Layout")
})

test_that("plot_anova_diagnostics() runs without error for aov model", {
  set.seed(1L)
  d  <- design_crd(c("A", "B", "C"), replicates = 5L, seed = 1L)
  df <- d$layout
  df$yield <- c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 11))
  m  <- analyze_design(d, response = "yield", data = df)
  expect_silent(plot_anova_diagnostics(m, which = 1L:2L))
})

# ============================================================================
# S3 Methods Tests
# ============================================================================

test_that("print.field_design() produces output", {
  d <- design_crd(c("A", "B", "C"), replicates = 3L, seed = 1L)
  expect_output(print(d), "CRD")
})

test_that("summary.field_design() produces output with error structure", {
  d <- design_split_plot(c("I+", "I-"), c("V1", "V2"), blocks = 3L, seed = 1L)
  expect_output(summary(d), "error")
})

test_that("print.agristat_aov() mentions design type", {
  set.seed(1L)
  d  <- design_crd(c("A", "B", "C"), replicates = 5L, seed = 1L)
  df <- d$layout; df$yield <- rnorm(15, 10, 2)
  m  <- analyze_design(d, response = "yield", data = df)
  expect_output(print(m), "CRD")
})

# ============================================================================
# Edge Cases
# ============================================================================

test_that("design_crd() with 2 treatments, 1 replicate works", {
  d <- design_crd(c("A", "B"), replicates = 1L, seed = 1L)
  expect_equal(nrow(d$layout), 2L)
})

test_that("design_rcbd() with 2 blocks works", {
  d <- design_rcbd(c("A", "B", "C"), blocks = 2L, seed = 1L)
  expect_equal(nrow(d$layout), 6L)
})

test_that("design_factorial() single replicate works", {
  d <- design_factorial(
    factor_names  = c("A", "B"),
    factor_levels = list(A = c("a1", "a2"), B = c("b1", "b2")),
    replicates    = 1L, seed = 1L
  )
  expect_equal(nrow(d$layout), 4L)
})

test_that("design_crd() with numeric treatments", {
  d <- design_crd(treatments = c(1, 2, 3), replicates = 2L, seed = 1L)
  expect_equal(nrow(d$layout), 6L)
})
