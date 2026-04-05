# tests/testthat/test-genetics.R
# Tests for Phase 4: Module 3 - Breeding & Genetics Analytics
# Steps 25-33: design functions, AMMI, heritability, QTL, GWAS

# ============================================================================
# Step 25: Alpha-Lattice and BIBD Design Tests
# ============================================================================

test_that("design_alpha_lattice: basic structure", {
  d <- design_alpha_lattice(20L, 4L, 3L, seed = 1L)
  expect_s3_class(d, "data.frame")
  expect_named(d, c("rep", "block", "plot_id", "genotype"))
  expect_equal(nrow(d), 60L)
})

test_that("design_alpha_lattice: every genotype appears in each replicate", {
  d <- design_alpha_lattice(12L, 3L, 2L, seed = 42L)
  for (r in 1:2) {
    genos_in_rep <- d[d == r]
    expect_equal(length(unique(genos_in_rep)), 12L)
    expect_equal(length(genos_in_rep), 12L)
  }
})

test_that("design_alpha_lattice: block sizes are correct", {
  d <- design_alpha_lattice(16L, 4L, 3L, seed = 5L)
  block_sizes <- table(d, d)
  expect_true(all(block_sizes == 4L))
})

test_that("design_alpha_lattice: plot_id is sequential", {
  d <- design_alpha_lattice(8L, 4L, 2L, seed = 7L)
  expect_equal(d, seq_len(nrow(d)))
})

test_that("design_alpha_lattice: reproducible with seed", {
  d1 <- design_alpha_lattice(20L, 4L, 2L, seed = 99L)
  d2 <- design_alpha_lattice(20L, 4L, 2L, seed = 99L)
  expect_identical(d1, d2)
})

test_that("design_alpha_lattice: error on non-divisible block size", {
  expect_error(design_alpha_lattice(10L, 3L, 2L), "evenly")
})

test_that("design_alpha_lattice: error on too-few genotypes", {
  expect_error(design_alpha_lattice(3L, 1L, 2L))
})

test_that("design_incomplete_block: returns field_design", {
  genos <- paste0("G", 1:7)
  d <- design_incomplete_block(genos, block_size = 3L, replicates = 3L,
                               seed = 1L)
  expect_s3_class(d, "field_design")
  expect_equal(d, "BIBD")
})

test_that("design_incomplete_block: layout has correct columns", {
  d <- design_incomplete_block(paste0("G", 1:7), 3L, 3L, seed = 2L)
  expect_named(d, c("block", "plot_id", "genotype"))
})

test_that("design_incomplete_block: each genotype appears r times", {
  genos <- paste0("G", 1:7)
  d <- design_incomplete_block(genos, block_size = 7L, replicates = 3L,
                               seed = 3L)
  counts <- table(d)
  expect_true(all(counts == 3L))
})

test_that("design_incomplete_block: parameters list is populated", {
  d <- design_incomplete_block(paste0("G", 1:7), 3L, 3L, seed = 4L)
  expect_equal(d, 7L)
  expect_equal(d, 3L)
  expect_equal(d, 3L)
})

test_that("design_incomplete_block: error when block_size >= length(genotypes)", {
  expect_error(design_incomplete_block(paste0("G", 1:5), 5L, 2L))
})

# ============================================================================
# Step 26: Augmented Design Tests
# ============================================================================

test_that("design_augmented: basic structure", {
  d <- design_augmented(paste0("T", 1:10), c("C1", "C2"), 4L, seed = 1L)
  expect_s3_class(d, "data.frame")
  expect_named(d, c("block", "plot_id", "genotype", "type"))
})

test_that("design_augmented: type column has only check/test", {
  d <- design_augmented(paste0("T", 1:8), c("C1", "C2"), 3L, seed = 2L)
  expect_true(all(d %in% c("check", "test")))
})

test_that("design_augmented: checks appear in every block", {
  d <- design_augmented(paste0("T", 1:12), c("C1", "C2", "C3"),
                        replicates_checks = 4L, seed = 3L)
  checks <- d[d == "check", ]
  for (b in 1:4) {
    block_checks <- checks[checks == b]
    expect_true(all(c("C1", "C2", "C3") %in% block_checks))
  }
})

test_that("design_augmented: test genotypes appear at most once", {
  tests <- paste0("T", 1:10)
  d <- design_augmented(tests, c("C1", "C2"), 4L, seed = 5L)
  test_counts <- table(d[d == "test"])
  expect_true(all(test_counts == 1L))
})

test_that("design_augmented: all test genotypes are included", {
  tests <- paste0("T", 1:15)
  d <- design_augmented(tests, c("C1"), 3L, seed = 6L)
  included_tests <- unique(d[d == "test"])
  expect_true(all(tests %in% included_tests))
})

test_that("design_augmented: error on empty test_genotypes", {
  expect_error(design_augmented(character(0), c("C1"), 2L))
})

# ============================================================================
# Step 27: GxE Analysis Tests
# ============================================================================

test_that("analyze_gxe: returns gxe_result S3 object", {
  set.seed(10L)
  df <- expand.grid(genotype    = paste0("G", 1:5),
                    environment = paste0("E", 1:4),
                    stringsAsFactors = FALSE)
  df <- rnorm(nrow(df), 5, 1)
  res <- analyze_gxe(df, response = "yield")
  expect_s3_class(res, "gxe_result")
})

test_that("analyze_gxe: output components present", {
  set.seed(11L)
  df <- expand.grid(genotype    = paste0("G", 1:4),
                    environment = paste0("E", 1:3),
                    stringsAsFactors = FALSE)
  df <- rnorm(nrow(df), 5, 1)
  res <- analyze_gxe(df, response = "yield")
  expect_true(all(c("genotype_effects", "env_effects",
                    "interaction_matrix", "stability",
                    "anova_table", "response") %in% names(res)))
})

test_that("analyze_gxe: genotype effects sum to approximately zero", {
  set.seed(12L)
  df <- expand.grid(genotype    = paste0("G", 1:5),
                    environment = paste0("E", 1:3),
                    stringsAsFactors = FALSE)
  df <- rnorm(nrow(df), 5, 1)
  res <- analyze_gxe(df, response = "yield")
  expect_true(abs(sum(res)) < 1e-9)
})

test_that("analyze_gxe: stability data frame has correct columns", {
  set.seed(13L)
  df <- expand.grid(genotype    = paste0("G", 1:3),
                    environment = paste0("E", 1:4),
                    stringsAsFactors = FALSE)
  df <- rnorm(nrow(df), 5, 1)
  res <- analyze_gxe(df, response = "yield")
  expect_named(res, c("genotype", "env_variance", "mean_yield"))
})

test_that("analyze_gxe: single environment does not error", {
  df <- data.frame(genotype    = c("G1", "G2", "G3"),
                   environment = "E1",
                   yield       = c(5.1, 4.8, 5.3))
  res <- analyze_gxe(df, response = "yield")
  expect_s3_class(res, "gxe_result")
})

test_that("analyze_gxe: print method works", {
  set.seed(14L)
  df <- expand.grid(genotype    = paste0("G", 1:3),
                    environment = paste0("E", 1:3),
                    stringsAsFactors = FALSE)
  df <- rnorm(nrow(df))
  res <- analyze_gxe(df, response = "yield")
  expect_output(print(res), "GxE Interaction")
})

# ============================================================================
# Step 28: AMMI Analysis Tests
# ============================================================================

test_that("fit_ammi: matrix input returns ammi_result", {
  mat <- matrix(c(10,12,11, 9,14,13, 8,11,10), nrow = 3,
                dimnames = list(paste0("G", 1:3), paste0("E", 1:3)))
  am  <- fit_ammi(mat)
  expect_s3_class(am, "ammi_result")
})

test_that("fit_ammi: variance_explained sums to 1", {
  mat <- matrix(rnorm(20), nrow = 4,
                dimnames = list(paste0("G", 1:4), paste0("E", 1:5)))
  am  <- fit_ammi(mat)
  expect_equal(sum(am), 1, tolerance = 1e-6)
})

test_that("fit_ammi: tidy data frame input works", {
  set.seed(20L)
  df <- expand.grid(genotype    = paste0("G", 1:4),
                    environment = paste0("E", 1:3),
                    stringsAsFactors = FALSE)
  df <- rnorm(nrow(df), 5, 1)
  am <- fit_ammi(df, response = "yield")
  expect_s3_class(am, "ammi_result")
  expect_equal(nrow(am),   4L)
  expect_equal(nrow(am), 3L)
})

test_that("fit_ammi: scores and loadings have correct IPCA column names", {
  mat <- matrix(rnorm(12), nrow = 3,
                dimnames = list(paste0("G", 1:3), paste0("E", 1:4)))
  am  <- fit_ammi(mat)
  expect_true(all(grepl("^IPCA", colnames(am))))
  expect_true(all(grepl("^IPCA", colnames(am))))
})

test_that("fit_ammi: interaction matrix dimensions correct", {
  mat <- matrix(rnorm(15), nrow = 3,
                dimnames = list(paste0("G", 1:3), paste0("E", 1:5)))
  am  <- fit_ammi(mat)
  expect_equal(dim(am), c(3L, 5L))
})

test_that("fit_ammi: SVD accuracy -- reconstruct interaction", {
  set.seed(22L)
  mat <- matrix(rnorm(20), nrow = 4,
                dimnames = list(paste0("G", 1:4), paste0("E", 1:5)))
  am  <- fit_ammi(mat)
  n_axes <- ncol(am)
  sv     <- am
  recon  <- am %*% diag(sv, n_axes, n_axes) %*% t(am)
  int_nona <- am
  int_nona[is.na(int_nona)] <- 0
  expect_equal(recon, int_nona, tolerance = 1e-8)
})

test_that("plot_ammi_biplot: returns ggplot object", {
  set.seed(30L)
  df <- expand.grid(genotype    = paste0("G", 1:5),
                    environment = paste0("E", 1:4),
                    stringsAsFactors = FALSE)
  df <- rnorm(nrow(df), 5, 1.5)
  am <- fit_ammi(df, response = "yield")
  p  <- plot_ammi_biplot(am)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ammi_biplot: error on out-of-range axis", {
  mat <- matrix(rnorm(9), nrow = 3,
                dimnames = list(paste0("G", 1:3), paste0("E", 1:3)))
  am  <- fit_ammi(mat)
  expect_error(plot_ammi_biplot(am, pc_x = 10L, pc_y = 2L), "exceed")
})

# ============================================================================
# Step 29: Heritability Calculator Tests
# ============================================================================

test_that("calc_heritability: known H2 value", {
  vc <- data.frame(source   = c("genotype", "residual"),
                   variance = c(2.5, 1.0))
  h  <- calc_heritability(vc, n_boot = 0L)
  expect_equal(h, 2.5 / 3.5, tolerance = 1e-9)
})

test_that("calc_heritability: returns heritability_result class", {
  vc <- data.frame(source   = c("genotype", "residual"),
                   variance = c(1.0, 1.0))
  h  <- calc_heritability(vc, n_boot = 0L)
  expect_s3_class(h, "heritability_result")
})

test_that("calc_heritability: h2 is NA when additive variance absent", {
  vc <- data.frame(source   = c("genotype", "residual"),
                   variance = c(1.5, 0.5))
  h  <- calc_heritability(vc, n_boot = 0L)
  expect_true(is.na(h))
})

test_that("calc_heritability: narrow-sense h2 when additive given", {
  vc <- data.frame(source   = c("genotype", "additive", "residual"),
                   variance = c(2.0, 1.0, 1.0))
  h  <- calc_heritability(vc, n_boot = 0L)
  expect_equal(h, 1.0 / 3.0, tolerance = 1e-9)
})

test_that("calc_heritability: accepts 'error' as synonym for 'residual'", {
  vc <- data.frame(source   = c("genotype", "error"),
                   variance = c(3.0, 1.0))
  h  <- calc_heritability(vc, n_boot = 0L)
  expect_equal(h, 0.75, tolerance = 1e-9)
})

test_that("calc_heritability: bootstrap CI is ordered and in [0,1]", {
  set.seed(40L)
  vc <- data.frame(source   = c("genotype", "residual"),
                   variance = c(2.0, 1.0))
  h  <- calc_heritability(vc, n_boot = 200L, seed = 40L)
  expect_true(h[1L] <= h[2L])
  expect_true(all(h >= 0))
  expect_true(all(h <= 1))
})

test_that("calc_heritability: interpretation is non-empty string", {
  vc <- data.frame(source   = c("genotype", "residual"),
                   variance = c(4.0, 1.0))
  h  <- calc_heritability(vc, n_boot = 0L)
  expect_type(h, "character")
  expect_true(nchar(h) > 0L)
})

test_that("calc_heritability: error when genotype source missing", {
  vc <- data.frame(source = c("env", "residual"), variance = c(1, 1))
  expect_error(calc_heritability(vc, n_boot = 0L), "genotype")
})

test_that("calc_heritability: print method works", {
  vc <- data.frame(source   = c("genotype", "residual"),
                   variance = c(2.0, 1.0))
  h  <- calc_heritability(vc, n_boot = 0L)
  expect_output(print(h), "Heritability")
})

# ============================================================================
# Step 30: QTL Mapping Tests
# ============================================================================

test_that("prepare_qtl_data: basic output structure", {
  set.seed(50L)
  n <- 30L; m <- 8L
  pheno <- data.frame(id    = paste0("ind", seq_len(n)),
                      yield = rnorm(n))
  markers <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                    dimnames = list(paste0("ind", seq_len(n)),
                                    paste0("M", seq_len(m))))
  res <- prepare_qtl_data(pheno, markers)
  expect_true(all(c("phenotypes", "genotypes", "map_info",
                    "n_individuals", "n_markers", "qc_report") %in%
                    names(res)))
})

test_that("prepare_qtl_data: n_markers equals number of marker columns", {
  set.seed(51L)
  n <- 20L; m <- 5L
  pheno <- data.frame(id = paste0("ind", seq_len(n)), trait = rnorm(n))
  markers <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                    dimnames = list(paste0("ind", seq_len(n)),
                                    paste0("M", seq_len(m))))
  res <- prepare_qtl_data(pheno, markers)
  expect_equal(res, m)
})

test_that("prepare_qtl_data: A/H/B coding is converted to 0/1/2", {
  set.seed(52L)
  n <- 10L; m <- 4L
  pheno <- data.frame(id = paste0("i", seq_len(n)), y = rnorm(n))
  markers <- matrix(
    sample(c("A", "H", "B"), n * m, replace = TRUE),
    nrow     = n,
    dimnames = list(paste0("i", seq_len(n)), paste0("M", seq_len(m)))
  )
  res <- prepare_qtl_data(pheno, markers)
  expect_true(is.numeric(res))
  expect_true(all(res %in% c(0L, 1L, 2L), na.rm = TRUE))
})

test_that("prepare_qtl_data: qc_report has correct columns", {
  set.seed(53L)
  n <- 15L; m <- 6L
  pheno   <- data.frame(id = paste0("ind", seq_len(n)), y = rnorm(n))
  markers <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                    dimnames = list(paste0("ind", seq_len(n)),
                                    paste0("M", seq_len(m))))
  res <- prepare_qtl_data(pheno, markers)
  expect_named(res,
               c("marker", "missing_rate", "minor_af", "monomorphic"))
})

test_that("plot_qtl_map: returns ggplot", {
  set.seed(55L)
  scan <- data.frame(chr = rep(1:3, each = 20L),
                     pos = rep(seq(0, 95, 5), 3L),
                     lod = abs(rnorm(60L)))
  p <- plot_qtl_map(scan, threshold = 2.5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_qtl_map: error when required columns missing", {
  df <- data.frame(chr = 1:3, pos = 1:3)
  expect_error(plot_qtl_map(df), "lod")
})

# ============================================================================
# Step 31: GWAS Pipeline Tests
# ============================================================================

test_that("run_basic_gwas: returns gwas_result class", {
  set.seed(60L)
  n <- 80L; m <- 20L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  g <- run_basic_gwas(pheno, gmat)
  expect_s3_class(g, "gwas_result")
})

test_that("run_basic_gwas: results data frame has required columns", {
  set.seed(61L)
  n <- 80L; m <- 20L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  g     <- run_basic_gwas(pheno, gmat)
  expect_true(all(c("marker", "chromosome", "position",
                    "p_value", "neg_log10_p", "beta", "se", "maf") %in%
                    names(g)))
})

test_that("run_basic_gwas: n_individuals is correct", {
  set.seed(62L)
  n <- 50L; m <- 10L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  g <- run_basic_gwas(pheno, gmat)
  expect_equal(g, n)
})

test_that("run_basic_gwas: Bonferroni threshold is 0.05/n_tested", {
  set.seed(63L)
  n <- 60L; m <- 30L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE,
                          prob = c(0.25, 0.5, 0.25)),
                  nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  g <- run_basic_gwas(pheno, gmat)
  expect_equal(g, 0.05 / g,
               tolerance = 1e-12)
})

test_that("run_basic_gwas: p-values are in (0, 1]", {
  set.seed(64L)
  n <- 60L; m <- 25L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  g     <- run_basic_gwas(pheno, gmat)
  pvals <- g[!is.na(g)]
  expect_true(all(pvals > 0 & pvals <= 1))
})

test_that("run_basic_gwas: error when phenotype length mismatch", {
  n <- 50L; m <- 10L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, (n + 5L) * m, replace = TRUE),
                  nrow = n + 5L)
  expect_error(run_basic_gwas(pheno, gmat), "length")
})

test_that("run_basic_gwas: binary trait option works", {
  set.seed(65L)
  n <- 80L; m <- 15L
  pheno <- as.numeric(sample(0:1, n, replace = TRUE))
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  g <- run_basic_gwas(pheno, gmat, trait_type = "binary")
  expect_s3_class(g, "gwas_result")
  expect_equal(g, "binary")
})

test_that("plot_manhattan: returns ggplot", {
  set.seed(66L)
  n <- 100L; m <- 50L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE),
                  nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  pos <- data.frame(marker     = paste0("SNP", seq_len(m)),
                    chromosome = rep(1:5, each = 10L),
                    position   = rep(seq(1, 100, 10), 5L))
  g <- run_basic_gwas(pheno, gmat, pos)
  p <- plot_manhattan(g)
  expect_s3_class(p, "ggplot")
})

test_that("plot_manhattan: point count matches non-NA p-values", {
  set.seed(67L)
  n <- 80L; m <- 40L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE,
                          prob = c(0.25, 0.5, 0.25)),
                  nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  g     <- run_basic_gwas(pheno, gmat)
  p     <- plot_manhattan(g)
  # Extract data used for geom_point
  built  <- ggplot2::ggplot_build(p)
  layers <- built
  point_layer <- layers[[which(vapply(layers, function(l)
    "size" %in% names(l), logical(1L)))[1L]]]
  expect_equal(nrow(point_layer), n_pts)
})

test_that("plot_qq_gwas: numeric p-value vector accepted", {
  set.seed(68L)
  pvals <- runif(200L)
  p     <- plot_qq_gwas(pvals)
  expect_s3_class(p, "ggplot")
})

test_that("plot_qq_gwas: gwas_result object accepted", {
  set.seed(69L)
  n <- 60L; m <- 20L
  pheno <- rnorm(n)
  gmat  <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                  dimnames = list(NULL, paste0("SNP", seq_len(m))))
  g <- run_basic_gwas(pheno, gmat)
  p <- plot_qq_gwas(g)
  expect_s3_class(p, "ggplot")
})

test_that("plot_qq_gwas: error on invalid input type", {
  expect_error(plot_qq_gwas("not a valid input"), "gwas_result")
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("analyze_gxe: single genotype does not error", {
  df <- data.frame(genotype    = "G1",
                   environment = c("E1", "E2", "E3"),
                   yield       = c(4.0, 5.0, 6.0))
  expect_s3_class(analyze_gxe(df, response = "yield"), "gxe_result")
})

test_that("run_basic_gwas: all-missing marker is skipped", {
  set.seed(70L)
  n <- 40L; m <- 5L
  pheno   <- rnorm(n)
  gmat    <- matrix(sample(0:2, n * m, replace = TRUE,
                             prob = c(0.25, 0.5, 0.25)),
                    nrow = n,
                    dimnames = list(NULL, paste0("SNP", seq_len(m))))
  gmat[, 1L] <- NA_real_
  g <- run_basic_gwas(pheno, gmat, min_maf = 0.05)
  expect_true(is.na(g[g == "SNP1"]))
})

test_that("calc_heritability: H2 = 0.5 for equal variances", {
  vc <- data.frame(source   = c("genotype", "residual"),
                   variance = c(1.0, 1.0))
  h  <- calc_heritability(vc, n_boot = 0L)
  expect_equal(h, 0.5, tolerance = 1e-9)
})

test_that("format_effect_size: formats correctly", {
  s <- format_effect_size(0.123456, 0.012345, digits = 4L)
  expect_type(s, "character")
  expect_true(grepl("0\.1235", s))
})

test_that("prepare_qtl_data: error with no matching IDs", {
  pheno   <- data.frame(id = c("a", "b"), y = 1:2)
  markers <- matrix(0L, 2, 2,
                    dimnames = list(c("x", "z"), c("M1", "M2")))
  expect_error(prepare_qtl_data(pheno, markers), "No matching IDs")
})
