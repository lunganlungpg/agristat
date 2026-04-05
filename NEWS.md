# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## agristat 1.0.0

*Phase 2 — Module 1: Data Ingestion & Pre-Processing*

### Added

#### Module 1 — Data Ingestion & Pre-Processing (Steps 8–14)

- `declare_hypothesis()` — S3 class `hypothesis` for tracking null and
  alternative hypotheses with significance level and declaration timestamp
- `print.hypothesis()`, `summary.hypothesis()` — S3 display methods
- `check_assumptions()` — Wrapper running Shapiro-Wilk normality test and
  (optionally) Levene's test for homogeneity of variance; returns S3 object
  `assumption_report` with pass/fail flags and plain-language recommendation
- `auto_transform()` — Tests log, sqrt, inverse, and Box-Cox transformations;
  selects the best based on Shapiro-Wilk p-value; returns comparison ggplot2 plot
- `suggest_nonparametric()` — Examines `assumption_report` and recommends
  Mann-Whitney U, Kruskal-Wallis, Welch's t-test, or Welch ANOVA
- `chi_square_segregation()` — Mendelian segregation chi-square test supporting
  any ratio (3:1, 9:3:3:1, 15:1, etc.)
- `import_field_data()` — Reads CSV, TSV, Excel, GeoJSON, and Shapefile with
  automatic type detection and missing-value standardisation
- `format_p_value()`, `format_statistic()` — Publication-ready formatting helpers
- Comprehensive test suite in `tests/testthat/test-data-ingestion.R`
- Vignette `vignettes/01-data-ingestion.Rmd` with complete workflow

#### Package Infrastructure (Phase 1)

- Complete `DESCRIPTION` with GPL-3 licence and all module dependencies
- NAMESPACE configured for roxygen2 auto-generation
- Package-level documentation in `R/agristat-package.R`
- Testing infrastructure using testthat 3rd edition
- 4 GitHub Actions workflows (R-CMD-check, coverage, pkgdown, lint)
- `_pkgdown.yml` with Bootstrap 5 site
- `README.md`, `NEWS.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`
