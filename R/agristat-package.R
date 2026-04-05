#' agristat: Comprehensive Biostatistical Suite for Agricultural and Genetic Research
#'
#' @description
#' The `agristat` package provides a comprehensive suite of biostatistical
#' tools designed for agricultural scientists and plant breeders. It covers
#' six modules:
#'
#' \enumerate{
#'   \item **Data Ingestion & Pre-Processing** — hypothesis tracking,
#'     assumption testing, data transformations, field data import.
#'   \item **Field Layout Engine** — experimental design generators (CRD,
#'     RCBD, factorial, split-plot, Latin square) and ANOVA analysis.
#'   \item **Breeding & Genetics Analytics** — AMMI, GxE analysis,
#'     heritability estimation, QTL mapping, GWAS.
#'   \item **Mean Separation & Post-Hoc Tests** — DMRT, Tukey HSD, Dunnett,
#'     orthogonal contrasts, publication-ready visualisation.
#'   \item **Phylogenetics** — sequence alignment, tree building,
#'     bootstrapping.
#'   \item **Thesis Export** — methods text generation, formatted tables,
#'     figures for academic manuscripts.
#' }
#'
#' @section Module 4 — Mean Separation & Post-Hoc Tests:
#' After a significant ANOVA F-test, use one of:
#' \itemize{
#'   \item [dmrt()] — Duncan's Multiple Range Test (maximises power)
#'   \item [tukey_hsd()] — Tukey's HSD (controls FWER)
#'   \item [dunnett_test()] — Dunnett's test (vs a control group)
#'   \item [test_contrasts()] — Orthogonal planned contrasts
#'   \item [mean_separation()] — Unified dispatcher
#' }
#' Visualise results with [plot_means()] and [plot_tukey()].
#'
#' @keywords internal
"_PACKAGE"
