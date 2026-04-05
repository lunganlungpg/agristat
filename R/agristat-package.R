#' @keywords internal
"_PACKAGE"

#' agristat: Comprehensive Biostatistical Suite for Agricultural and Genetic Research
#'
#' @description
#' The `agristat` package provides a comprehensive suite of biostatistical tools
#' tailored for agricultural and genetic research. It implements six integrated
#' modules covering the full analytical workflow from data ingestion to thesis-ready
#' output.
#'
#' ## Modules
#'
#' The package is organised into six modules:
#'
#' ### Module 1 — Data Ingestion & Pre-Processing (`data-ingestion`)
#' Functions for importing field trial data, declaring and tracking hypotheses,
#' testing statistical assumptions (normality, homogeneity of variance), and
#' suggesting appropriate transformations or non-parametric alternatives.
#'
#' ### Module 2 — Field Layout Engine (`field-layouts`)
#' Tools for constructing and visualising agricultural experimental designs
#' including Completely Randomised Design (CRD), Randomised Complete Block Design
#' (RCBD), Latin Square, Split-Plot, and augmented designs.
#'
#' ### Module 3 — Breeding & Genetics Analytics (`genetics`)
#' Implementation of AMMI (Additive Main Effects and Multiplicative Interaction)
#' analysis, GGE biplots, heritability estimation, combining ability (GCA/SCA),
#' and Mendelian segregation tests.
#'
#' ### Module 4 — Mean Separation & Post-Hoc Tests (`posthoc`)
#' A unified interface to multiple comparison procedures including Tukey HSD,
#' Duncan's Multiple Range Test, Fisher's LSD, Scott-Knott, and non-parametric
#' alternatives.
#'
#' ### Module 5 — Phylogenetics & Evolutionary Mapping (`phylogenetics`)
#' Wrappers around `ape` and `phangorn` for constructing phylogenetic trees,
#' computing genetic distances, and mapping traits onto evolutionary trees.
#'
#' ### Module 6 — Thesis Assembly Export Suite (`thesis-export`)
#' Utilities to generate publication-ready tables, figures, and structured
#' R Markdown / Word documents suitable for thesis chapters and journal
#' submissions.
#'
#' ## Getting Started
#'
#' ```r
#' library(agristat)
#' ```
#'
#' See the package vignettes for worked examples:
#' `browseVignettes("agristat")`
#'
#' @author Lung An Lung Pg \email{lunganlungpg@@example.com}
#'
#' @docType package
#' @name agristat-package
#' @aliases agristat
NULL
