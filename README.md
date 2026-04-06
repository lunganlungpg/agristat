# agristat

<!-- badges: start -->
[![R-CMD-check](https://github.com/lunganlungpg/agristat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lunganlungpg/agristat/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**agristat** is a comprehensive R package providing a complete biostatistical
workflow for agricultural and genetic research — from experimental design
through publication-ready output.

## Modules

| # | Module | Key Functions |
|---|---|---|
| 1 | **Data Ingestion & Pre-Processing** | `import_field_data()`, `declare_hypothesis()`, `check_assumptions()`, `auto_transform()`, `chi_square_segregation()` |
| 2 | **Field Layout Engine** | `design_crd()`, `design_rcbd()`, `design_factorial()`, `design_split_plot()`, `analyze_design()`, `plot_field_layout()` |
| 3 | **Breeding & Genetics Analytics** | `fit_ammi()`, `analyze_gxe()`, `calc_heritability()`, `run_basic_gwas()`, `plot_manhattan()` |
| 4 | **Mean Separation & Post-Hoc Tests** | `dmrt()`, `tukey_hsd()`, `dunnett_test()`, `mean_separation()`, `plot_means()` |
| 5 | **Phylogenetics & Evolutionary Mapping** | `read_fasta()`, `align_sequences()`, `build_tree_nj()`, `build_tree_ml()`, `bootstrap_tree()`, `plot_phylo()` |
| 6 | **Thesis Assembly & Export** | `generate_methods_text()`, `export_results_table()`, `export_figure()`, `compile_report()`, `thesis_template_methods()` |

## Installation

``` r
# install.packages("devtools")
devtools::install_github("lunganlungpg/agristat")
```

## Quick Start

``` r
library(agristat)

# ── Module 1: Declare a hypothesis ──────────────────────────────────────────
h <- declare_hypothesis(
  h0    = "No significant difference in yield between varieties",
  h1    = "At least one variety differs in yield",
  alpha = 0.05
)
print(h)

# ── Module 2: Design a field experiment ─────────────────────────────────────
design <- design_rcbd(
  treatments = c("Control", "T1", "T2", "T3"),
  blocks     = 4L,
  seed       = 42L
)
plot_field_layout(design)

# ── Module 2: Analyse the design ─────────────────────────────────────────────
df <- design$layout
df$yield <- rnorm(nrow(df), mean = 5, sd = 0.8)
aov_res  <- analyze_design(design, response = "yield", data = df)
summary(aov_res)

# ── Module 4: Post-hoc mean separation ──────────────────────────────────────
trt_means <- tapply(df$yield, df$treatment, mean)
ph <- dmrt(trt_means, mse = 0.64, df = 9, n = 4)
print(ph)

# ── Module 6: Export results ─────────────────────────────────────────────────
mt  <- generate_methods_text(design, aov_res, ph)
cat(mt$text)
cat(export_results_table(aov_res, format = "markdown"))
```

## Vignettes

- [01 Data Ingestion & Pre-Processing](vignettes/01-data-ingestion.Rmd)
- [02 Field Designs](vignettes/02-field-designs.Rmd)
- [03 Genetics & GWAS](vignettes/03-genetics-gwas.Rmd)
- [04 Mean Separation](vignettes/04-mean-separation.Rmd)
- [05 Phylogenetics](vignettes/05-phylogenetics.Rmd)
- [06 Thesis Assembly & Export](vignettes/06-thesis-export.Rmd)

## License

GPL (>= 3).
