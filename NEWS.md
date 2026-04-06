# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## agristat 1.0.0

### Added

#### Module 6 — Thesis Assembly & Export (Phase 7, Steps 50–57)

- `generate_methods_text()` — auto-generates a publication-quality methods
  paragraph from design, analysis, and post-hoc objects; returns S3 class
  `methods_text` with `$text`, `$sections`, and `$references`
- `export_results_table()` — renders ANOVA, means, or GWAS results as
  Markdown, HTML, LaTeX, or CSV; supports writing directly to a file
- `export_figure()` — saves ggplot2 or base-R graphics to PNG, PDF, SVG, or
  EPS at configurable DPI (default 300); returns file path, size, and format
- `generate_appendix()` — produces a structured appendix with descriptive
  statistics, Shapiro-Wilk residual normality test, and raw data table;
  returns S3 class `appendix`
- `compile_report()` — assembles methods text, results tables, figures, and
  appendix into an Rmd document; renders via `rmarkdown::render()` when
  available, otherwise writes a plain Markdown file
- `thesis_template_methods()`, `thesis_template_results()`,
  `thesis_template_discussion()` — opinionated academic section templates
  with `{placeholder}` tokens; return S3 class `thesis_template`
- `fill_template()` — substitutes `{placeholder}` tokens with supplied values
- `print.methods_text()`, `print.appendix()`, `print.thesis_template()` —
  S3 display methods
- Test suite in `tests/testthat/test-thesis-export.R` (45 tests)
- Vignette `vignettes/06-thesis-export.Rmd` with end-to-end rice yield trial

#### Module 5 — Phylogenetics & Evolutionary Mapping (Phase 6, Steps 34–49)

- `read_fasta()` — parses FASTA files; auto-detects nucleotide / protein
  type; returns S3 class `dna_sequences` or `protein_sequences`
- `validate_sequence_format()` — validates character alphabet and length
- `align_sequences()` — performs multiple sequence alignment; returns S3
  class `aligned_sequences`
- `compute_distance_matrix()` — computes pairwise evolutionary distances
- `build_tree_nj()` — Neighbour-Joining tree construction; returns S3 class
  `phylo_tree_nj`
- `build_tree_ml()` — Maximum-Likelihood tree construction; returns S3 class
  `phylo_tree_ml`
- `bootstrap_tree()` — bootstrap resampling for branch support values
- `extract_bootstrap_values()` — extracts bootstrap support from a `phylo`
  object
- `plot_phylo()` — plots a rooted phylogenetic tree with bootstrap labels
- `plot_tree_unrooted()` — unrooted (circular) tree plot
- `reroot_tree()` — midpoint or outgroup re-rooting
- `write_tree()` — exports tree in Newick or Nexus format
- Test suite in `tests/testthat/test-phylogenetics.R`
- Vignette `vignettes/05-phylogenetics.Rmd`

#### Module 4 — Mean Separation & Post-Hoc Tests (Phase 5, Steps 26–33)

- `dmrt()` — Duncan's Multiple Range Test; returns S3 class `dmrt_result`
  with means table and critical ranges; supports unequal replication via
  harmonic-mean correction
- `tukey_hsd()` — Tukey's Honestly Significant Difference; returns S3 class
  `tukey_result` with pairwise differences and confidence intervals
- `dunnett_test()` — Dunnett's many-to-one comparison with a control;
  returns S3 class `dunnett_result`
- `test_contrasts()` — user-defined linear contrasts; returns S3 class
  `contrast_result`
- `mean_separation()` — dispatcher that routes to DMRT, Tukey, or Dunnett
  based on the `method` argument
- `plot_means()` — bar or point plot of treatment means with error bars and
  significance letters
- `plot_tukey()` — compact letter display plot for Tukey results
- Helpers: `assign_letters()`, `assign_letters_from_matrix()`,
  `validate_contrasts()`, `format_contrast_output()` (exported)
- Test suite in `tests/testthat/test-posthoc.R`
- Vignette `vignettes/04-mean-separation.Rmd`

#### Module 3 — Breeding & Genetics Analytics (Phase 4, Steps 18–25)

- `design_alpha_lattice()` — generates alpha-lattice incomplete block designs
- `design_incomplete_block()` — balanced incomplete block design generator
- `design_augmented()` — augmented design with checks and unreplicated entries
- `analyze_gxe()` — genotype-by-environment interaction analysis; returns S3
  class `gxe_result` with stability statistics (Eberhart-Russell)
- `fit_ammi()` — AMMI (Additive Main Effects and Multiplicative Interaction)
  model; returns S3 class `ammi_result` with IPCAs and variance partitioning
- `plot_ammi_biplot()` — biplot of genotype and environment scores
- `calc_heritability()` — broad-sense heritability via ANOVA components;
  returns S3 class `heritability_result`
- `prepare_qtl_data()` — data preparation for QTL mapping
- `plot_qtl_map()` — visualises QTL positions on a linkage map
- `run_basic_gwas()` — simple marker-trait association scan; returns S3 class
  `gwas_result`
- `plot_manhattan()` — Manhattan plot with genome-wide significance threshold
- `plot_qq_gwas()` — QQ plot for GWAS p-values
- Test suite in `tests/testthat/test-genetics.R`
- Vignette `vignettes/03-genetics-gwas.Rmd`

#### Module 2 — Field Layout Engine (Phase 3, Steps 8–17)

- `design_crd()` — Completely Randomised Design layout with randomisation
- `design_rcbd()` — Randomised Complete Block Design
- `design_factorial()` — two-way or multi-way factorial design
- `design_split_plot()` — split-plot design with whole-plot and sub-plot factors
- `design_split_split_plot()` — three-level split-split-plot design
- `design_strip_plot()` — strip-plot (criss-cross) design
- `design_latin_square()` — Latin Square design
- All design functions return S3 class `field_design` via `new_field_design()`
- `analyze_design()` — fits the appropriate ANOVA model for a `field_design`;
  returns S3 class `agristat_aov` with `$model`, `$design`, `$response`
- `plot_field_layout()` — ggplot2 tile map of the randomised field layout
- `compute_error_terms()` — extracts degrees of freedom for each error stratum
- `get_error_structure()` — returns the error structure for a design type
- `extract_variance_components()` — variance component estimation
- Test suite in `tests/testthat/test-field-layouts.R`
- Vignette `vignettes/02-field-designs.Rmd`

#### Module 1 — Data Ingestion & Pre-Processing (Phase 2, Steps 1–7)

- `declare_hypothesis()` — S3 class `hypothesis` for tracking null and
  alternative hypotheses with significance level and declaration timestamp
- `print.hypothesis()`, `summary.hypothesis()` — S3 display methods
- `check_assumptions()` — Shapiro-Wilk normality and Levene's homogeneity of
  variance tests; returns S3 class `assumption_report` with pass/fail flags
- `auto_transform()` — evaluates log, sqrt, inverse, and Box-Cox
  transformations; selects best by Shapiro-Wilk p-value
- `suggest_nonparametric()` — recommends Mann-Whitney U, Kruskal-Wallis,
  Welch's t-test, or Welch ANOVA based on assumption failures
- `chi_square_segregation()` — Mendelian segregation chi-square test for any
  ratio (3:1, 9:3:3:1, 15:1, etc.)
- `import_field_data()` — reads CSV, TSV, Excel, GeoJSON, and Shapefile with
  automatic type detection and missing-value standardisation
- `format_p_value()`, `format_statistic()`, `format_effect_size()` —
  publication-ready formatting helpers
- Test suite in `tests/testthat/test-data-ingestion.R`
- Vignette `vignettes/01-data-ingestion.Rmd`

#### Package Infrastructure (Phase 1)

- Complete `DESCRIPTION` with GPL-3 licence and all module dependencies
- NAMESPACE configured for roxygen2 auto-generation
- Package-level documentation in `R/agristat-package.R`
- Testing infrastructure using testthat 3rd edition
- 4 GitHub Actions workflows (R-CMD-check, coverage, pkgdown, lint)
- `_pkgdown.yml` with Bootstrap 5 pkgdown site
- `README.md`, `NEWS.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`
