# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## agristat 1.0.0

*Initial release — foundational package infrastructure (Phase 1)*

### Added

#### Package Infrastructure
- Complete DESCRIPTION file with all six-module dependencies
- NAMESPACE configured for roxygen2 auto-generation
- Package-level documentation in `R/agristat-package.R`
- Testing infrastructure using testthat 3rd edition

#### CI/CD Pipelines (.github/workflows/)
- **R-CMD-check.yaml**: Multi-platform checks (Ubuntu R-release/R-devel,
  macOS, Windows)
- **test-coverage.yaml**: Code coverage tracking with Codecov integration
- **pkgdown.yaml**: Automatic documentation deployment to GitHub Pages
- **lint.yaml**: Code style checking with lintr

#### Documentation
- `_pkgdown.yml`: Bootstrap 5 pkgdown site with all six module reference sections
- `README.Rmd`: Installation instructions, module overview, and quick-start examples
- `CONTRIBUTING.md`: Contribution guidelines
- `CODE_OF_CONDUCT.md`: Contributor Covenant Code of Conduct
- `inst/WORDLIST`: Domain-specific spell-check dictionary

#### Module Stubs (to be fully implemented in Phases 2–7)
- Module 1: Data Ingestion & Pre-Processing (`R/data-ingestion.R`)
- Module 2: Field Layout Engine (`R/field-layouts.R`)
- Module 3: Breeding & Genetics Analytics (`R/genetics.R`)
- Module 4: Mean Separation & Post-Hoc Tests (`R/posthoc.R`)
- Module 5: Phylogenetics & Evolutionary Mapping (`R/phylogenetics.R`)
- Module 6: Thesis Assembly Export Suite (`R/thesis-export.R`)
- Utility functions (`R/utils.R`)
- Dataset documentation (`R/data.R`)
