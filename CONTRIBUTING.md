# Contributing to agristat

Thank you for your interest in contributing to **agristat**! This document
provides guidelines for contributing code, documentation, and bug reports.

## Code of Conduct

Please review our
[Code of Conduct](CODE_OF_CONDUCT.md) before contributing.
All participants are expected to uphold these standards.

## How to Contribute

### Reporting Bugs

1. Check the [existing issues](https://github.com/lunganlungpg/agristat/issues)
   to avoid duplicates.
2. Open a new issue using the **Bug Report** template.
3. Include a minimal reproducible example (`reprex::reprex()`).
4. Describe expected vs. actual behaviour.

### Suggesting Features

1. Open a new issue using the **Feature Request** template.
2. Describe the proposed feature and its use case in agricultural/genetic research.
3. Reference relevant literature or methodology where appropriate.

### Submitting Pull Requests

1. Fork the repository and create a branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```
2. Follow the coding style described below.
3. Write tests for new functions using `testthat`.
4. Add roxygen2 documentation for all exported functions.
5. Run `devtools::check()` and ensure there are no errors or warnings.
6. Open a pull request targeting `main`.

## Coding Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/).
- Use `snake_case` for function and variable names.
- Prefer explicit namespace calls (`pkg::fun()`) for non-base functions in
  internal helpers.
- Keep lines to ≤ 80 characters.
- Run `lintr::lint_package()` before submitting a PR.

## Documentation

- All exported functions **must** have roxygen2 documentation.
- Include `@examples` blocks that are self-contained and run without errors.
- Vignettes should follow the module structure (one per module).

## Testing

- Aim for ≥ 80 % code coverage.
- Tests live in `tests/testthat/` and are named `test-<module>.R`.
- Use `testthat::test_that()` with descriptive test names.
- Run `devtools::test()` locally before pushing.

## Development Setup

```r
# Install development dependencies
install.packages("devtools")
devtools::install_dev_deps()

# Run checks
devtools::check()

# Run tests
devtools::test()

# Generate documentation
devtools::document()

# Build pkgdown site locally
pkgdown::build_site()
```

## Questions?

Open a
[GitHub Discussion](https://github.com/lunganlungpg/agristat/discussions)
or email the maintainer at <lunganlungpg@example.com>.
