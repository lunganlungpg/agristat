# Contributing to agristat

Thank you for your interest in contributing to **agristat**!

## Code of Conduct

Please review our [Code of Conduct](CODE_OF_CONDUCT.md) before contributing.

## Reporting Bugs

1. Check [existing issues](https://github.com/lunganlungpg/agristat/issues).
2. Open a new issue with a minimal reproducible example (`reprex::reprex()`).
3. Describe expected vs. actual behaviour.

## Submitting Pull Requests

1. Fork the repository and create a branch from `main`.
2. Follow the [tidyverse style guide](https://style.tidyverse.org/).
3. Write tests for new functions using `testthat`.
4. Add roxygen2 documentation for all exported functions.
5. Run `devtools::check()` — no errors or warnings.
6. Open a pull request targeting `main`.

## Development Setup

```r
install.packages("devtools")
devtools::install_dev_deps()
devtools::check()
devtools::test()
devtools::document()
```
