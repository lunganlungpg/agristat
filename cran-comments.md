## R CMD CHECK results

0 errors | 0 warnings | 0 notes

## Test environments

- Local: Ubuntu 24.04, R 4.3.3
- GitHub Actions: ubuntu-latest (R release, R devel), macos-latest (R release), windows-latest (R release)

## Downstream dependencies

There are no existing reverse dependencies on CRAN (first submission).

## Notes

- This is the first submission of the agristat package (v1.0.0).
- The package implements a comprehensive biostatistical suite for
  agricultural and genetic research, covering six analytical modules:
  data ingestion, field layout design, breeding/genetics analytics,
  mean separation/post-hoc tests, phylogenetics, and thesis export.
- All examples in the documentation are self-contained and do not
  require any external data files beyond the built-in datasets
  included in the package.
- The `sf` package is listed in Imports because spatial field layouts
  are a core feature. If `sf` is unavailable, relevant functions
  provide informative error messages.
