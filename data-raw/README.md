# data-raw/

This directory contains scripts used to generate and clean the datasets
bundled with the `agristat` package (stored as `.rda` files in `data/`).

Each dataset has a corresponding `<dataset-name>.R` script here that
documents the data provenance and transformation steps.

Planned scripts:
- `rice_yield.R`     — Multi-environment rice yield trial data
- `maize_trials.R`   — CIMMYT-style maize performance trials
- `wheat_pedigree.R` — Example wheat pedigree for genetic analyses
- `soil_profile.R`   — Soil characteristics linked to field layouts
