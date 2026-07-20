# Length at Maturity — DATRAS

Estimates of length at maturity (L50, L25, L75) for Northeast Atlantic fish and elasmobranch species, derived from the ICES DATRAS trawl survey database. Models are fitted at multiple spatial and temporal scales, quality-filtered, and checked against FishBase reference values.

---

## Overview

Maturity status is recorded during bottom-trawl surveys but the data are rarely used at scale to produce comparable L50 estimates across species and regions. This repository builds a reproducible pipeline that:

1. Downloads and harmonises raw maturity observations from DATRAS across all available surveys
2. Assigns each observation to its ICES statistical area and filters to the spawning / pre-spawning period using two independent sources (GO-FISH and WKMAT07)
3. Decodes heterogeneous maturity scales (SMSF, Roman, numeric, R1/R2, M6, I/M) into a common binary mature/immature variable
4. Fits logistic regression models with delta-method confidence intervals to estimate L25, L50 and L75 at global, stock, area, and 5-year period levels
5. Quality-filters the estimates, reports per-model fit (AUC, Tjur R²), checks them against FishBase, and exports a tidy CSV for archiving

---

## Repository structure

```
length_at_maturity/
│
├── R/
│   ├── 00_download.R              # download DATRAS surveys, ICES shapefile, GO-FISH
│   ├── 01_get_species_info.R      # fetch species metadata from FishBase / WoRMS
│   ├── 02_create_spawning_data.R  # build spawning-season lookup (GO-FISH + WKMAT07)
│   ├── 03_assign_maturity.R       # decode maturity scales, spatial join, season filter
│   ├── 04_fit_maturity_models.R   # fit logistic GLMs, delta-method L25/L50/L75
│   ├── 05_clean_estimates.R       # quality-filter estimates, export CSV
│   ├── 06_plot_workflow.R         # workflow diagram
│   ├── 07_plot_map.R              # sampling coverage map
│   ├── 08_plot_l50.R              # main figure (dot plot, stock map, time series, ogive)
│   ├── 09_plot_performance.R      # AUC / Tjur R² by aggregation level
│   └── 10_plot_fishbase.R         # global L50 vs FishBase median
│
├── tests/
│   └── bits_sensitivity.R         # Baltic Roman stage-II sensitivity analysis
│
├── data/
│   ├── raw/                       # raw DATRAS downloads (not tracked)
│   ├── metadata/                  # ICES areas, GO-FISH shapefiles, WKMAT07, species info
│   ├── intermediate/              # cleaned maturity data, raw L50 estimates
│   └── final/                     # validated estimates (l50_clean.rds, l50_estimates.csv)
│
├── outputs/
│   ├── main/                      # primary figures
│   └── supp/                      # supplementary figures
│
└── length_at_maturity.Rproj
```

---

## Pipeline

Each script is numbered and intended to be run in order. Scripts read from and write to `data/` subfolders via `here::here()`.

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 0 | `00_download.R` | — | `data/raw/raw_datras.rds`, ICES shapefile, GO-FISH |
| 1 | `01_get_species_info.R` | raw DATRAS | `data/metadata/species_info.rds` |
| 2 | `02_create_spawning_data.R` | species info, GO-FISH, WKMAT07 | `data/metadata/spawning_lookup.rds` |
| 3 | `03_assign_maturity.R` | raw DATRAS, spawning lookup | `data/intermediate/maturity_clean.rds` |
| 4 | `04_fit_maturity_models.R` | cleaned maturity data | `data/intermediate/l50_raw.rds` |
| 5 | `05_clean_estimates.R` | raw L50 estimates | `data/final/l50_estimates.csv` |
| 6–10 | `06`–`10_plot_*.R` | cleaned data, final estimates | figures in `outputs/` |

---

## Spawning season filter

To ensure observations reflect reproductive condition, each record is matched to a species × ICES-area spawning window using:

- **GO-FISH** — spatially explicit monthly spawning habitat maps
- **WKMAT07** — ICES expert-elicited spawning seasons by area

Areas are resolved through the full ICES hierarchy (e.g. `27.4.b` falls back to `27.4`, then `27`). The active window is extended by 3 months before the first spawning month to capture pre-spawning fish. Records with no spawning information are retained rather than dropped.

---

## Maturity scales

DATRAS surveys use several incompatible maturity coding schemes (SMSF, Roman, numeric, R1/R2, M6, I/M). The pipeline detects the scheme in use for each country × survey × year and converts it to a binary immature/mature variable. "Roman" and "Numeric" are not distinct biological scales — they refer to how stages are encoded in DATRAS (Roman numerals vs. integers) and may correspond to different underlying scales depending on the survey. The full per-survey coding is written to `data/metadata/maturity_scale_summary.{csv,tex}`.

For the Baltic surveys, the ambiguous Roman stage II is dropped by default, with documented exceptions (Russia/Lithuania, the SE-SOUND survey, and herring in the Polish and Swedish BITS). The sensitivity of L50 to this decision is explored in `tests/bits_sensitivity.R`.

---

## Models

A binomial GLM (`mature ~ lngt_cm`, logit link) is fitted for each species × grouping combination. L25, L50 and L75 are derived analytically from the coefficients, with delta-method standard errors and symmetric 95% confidence intervals. Model fit is reported via AUC and Tjur R².

Grouping levels:

- **Global** — all areas and years combined
- **Stock** — ICES stock unit (`.27.` notation)
- **Area** — individual ICES statistical area
- **Period** — 5-year bins (2000–2004, …, 2020–2024)

Each level is run for combined sexes and male/female separately where sample sizes allow (≥ 15 immature and ≥ 15 mature).

---

## Key dependencies

```r
# CRAN
tidyverse, here, sf, janitor, jsonlite, rfishbase, worrms, taxize, readxl,
pROC, performance, knitr, kableExtra, rnaturalearth, patchwork, ggrepel,
paletteer, geomtextpath, ggstats, DiagrammeR, DiagrammeRsvg, rsvg

# GitHub
remotes::install_github("tokami/DATRASextra")
```

---

## Output

`data/final/l50_estimates.csv` — one row per species × sex × spatial scale × time period, with columns:

`aphia_id`, `species`, `sex`, `spatial_scale`, `ices_stock`, `ices_area`, `time_period`, `sample_size`, `l25_cm`, `l25_se_cm`, `l25_ci_lower_cm`, `l25_ci_upper_cm`, `l50_cm`, `l50_se_cm`, `l50_ci_lower_cm`, `l50_ci_upper_cm`, `l75_cm`, `l75_se_cm`, `l75_ci_lower_cm`, `l75_ci_upper_cm`, `model_auc`, `model_tjur_r2`
