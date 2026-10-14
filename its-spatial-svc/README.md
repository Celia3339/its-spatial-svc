# Spatiotemporal Bayesian ITS Model with Spatially Varying Coefficients

**Article:** *A spatiotemporal Bayesian interrupted time series model with spatially varying coefficients: application to fentanyl prescription regulation in the Valencian Community*  
**Journal:** Statistics in Medicine (under review)  
**Authors:** Celia Robles · Adina Iftimi  
**Affiliation:** FISABIO / Universitat de València

---

## Overview

This repository contains the full R code for a Bayesian spatiotemporal interrupted time series (ITS) model fitted with R-INLA. The model combines:

- A **Leroux CAR** prior for the baseline spatial random effect
- An **AR(1)** process for the temporal random effect
- **Spatially varying coefficients** (SVCs, Besag model) for the immediate level shifts and the post-intervention trajectory spline bases at each of the three interventions

The analysis covers 241 basic health areas of the Valencian Community over 209 weeks (2015–2018), with opioid prescription counts as the outcome.

---

## Repository structure

```
its-spatial-svc/
├── R/
│   ├── 00_setup.R               # packages, paths, study constants
│   ├── functions_model.R        # data preparation, priors, formulas, two-stage INLA fit
│   ├── functions_trajectories.R # posterior trajectory summaries
│   └── functions_simulation.R   # data-generating mechanism and ADEMP metrics
├── scripts/
│   ├── 01_fit_models.R          # fit Poisson and NB models; model comparison table
│   ├── 02_trajectories.R        # population-average and area-level trajectories
│   ├── 03_simulation.R          # simulation study (dependence + effect-size scenarios)
│   ├── 04_diagnostics.R         # mid-PIT, Pearson residual ACF, Moran's I
│   └── 05_robustness_checks.R   # overdispersion, Poisson vs NB, temporal structure
├── data/                        # ← place your input files here (not distributed)
├── output/                      # fitted models, summaries and figures (git-ignored)
├── its-spatial-svc.Rproj        # open this in RStudio
├── CITATION.cff
└── LICENSE
```

---

## Data

The health data cannot be distributed. Place the following files in `data/` before running any script:

| File | Description |
|------|-------------|
| `analysis_data.rds` | Long-format data frame, 241 areas × 209 weeks. Required columns: `id_zona` (integer, 1–241), `semana` (integer, 1–209), `suma_zona` (integer count), `poblacion` (population offset), `int1 / int2 / int3` (step indicators), `postint1 / postint2 / postint3` (elapsed-time counters), `mayor65_porc`, `homb_porc`, `tipo_suelo` |
| `adjacency.RData` | Contains object `adj_dep`: a 241 × 241 adjacency matrix of the areas |

---

## Reproducibility

Open `its-spatial-svc.Rproj` in RStudio, then run the scripts in order:

```r
source("scripts/01_fit_models.R")       # ~2 h
source("scripts/02_trajectories.R")     # ~20 min
source("scripts/03_simulation.R")       # ~24 h (200 × 3 scenarios)
source("scripts/04_diagnostics.R")      # <5 min
source("scripts/05_robustness_checks.R")# ~4 h
```

Set `TEST_MODE <- TRUE` in `scripts/03_simulation.R` to run only 5 replicates and verify the pipeline.

---

## Software

| Package | Version tested |
|---------|---------------|
| R | 4.3.x |
| INLA | 24.x |
| Matrix | ≥ 1.6 |
| splines | (base) |
| dplyr / tidyr | ≥ 1.3 / 1.3 |
| spdep | ≥ 1.3 |
| igraph | ≥ 2.0 |

Install INLA:
```r
install.packages("INLA", repos = c(getOption("repos"),
  INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

---

## Citation

If you use this code, please cite:

> Robles C, Iftimi A. A spatiotemporal Bayesian interrupted time series model with spatially varying coefficients. *Statistics in Medicine* (under review).

A Zenodo DOI will be added upon publication.

---

## License

GPL-3 © 2026 Celia Robles, Adina Iftimi — see `LICENSE`.
