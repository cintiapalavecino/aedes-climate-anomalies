# Climate anomalies and *Aedes aegypti* oviposition

R code for analysing how weekly temperature and precipitation anomalies relate to *Aedes aegypti* egg-laying activity in **Santo Tomé (Santa Fe, Argentina)** across three dengue seasons: 2021–22 (non-epidemic), 2022–23 (epidemic) and 2023–24 (epidemic).

> Associated manuscript currently under review.

## Questions

1. How did weekly minimum and maximum temperatures and precipitation during each transmission season (October–May) depart from the 1991–2020 climatology?
2. How did oviposition (weekly eggs per active ovitrap) vary within and between seasons?
3. Do temperature anomalies have delayed, non-linear effects on oviposition?

## Workflow

| Script | Description |
|---|---|
| `R/01_climate_anomalies.R` | Builds weekly means per epidemiological week, computes anomalies against the 1991–2020 baseline, and fits season-specific **GAMs** (Tmin, Tmax) and a **quantile GAM** (median precipitation). Produces Fig. 2. |
| `R/02_oviposition_time_series.R` | Aggregates ovitrap egg counts by season and epidemiological week and plots the seasonal time series (Fig. 3). |
| `R/03_dlnm_climate_oviposition.R` | Fits **distributed lag non-linear models** (`dlnm` cross-bases, lags 0–7 weeks) inside negative binomial models (`glmmTMB`, NB1) with an offset for sampling effort. Includes AIC-based model selection, VIF and `DHARMa` diagnostics (dispersion, temporal autocorrelation), and exposure–response, lag–response and contour plots (Fig. 4). |

## Methods at a glance

- Climate anomaly = weekly observed mean − 1991–2020 mean for the same epidemiological week
- GAMs with season-specific smooths (`mgcv`, REML); quantile GAM for skewed precipitation (`qgam`)
- DLNM: linear exposure for Tmin, natural spline for Tmax, natural spline lag structure (3 df)
- Count model: `total_eggs ~ cb_tmin + cb_tmax + Season + ns(week, 4) + offset(log(active traps))`

## Data

The data are not included in this repository. To run the scripts, place these files in `data/`:

- `Santo_Tome_OVI_weekly.csv`: weekly ovitrap counts (`Date`, `Year`, `EW`, `TOTAL`, …)
- `Santo_Tome_climate_2021-2024.csv`: daily `tmax`, `tmin`, `prcp`
- `Santo_Tome_climate_1991-2020.csv`: historical daily climate (baseline)

Files use `;` as separator and `dd/mm/yyyy` dates. Run scripts from the repository root; outputs are written to `output/`.

## Main packages

`dplyr`, `lubridate`, `ggplot2`, `cowplot`, `mgcv`, `qgam`, `dlnm`, `splines`, `glmmTMB`, `MASS`, `MuMIn`, `DHARMa`, `car`, `writexl`

## Author

Cintia C. Palavecino, CONICET / Universidad Nacional del Litoral, Argentina
