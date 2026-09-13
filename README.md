# CESR : Backtesting Conditional Portfolio Expected Shortfall with D-Vine Copula Regression

Replication code for the manuscript *"Backtesting Conditional Portfolio Expected
Shortfall with D-Vine Copula Regression"* by Fatima Zahrae Chaayra
(INSEA) and Khalil Said (INSEA / UPPA).

## Overview

The Conditional Expected Shortfall Regression (CESR) test evaluates the
calibration of conditional Expected Shortfall forecasts under external
stress events with strictly positive probability. The test combines a
Fissler–Ziegel identification residual with a Wald-type pooled regression
across a grid of stress-probability levels.

This repository contains the R code used to produce all empirical results,
tables, and figures reported in the manuscript. The application is on
nine large-cap Spanish equities from the MSCI Spain index (2015–2025),
with the S&P 500 and Euro Stoxx 50 as external stress factors.

## Repository structure

```
cesr-portvine/
├── README.md               # this file
├── LICENSE                 # MIT license
├── .gitignore              # ignores large .RData files
├── patch_etape6.R          # CESR test implementation (Z1, Z2, Z3 with event weighting)
├── 01_download_data.R      # downloads Yahoo Finance daily returns
├── 02_estimate_portvine.R  # fits ARMA-GARCH marginals + D-vine copula (rolling windows)
├── 03_simulate_forecasts.R # generates conditional VaR/ES forecasts under stress events
├── 04_run_cesr_tests.R     # applies CESR tests, produces main results
├── 05_robustness.R         # four robustness checks (alpha, sub-periods, T_out)
└── 06_figures.R            # regenerates all figures used in the manuscript
```

## Software requirements

- **R version:** 4.5.1 or higher (tested on R 4.5.1, Windows 10)
- **Required R packages:**
  - `quantmod` (data download from Yahoo Finance)
  - `rugarch` (ARMA-GARCH marginal filtering)
  - `rvinecopulib` (D-vine copula fitting and simulation)
  - `portvine` (rolling window portfolio risk estimation)
  - `ggplot2` (figures)
  - `moments` (descriptive statistics)

Install missing packages with:
```r
install.packages(c("quantmod", "rugarch", "rvinecopulib",
                   "portvine", "ggplot2", "moments"))
```

## How to reproduce the results

Scripts should be executed in numerical order. Each script assumes that
the previous ones have completed and saved their outputs in the working
directory.

```r
setwd("~/reel_data")   # or your local path

source("01_download_data.R")        # ~30 seconds
source("02_estimate_portvine.R")    # ~3 hours (rolling estimation)
source("03_simulate_forecasts.R")   # ~40 minutes
source("04_run_cesr_tests.R")       # ~10 seconds
source("05_robustness.R")           # ~20 minutes
source("06_figures.R")              # ~1 minute
```

Total runtime: approximately 4 hours on a standard laptop, dominated by
the rolling estimation of the D-vine copula in step 02.

## Key parameters (documented in each script)

- **Sample period:** 2015-05-07 to 2025-04-30 (T = 2,446 daily returns)
- **Out-of-sample window:** T_out = 2,196 dates
- **Marginal rolling window:** Ψ = 250 observations, refit every 200
- **D-vine rolling window:** Ψv = 200 observations, refit every 100
- **Tail probability:** α = 0.025
- **Stress-probability grid:** α_I ∈ {0.05, 0.10, 0.25, 0.50}
- **Conditional simulations per date:** S = 1,000
- **Random seed:** 20260817

## Data availability

The daily return data are obtained from Yahoo Finance using the R package
`quantmod`. Ticker symbols used: `IBE.MC`, `SAN.MC`, `BBVA.MC`, `ITX.MC`,
`CLNX.MC`, `AMS.MC`, `TEF.MC`, `REP.MC`, `FER.MC`, `^GSPC`, `^STOXX50E`.
All data are publicly available.

## Correspondence

For questions about the code or the manuscript, please contact the
corresponding author.

## License

MIT License — see `LICENSE` file for details.

## Citation

If you use this code, please cite:

> Chaayra, F. Z. and Said, K. (2026). Backtesting Conditional Portfolio
> Expected Shortfall with D-Vine Copula Regression. *Working paper*.
