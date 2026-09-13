## =============================================================================
##  02_estimate_portvine.R
##
##  Fits rolling-window ARMA-GARCH marginals and D-vine copulas via the
##  portvine package, separately for each of the two external stress factors
##  (S&P 500 and Euro Stoxx 50).
##
##  Produces cesr_portvine_2016_2025_fits.RData containing:
##    - pv_sp : portvine object for the S&P 500 stress factor
##    - pv_es : portvine object for the Euro Stoxx 50 stress factor
##
##  This is the computationally heaviest step of the pipeline.
##
##  Note on marginal specification : the theoretical model in the manuscript
##  is ARMA(1,1)-GARCH(1,1); the empirical implementation uses a constant
##  conditional mean because the ARMA(1,1) fails to converge on several
##  rolling windows during volatile periods (near-cancelling AR/MA roots).
##  For daily returns the conditional mean is well-approximated by a
##  constant, and all dynamics of interest live in the variance equation.
## =============================================================================

library(portvine)
library(rugarch)

set.seed(20260728)

## --- Load data --------------------------------------------------------------
if (!file.exists("cesr_data_2015_2025.RData"))
  stop("Data file not found. Run 01_download_data.R first.")
load("cesr_data_2015_2025.RData")

## --- Parameters (matched to manuscript) -------------------------------------
noms_actifs <- c("IBE","SAN","BBVA","ITX","CLNX","AMS","TEF","REP","FER")
weights     <- setNames(rep(1/9, 9), noms_actifs)

alpha_ES   <- 0.025    # tail probability
T_in_marg  <- 250      # marginal training window size (Psi)
refit_marg <- 200      # marginal refit frequency
T_in_vine  <- 200      # D-vine training window size (Psi_v)
refit_vine <- 100      # D-vine refit frequency

## portvine constraints : refit_marg divisible by refit_vine (200/100 OK),
##                        T_in_vine <= T_in_marg (200 <= 250 OK).

fichier_portvine <- "cesr_portvine_2016_2025_fits.RData"

## --- Fitting function -------------------------------------------------------
##  Note : rugarch is numerically unstable on decimal returns (omega ~ 4e-6
##  makes the Hessian non-invertible). We fit in PERCENT (x100). The PIT is
##  scale-invariant; mu and sigma extracted later are in %, and the final
##  simulation step in 03 converts back to decimals.
run_fit <- function(stress_name) {
  dat <- as.data.frame(coredata(100 * ret[, c(noms_actifs, stress_name)]))
  colnames(dat) <- c(noms_actifs, stress_name)
  w   <- c(weights, setNames(0, stress_name))

  spec_robuste <- default_garch_spec(ar = 0, ma = 0)  # constant mean

  estimate_risk_roll(
    data              = dat,
    weights           = w,
    marginal_settings = marginal_settings(train_size   = T_in_marg,
                                          refit_size   = refit_marg,
                                          default_spec = spec_robuste),
    vine_settings     = vine_settings(train_size = T_in_vine,
                                      refit_size = refit_vine,
                                      family_set = "parametric",
                                      vine_type  = "dvine"),
    alpha             = alpha_ES,
    risk_measures     = c("VaR"),
    n_samples         = 100,        # minimal: these forecasts are discarded
    cond_vars         = stress_name,
    cond_u            = 0.05,       # single level: only fitted objects matter
    trace             = TRUE)
}

## --- Fit both factors -------------------------------------------------------
if (file.exists(fichier_portvine)) {
  cat("Fit file already exists:", fichier_portvine, "\n")
  cat("Delete it and re-run this script to refit (~3h).\n")
} else {
  cat("\n=== Fitting portvine for S&P 500 (long) ===\n")
  pv_sp  <- run_fit("SP500")
  cat("\n=== Fitting portvine for Euro Stoxx 50 (long) ===\n")
  pv_es  <- run_fit("ESTOXX")
  save(pv_sp, pv_es,
       T_in_marg, refit_marg, T_in_vine, refit_vine, alpha_ES,
       file = fichier_portvine)
  cat("\nSaved:", fichier_portvine, "\n")
}
