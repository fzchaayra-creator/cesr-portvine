## =============================================================================
##  05_robustness.R
##
##  Four robustness checks on the pooled Z2 statistic, all based on the
##  weighted event-based estimand of patch_etape6.R:
##
##    rob1 : effect of the tail probability level alpha (three values)
##           - requires re-simulating VaR/ES with alpha = 0.01 and 0.05,
##             using the same fitted vines and marginals as step 03
##           - runtime : ~15 minutes
##
##    rob2 : two equal sub-periods of the out-of-sample window
##           - reuses the baseline forecasts, no re-simulation
##           - runtime : < 1 minute
##
##    rob3 : level-specific Z1 tests (already reported in the main results)
##           - reuses the baseline forecasts, no re-simulation
##           - runtime : < 1 minute
##
##    rob4 : effect of the out-of-sample window size T_out
##           - reuses the baseline forecasts, no re-simulation
##           
##
##  Produces cesr_robustness.RData containing rob1..rob4 data frames.
## =============================================================================

setwd("~/reel_data")

if (!file.exists("cesr_forecasts_event_2016_2025.RData"))
  stop("Forecasts not found. Run 03_simulate_forecasts.R first.")
if (!file.exists("cesr_portvine_2016_2025_fits.RData"))
  stop("Portvine fits not found. Run 02_estimate_portvine.R first.")

load("cesr_forecasts_event_2016_2025.RData")
load("cesr_portvine_2016_2025_fits.RData")
source("patch_etape6.R")

library(portvine)
library(rvinecopulib)
library(rugarch)

noms_actifs <- c("IBE","SAN","BBVA","ITX","CLNX","AMS","TEF","REP","FER")
weights     <- setNames(rep(1/9, 9), noms_actifs)

## --- rob1 : effect of alpha in {0.01, 0.025, 0.05} --------------------------
##  Re-uses cached vines/marginals from step 02.
extract_marginals <- function(pv, stress_name) {
  vars <- c(noms_actifs, stress_name); out <- list()
  for (v in vars) out[[v]] <- as.data.frame(
    fitted_marginals(pv)[[v]], which = "density")
  out
}
marg_sp <- extract_marginals(pv_sp, "SP500")
marg_es <- extract_marginals(pv_es, "ESTOXX")
vines_sp <- fitted_vines(pv_sp)
vines_es <- fitted_vines(pv_es)
vine_index_for_t <- function(t) 1L + (t - 1L) %/% refit_vine

n_samples_rob1 <- 250  # lighter than main run for speed
sim_forecast <- function(vines, marg, u_I_series, alpha, alpha_I_grid) {
  T <- length(u_I_series); M <- length(alpha_I_grid); d <- 9
  VaR <- ES <- matrix(NA_real_, T, M)
  for (t in seq_len(T)) {
    vi <- vine_index_for_t(t)
    for (m in seq_len(M)) {
      aI  <- alpha_I_grid[m]
      W   <- matrix(runif(n_samples_rob1 * d), ncol = d)
      U   <- cbind(W, runif(n_samples_rob1, 0, aI))
      Usim<- rvinecopulib::inverse_rosenblatt(U, vines[[vi]])
      r <- matrix(0, n_samples_rob1, d)
      for (j in 1:d) {
        p <- marg[[noms_actifs[j]]]
        z <- rugarch::qdist("sstd", p = Usim[, j], mu = 0, sigma = 1,
                            skew = p$Skew[t], shape = p$Shape[t])
        r[, j] <- p$Mu[t] + p$Sigma[t] * z
      }
      draws <- drop(r %*% weights) / 100
      q <- quantile(draws, probs = alpha, type = 1, names = FALSE)
      VaR[t, m] <- q; ES[t, m] <- mean(draws[draws <= q])
    }
    if (t %% 500 == 0) cat(sprintf("    t = %d / %d\n", t, T))
  }
  list(VaR = VaR, ES = ES)
}

set.seed(20260817)
alphas <- c(0.01, 0.025, 0.05)
rob1 <- data.frame(alpha = alphas, p_SP500 = NA, p_ESTOXX = NA)
for (i in seq_along(alphas)) {
  a <- alphas[i]; cat(sprintf("\n[rob1] alpha = %.3f\n", a))
  cat("  -- S&P 500 --\n")
  fc_sp <- sim_forecast(vines_sp, marg_sp, u_I_sp, a, alpha_I_grid)
  W_sp  <- sapply(1:4, function(j) as.numeric(u_I_sp <= alpha_I_grid[j]))
  XI_sp <- sapply(1:4, function(j) compute_xi(fc_sp$VaR[, j], fc_sp$ES[, j],
                                              port_oos, alpha = a))
  rob1$p_SP500[i]  <- cesr_pooled(XI_sp, 100*alpha_I_grid, W = W_sp,
                                   spec = "Z2")$p_value
  cat("  -- Euro Stoxx 50 --\n")
  fc_es <- sim_forecast(vines_es, marg_es, u_I_es, a, alpha_I_grid)
  W_es  <- sapply(1:4, function(j) as.numeric(u_I_es <= alpha_I_grid[j]))
  XI_es <- sapply(1:4, function(j) compute_xi(fc_es$VaR[, j], fc_es$ES[, j],
                                              port_oos, alpha = a))
  rob1$p_ESTOXX[i] <- cesr_pooled(XI_es, 100*alpha_I_grid, W = W_es,
                                   spec = "Z2")$p_value
}
print(rob1)

## --- Helper for rob2 and rob4 : test on a date sub-set ---------------------
test_subperiod <- function(idx, fc, uI, y, alpha_I_grid) {
  W  <- sapply(seq_along(alpha_I_grid), function(j)
    as.numeric(uI[idx] <= alpha_I_grid[j]))
  XI <- sapply(seq_along(alpha_I_grid), function(j)
    compute_xi(fc$VaR[idx, j], fc$ES[idx, j], y[idx], alpha = 0.025))
  cesr_pooled(XI, 100*alpha_I_grid, W = W, spec = "Z2")$p_value
}

## --- rob2 : sub-periods (two equal halves) ---------------------------------
mid <- T_out %/% 2
rob2 <- data.frame(
  factor    = c("SP500","SP500","EuroStoxx","EuroStoxx"),
  subperiod = c("1st half","2nd half","1st half","2nd half"),
  n         = c(mid, T_out - mid, mid, T_out - mid),
  p_Z2      = c(
    test_subperiod(1:mid,           fc_sp500,  u_I_sp, port_oos, alpha_I_grid),
    test_subperiod((mid+1):T_out,   fc_sp500,  u_I_sp, port_oos, alpha_I_grid),
    test_subperiod(1:mid,           fc_estoxx, u_I_es, port_oos, alpha_I_grid),
    test_subperiod((mid+1):T_out,   fc_estoxx, u_I_es, port_oos, alpha_I_grid)
  )
)
print(rob2)

## --- rob3 : level-specific Z1 tests ----------------------------------------
W_sp <- sapply(1:4, function(j) as.numeric(u_I_sp <= alpha_I_grid[j]))
W_es <- sapply(1:4, function(j) as.numeric(u_I_es <= alpha_I_grid[j]))
XI_sp <- sapply(1:4, function(j) compute_xi(fc_sp500$VaR[, j],
                                            fc_sp500$ES[, j], port_oos, 0.025))
XI_es <- sapply(1:4, function(j) compute_xi(fc_estoxx$VaR[, j],
                                            fc_estoxx$ES[, j], port_oos, 0.025))
rob3 <- data.frame(
  factor  = rep(c("SP500","EuroStoxx"), each = 4),
  alpha_I = rep(alpha_I_grid, 2),
  p_Z1    = c(sapply(1:4, function(j) cesr_Z1(XI_sp[, j], w = W_sp[, j])$p_value),
              sapply(1:4, function(j) cesr_Z1(XI_es[, j], w = W_es[, j])$p_value))
)
print(rob3)

## --- rob4 : effect of T_out ------------------------------------------------
T_grid <- c(500, 1000, 1500, 2000)
rob4 <- data.frame(Tout = T_grid, p_SP500 = NA, p_ESTOXX = NA)
for (i in seq_along(T_grid)) {
  T_i <- T_grid[i]
  rob4$p_SP500[i]  <- test_subperiod(1:T_i, fc_sp500,  u_I_sp, port_oos, alpha_I_grid)
  rob4$p_ESTOXX[i] <- test_subperiod(1:T_i, fc_estoxx, u_I_es, port_oos, alpha_I_grid)
}
print(rob4)

save(rob1, rob2, rob3, rob4, file = "cesr_robustness.RData")
cat("\n=== DONE. Robustness results saved to cesr_robustness.RData ===\n")
